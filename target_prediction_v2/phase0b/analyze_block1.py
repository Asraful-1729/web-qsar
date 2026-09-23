"""
Phase 0b, Block 1 + reliability/popularity/threshold-sweep analysis
(P0b-3, P0b-4, P0b-5, P0b-6, P0b-13, P0b-14, P0b-18).

Consumes capture.py's per-query JSONL (results/capture_{A,B,D}.jsonl) and
capture_query_potency.py's data/query_own_potency.json -- no re-querying
the index; every metric here is a fast in-memory recomputation over saved
per-query neighbour/result data.

MACRO vs POOLED (P0b-3, rev5 Section 1.1's "averaging trap"): computed
explicitly, side by side, for every precision/recall table in this script.
Macro = mean of per-query precision (over queries with >=1 predicted
target; precision is undefined, not zero, for a query with no predictions
-- excluded from the mean, not counted as 0) and per-query recall (every
query here has >=1 true annotated target by construction, so recall is
always defined). Pooled = pooled TP / pooled (TP+FP) and TP / (TP+FN)
across all queries combined. This script NEVER reports one without the
other where both apply -- rev5's E1/E2 errata exist precisely because a
single number without its convention stated is not interpretable.

Target universe size for MCC: fixed at 4,658 (v1's frozen global count),
matching Phase 0's convention -- leakage removal changes this by at most a
handful of compounds' worth of targets per query, immaterial to MCC's true-
negative-dominated denominator. Disclosed, not hidden.

Usage: python3 analyze_block1.py
"""
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
DATA_DIR = os.path.join(HERE, "data")

N_TARGETS_UNIVERSE = 4658
THRESHOLDS = [0.2, 0.3, 0.4, 0.5, 0.6]
WEAK_POTENCY_10uM = 5.0
WEAK_POTENCY_1uM = 6.0


def load_capture(setname):
    path = os.path.join(RESULTS_DIR, f"capture_{setname}.jsonl")
    with open(path) as f:
        return [json.loads(line) for line in f]


def load_query_own_potency():
    with open(os.path.join(DATA_DIR, "query_own_potency.json")) as f:
        return json.load(f)


# ---------- P0b-3: macro vs pooled, best_similarity ranking, threshold sweep ----------

def precision_recall_mcc_per_query(predicted_targets, true_targets, n_universe=N_TARGETS_UNIVERSE):
    predicted = set(predicted_targets)
    true = set(true_targets)
    tp = len(predicted & true)
    fp = len(predicted - true)
    fn = len(true - predicted)
    tn = n_universe - tp - fp - fn
    precision = (tp / len(predicted)) if predicted else None
    recall = (tp / len(true)) if true else None
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = ((tp * tn - fp * fn) / denom) if denom > 0 else None
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "precision": precision, "recall": recall, "mcc": mcc,
            "n_predicted": len(predicted), "n_true": len(true)}


def macro_and_pooled(per_query_stats):
    precisions = [s["precision"] for s in per_query_stats if s["precision"] is not None]
    recalls = [s["recall"] for s in per_query_stats if s["recall"] is not None]
    mccs = [s["mcc"] for s in per_query_stats if s["mcc"] is not None]
    tp_sum = sum(s["tp"] for s in per_query_stats)
    fp_sum = sum(s["fp"] for s in per_query_stats)
    fn_sum = sum(s["fn"] for s in per_query_stats)
    n_pred_sum = sum(s["n_predicted"] for s in per_query_stats)
    n_true_sum = sum(s["n_true"] for s in per_query_stats)
    return {
        "n_queries": len(per_query_stats),
        "macro_precision": round(statistics.mean(precisions), 4) if precisions else None,
        "macro_precision_n": len(precisions),
        "macro_recall": round(statistics.mean(recalls), 4) if recalls else None,
        "macro_recall_n": len(recalls),
        "macro_mcc": round(statistics.mean(mccs), 4) if mccs else None,
        "pooled_precision": round(tp_sum / n_pred_sum, 4) if n_pred_sum else None,
        "pooled_recall": round(tp_sum / n_true_sum, 4) if n_true_sum else None,
        "relative_diff_precision_pct": (round(100 * (statistics.mean(precisions) - tp_sum / n_pred_sum) / (tp_sum / n_pred_sum), 1)
                                         if precisions and n_pred_sum else None),
        "relative_diff_recall_pct": (round(100 * (statistics.mean(recalls) - tp_sum / n_true_sum) / (tp_sum / n_true_sum), 1)
                                      if recalls and n_true_sum else None),
        "mean_predicted_targets_per_query": round(n_pred_sum / len(per_query_stats), 2) if per_query_stats else None,
    }


def p0b3_macro_vs_pooled(records, tier="annotated_targets"):
    out = {}
    for th in THRESHOLDS:
        per_query = []
        for r in records:
            predicted = [t for t, sim in r["bestsim_results"] if sim >= th]
            per_query.append(precision_recall_mcc_per_query(predicted, r[tier]))
        out[str(th)] = macro_and_pooled(per_query)
    return out


# ---------- P0b-4 / P0b-18: unweighted 10-NN baseline + L-score reliability ----------

def p0b4_and_lscore(records, tier="annotated_targets"):
    any_ranks, primary_ranks = [], []
    l_buckets = {}  # l (0-10) -> list of (is_true_positive) for that l's TOP prediction? Actually rev5's L
    # table is per-PREDICTION, not per-query: for every (query, predicted target) pair with vote count l,
    # is that predicted target actually true? mean precision per l level.
    l_records = {l: {"n_predictions": 0, "n_correct": 0, "queries": set()} for l in range(0, 11)}

    for qi, r in enumerate(records):
        scores = S.score_query(r["neighbours"], k=10, alpha=0)  # alpha=0 -> plain vote count = "l"
        ranked = S.rank_from_scores(scores)
        annotated = set(r["annotated_targets"])
        primary = set(r["primary_targets"])
        any_ranks.append(S.best_rank(ranked, annotated))
        if primary:
            primary_ranks.append(S.best_rank(ranked, primary))

        for tcid, l in scores.items():
            l_int = int(round(l))
            if l_int < 0 or l_int > 10:
                continue
            l_records[l_int]["n_predictions"] += 1
            l_records[l_int]["queries"].add(qi)
            if tcid in annotated:
                l_records[l_int]["n_correct"] += 1

    def topk_mrr(ranks, ks=(1, 5, 10)):
        n = len(ranks)
        if n == 0:
            return {"n": 0}
        out = {"n": n}
        for k in ks:
            out[f"top_{k}"] = round(sum(1 for x in ranks if x is not None and x <= k) / n, 4)
        out["mrr"] = round(sum((1.0 / x) if x is not None else 0.0 for x in ranks) / n, 4)
        return out

    lscore_table = {}
    for l in range(1, 11):  # l=0 (no votes at all) is not a "prediction", skip
        rec = l_records[l]
        if rec["n_predictions"] == 0:
            continue
        lscore_table[f"L={l/10:.1f}"] = {
            "l": l,
            "n_queries_with_this_l": len(rec["queries"]),
            "n_predictions": rec["n_predictions"],
            "mean_precision": round(rec["n_correct"] / rec["n_predictions"], 4),
        }

    return {
        "baseline_10nn_any_annotated": topk_mrr(any_ranks),
        "baseline_10nn_primary": topk_mrr(primary_ranks),
        "lscore_reliability_table": lscore_table,
    }


# ---------- P0b-5: index-filter vs ground-truth-filter potency ----------

def p0b5_potency_scope(records, query_potency, tier="annotated_targets"):
    out = {"index_filter": {}, "ground_truth_filter": {}}

    for label, floor in [("baseline_no_filter", None), ("10uM_5.0", 5.0), ("1uM_6.0", 6.0)]:
        any_ranks = []
        for r in records:
            pot_fn = None if floor is None else (lambda p, fl=floor: 1.0 if (p is None or p >= fl) else 0.0)
            scores = S.score_query(r["neighbours"], k=25, alpha=1.0, potency_fn=pot_fn)
            ranked = S.rank_from_scores(scores)
            any_ranks.append(S.best_rank(ranked, r[tier]))
        n = len(any_ranks)
        out["index_filter"][label] = {
            "top_1": round(sum(1 for x in any_ranks if x is not None and x <= 1) / n, 4),
            "top_10": round(sum(1 for x in any_ranks if x is not None and x <= 10) / n, 4),
            "mrr": round(sum((1.0 / x) if x is not None else 0.0 for x in any_ranks) / n, 4),
        }

    for label, floor in [("baseline_no_filter", None), ("10uM_5.0", 5.0), ("1uM_6.0", 6.0)]:
        any_ranks = []
        for r in records:
            own_pot = query_potency.get(r["smiles"], {})
            if floor is None:
                true_targets = r[tier]
            else:
                true_targets = [t for t in r[tier] if (own_pot.get(t) is None or own_pot.get(t) >= floor)]
            scores = S.score_query(r["neighbours"], k=25, alpha=1.0)
            ranked = S.rank_from_scores(scores)
            any_ranks.append(S.best_rank(ranked, true_targets) if true_targets else None)
        valid = [x for x in any_ranks]
        n_with_gt = sum(1 for r, own in [(r, query_potency.get(r["smiles"], {})) for r in records]
                         if floor is None or any((own.get(t) is None or own.get(t) >= floor) for t in r[tier]))
        n = len(valid)
        out["ground_truth_filter"][label] = {
            "n_queries_with_qualifying_true_target": n_with_gt,
            "top_1": round(sum(1 for x in valid if x is not None and x <= 1) / n, 4) if n else None,
            "top_10": round(sum(1 for x in valid if x is not None and x <= 10) / n, 4) if n else None,
        }
    return out


# ---------- P0b-6: weak-potency fraction among VOTING neighbours ----------

def p0b6_weak_voting_neighbours(records, k=10):
    n_pairs = 0
    n_weak_10uM = 0
    n_weak_1uM = 0
    n_no_pchembl = 0
    for r in records:
        for smi, sim, targets in r["neighbours"][:k]:
            for tcid, pchembl in targets:
                n_pairs += 1
                if pchembl is None:
                    n_no_pchembl += 1
                    continue
                if pchembl < WEAK_POTENCY_10uM:
                    n_weak_10uM += 1
                if pchembl < WEAK_POTENCY_1uM:
                    n_weak_1uM += 1
    return {
        "n_voting_compound_target_pairs": n_pairs,
        "n_no_pchembl": n_no_pchembl,
        "pct_weaker_than_10uM": round(100 * n_weak_10uM / n_pairs, 2) if n_pairs else None,
        "pct_weaker_than_1uM": round(100 * n_weak_1uM / n_pairs, 2) if n_pairs else None,
    }


# ---------- P0b-13: popularity baseline ----------

def p0b13_popularity(records, popularity, tier="annotated_targets"):
    ranked_targets = sorted(popularity.items(), key=lambda kv: -kv[1])
    rank_map = {t: i + 1 for i, (t, _) in enumerate(ranked_targets)}
    ranks = [S.best_rank([(t, 0) for t, _ in ranked_targets], r[tier]) for r in records]
    # best_rank expects a ranked list of (target, score); reuse rank_map directly instead, faster:
    ranks = []
    for r in records:
        rr = [rank_map[t] for t in r[tier] if t in rank_map]
        ranks.append(min(rr) if rr else None)
    n = len(ranks)
    return {
        "n_queries": n,
        "top_1": round(sum(1 for x in ranks if x is not None and x <= 1) / n, 4),
        "top_10": round(sum(1 for x in ranks if x is not None and x <= 10) / n, 4),
        "mrr": round(sum((1.0 / x) if x is not None else 0.0 for x in ranks) / n, 4),
    }


# ---------- P0b-14: similarity-threshold sweep, rev5 Section 1.6 table shape ----------

def p0b14_threshold_sweep(records, tier="annotated_targets"):
    out = {}
    n = len(records)
    for th in THRESHOLDS + ["top10_no_cutoff"]:
        per_query = []
        n_no_prediction = 0
        for r in records:
            if th == "top10_no_cutoff":
                predicted = [t for t, sim in r["bestsim_results"][:10]]
            else:
                predicted = [t for t, sim in r["bestsim_results"] if sim >= th]
            if not predicted:
                n_no_prediction += 1
            per_query.append(precision_recall_mcc_per_query(predicted, r[tier]))
        agg = macro_and_pooled(per_query)
        mccs = [s["mcc"] for s in per_query if s["mcc"] is not None]
        out[str(th)] = {
            "n_queries_no_prediction": n_no_prediction,
            "pct_queries_no_prediction": round(100 * n_no_prediction / n, 1),
            "mean_predicted_targets_per_query": agg["mean_predicted_targets_per_query"],
            "macro_mcc": round(statistics.mean(mccs), 4) if mccs else None,
            "macro_precision": agg["macro_precision"],
            "macro_recall": agg["macro_recall"],
        }
    return out


def main():
    query_potency = load_query_own_potency()
    report = {}
    for setname in ["A", "B", "D"]:
        records = load_capture(setname)
        print(f"Set {setname}: {len(records)} captured queries", flush=True)
        report[setname] = {
            "p0b3_macro_vs_pooled_any_annotated": p0b3_macro_vs_pooled(records, "annotated_targets"),
            "p0b3_macro_vs_pooled_primary": p0b3_macro_vs_pooled(records, "primary_targets"),
            "p0b4_and_p0b18": p0b4_and_lscore(records, "annotated_targets"),
            "p0b5_potency_scope": p0b5_potency_scope(records, query_potency, "annotated_targets"),
            "p0b6_weak_voting_neighbours": p0b6_weak_voting_neighbours(records),
            "p0b14_threshold_sweep": p0b14_threshold_sweep(records, "annotated_targets"),
        }

    # popularity baseline needs the FULL index's global target frequency -- reuse phase0's approach
    sys.path.insert(0, os.path.join(HERE, "..", "phase0"))
    import os as _os
    _os.environ.setdefault("TARGET_FISHING_INDEX_DIR", _os.path.abspath(
        _os.path.join(HERE, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))
    sys.path.insert(0, _os.path.join(HERE, "..", "..", "target_fishing_v1_freeze", "code"))
    import target_fishing as TF
    _, _, full_df = TF._load()
    popularity = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size().to_dict()

    all_records = []
    for setname in ["A", "B", "D"]:
        all_records.extend(load_capture(setname))
    report["popularity_baseline_combined"] = p0b13_popularity(all_records, popularity, "annotated_targets")

    out_path = os.path.join(RESULTS_DIR, "phase0b_block1_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}", flush=True)
    print(json.dumps(report, indent=2, default=str)[:3000], flush=True)


if __name__ == "__main__":
    main()
