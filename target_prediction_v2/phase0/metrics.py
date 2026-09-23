"""
Phase 0 -- Section 3.2 metrics, implemented against whatever a search
variant returns (a ranked list of {target_chembl, best_similarity,
evidence_score, ...} dicts, exactly target_fishing.py's own result shape).

Deliberately excludes: family-level Top-k, within-family discrimination
(need a target-family/protein-classification mapping not yet pulled -- see
PHASE0_AUTOPSY.md's "not implemented" section) and calibration metrics
(Brier/ECE -- Phase 0 has no calibrated scorer to evaluate yet; that is
Phase 4's job per the procedure doc, Section 6). Implementing everything
else Section 3.2 lists that IS computable now: any-annotated Top-k,
primary-target Top-k/MRR/median-rank, precision@k/recall/targets-per-query,
and the precision-coverage curve.
"""
import statistics


def rank_of_best(results, target_set, rank_field="best_similarity"):
    """Best (lowest) rank among any target in target_set, or None if none
       of them appear in results at all. `results` must already be sorted
       by the ranking key the caller wants scored (v1's own default is
       best_similarity descending -- see target_fishing.py's search())."""
    if not target_set:
        return None
    best = None
    for i, r in enumerate(results):
        if r["target_chembl"] in target_set:
            rank = i + 1
            if best is None or rank < best:
                best = rank
    return best


def any_annotated_topk_hit(results, annotated_targets, k):
    rank = rank_of_best(results, set(annotated_targets))
    return rank is not None and rank <= k


def primary_metrics(per_query_primary_ranks, ks=(1, 5, 10)):
    """per_query_primary_ranks: list of (rank or None) -- the best rank
       among a query's primary/mechanism targets, one entry per query that
       HAD at least one primary target (queries without one are excluded
       upstream, not counted as failures -- Section 3.1: 'never score a
       legitimate secondary target as a failure merely because it is not
       the primary one' extends to 'don't penalize a query for lacking
       primary-target annotation at all, that is a coverage gap, not a
       method failure')."""
    n = len(per_query_primary_ranks)
    if n == 0:
        return {"n_queries_with_primary_target": 0}
    out = {"n_queries_with_primary_target": n}
    for k in ks:
        hit = sum(1 for r in per_query_primary_ranks if r is not None and r <= k)
        out[f"top_{k}"] = round(hit / n, 4)
    mrr = sum((1.0 / r) if r is not None else 0.0 for r in per_query_primary_ranks) / n
    out["mrr"] = round(mrr, 4)
    found_ranks = [r for r in per_query_primary_ranks if r is not None]
    out["median_rank_when_found"] = statistics.median(found_ranks) if found_ranks else None
    out["pct_never_found"] = round(100 * (n - len(found_ranks)) / n, 1)
    return out


def any_annotated_metrics(per_query_ranks, ks=(1, 5, 10)):
    n = len(per_query_ranks)
    if n == 0:
        return {"n_queries": 0}
    out = {"n_queries": n}
    for k in ks:
        hit = sum(1 for r in per_query_ranks if r is not None and r <= k)
        out[f"top_{k}"] = round(hit / n, 4)
    mrr = sum((1.0 / r) if r is not None else 0.0 for r in per_query_ranks) / n
    out["mrr"] = round(mrr, 4)
    return out


def precision_recall_at_threshold(results, annotated_targets):
    """At a fixed similarity threshold, v1's search() already returns
       every target clearing that threshold -- so `results` IS the
       predicted set. Precision = predicted-and-annotated / predicted;
       recall = predicted-and-annotated / annotated. MCC needs a defined
       negative class (all OTHER targets in the reference universe), which
       run_phase0.py supplies as n_targets_universe."""
    predicted = set(r["target_chembl"] for r in results)
    annotated = set(annotated_targets)
    tp = len(predicted & annotated)
    n_pred = len(predicted)
    n_true = len(annotated)
    precision = tp / n_pred if n_pred else None
    recall = tp / n_true if n_true else None
    return {
        "n_predicted": n_pred,
        "n_annotated": n_true,
        "tp": tp,
        "precision": round(precision, 4) if precision is not None else None,
        "recall": round(recall, 4) if recall is not None else None,
    }


def mcc(tp, fp, fn, tn):
    num = tp * tn - fp * fn
    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return round(num / denom, 4) if denom > 0 else None


def aggregate_precision_recall(per_query_pr, n_targets_universe):
    """Micro-averaged precision/recall/MCC + mean predicted-targets-per-
       query, across a list of precision_recall_at_threshold() outputs.
       MCC treats each query's non-predicted, non-annotated targets (out
       of the full target universe) as true negatives -- this is what
       makes MCC meaningful here (a huge, mostly-irrelevant target
       universe per query), matching [1]'s own use of MCC for this task."""
    tp = sum(x["tp"] for x in per_query_pr)
    n_pred_total = sum(x["n_predicted"] for x in per_query_pr)
    n_true_total = sum(x["n_annotated"] for x in per_query_pr)
    fp = n_pred_total - tp
    fn = n_true_total - tp
    n_queries = len(per_query_pr)
    tn = n_queries * n_targets_universe - tp - fp - fn
    return {
        "n_queries": n_queries,
        "micro_precision": round(tp / n_pred_total, 4) if n_pred_total else None,
        "micro_recall": round(tp / n_true_total, 4) if n_true_total else None,
        "mean_predicted_targets_per_query": round(n_pred_total / n_queries, 2) if n_queries else None,
        "mcc": mcc(tp, fp, fn, tn),
    }


def precision_coverage_curve(per_query_results_by_threshold, annotated_by_query, n_targets_universe, thresholds):
    """per_query_results_by_threshold: {threshold: {query_key: results_list}}
       Returns [{threshold, micro_precision, micro_recall,
                 mean_predicted_targets_per_query}, ...] -- the primary
       plot Section 3.2 calls out ('to prevent recall inflation by
       breadth')."""
    curve = []
    for th in thresholds:
        per_query = per_query_results_by_threshold[th]
        prs = [precision_recall_at_threshold(res, annotated_by_query[qk]) for qk, res in per_query.items()]
        agg = aggregate_precision_recall(prs, n_targets_universe)
        agg["threshold"] = th
        curve.append(agg)
    return curve
