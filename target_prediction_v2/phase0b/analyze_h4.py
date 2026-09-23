"""
Phase 0b, Block 2 (P0b-7, P0b-8, P0b-9) -- re-check the H4 claim (weighted
k-NN beats best_similarity) under nested CV, a scaffold-clustered paired
BCa bootstrap CI, and matched-breadth comparison.

METHODOLOGICAL CHOICES MADE HERE, disclosed rather than silent (rev5 does
not hand down exact fold counts / grids -- Section 8 states the protocol
class, not the numbers):
  - Outer: 5-fold, grouped by Bemis-Murcko scaffold (GroupKFold, hand-rolled).
  - Inner: 3-fold, same grouping, for hyperparameter selection within each
    outer-training fold.
  - Grid: k in {5,10,15,25,50,75,100}, alpha in {0.0,0.25,0.5,0.75,1.0,1.5,2.0}.
    alpha=0.0 is INCLUDED deliberately -- it is exactly the unweighted-vote
    baseline B, so nested CV can select "no weighting helps" if that is
    what the data says, rather than assuming weighting wins by construction.
  - Selection metric: mean MRR on ANY-ANNOTATED targets, on inner-validation
    folds. Chosen over primary-target MRR for tuning stability -- only 64
    (Set A) / 57 (Set B) queries have a primary target at all, and tuning
    on that few points per inner fold (~13-19) risks selecting noise.
    Primary-target metrics are still MEASURED (not used to tune) and
    reported as the headline outcome, since that is what Phase 0 actually
    claimed the gain on.
  - Every outer-test query's rank comes from a (k, alpha) chosen WITHOUT
    that query ever being in the fold that selected it -- this is what
    makes the resulting numbers an honest out-of-sample estimate, unlike
    Phase 0's original H4 grid sweep which reported the best grid cell's
    in-sample performance.

Sample-size caveat, stated up front and again in the output: Set A's
primary-target tier is only 64 queries; a 5-fold split puts ~13 per outer
fold. The BCa CI on primary-target deltas will be WIDE. This is a real
constraint of the data, not a bug -- reported honestly per rev5 Section 8.9
("a mean precision computed on 81 points needs its n printed beside it").
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

K_GRID = [5, 10, 15, 25, 50, 75, 100]
ALPHA_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
N_OUTER = 5
N_INNER = 3
SEED = 42


def load_capture(setname):
    with open(os.path.join(RESULTS_DIR, f"capture_{setname}.jsonl")) as f:
        return [json.loads(line) for line in f]


def load_scaffolds(records):
    phase0_dir = os.path.join(HERE, "..", "phase0")
    sys.path.insert(0, phase0_dir)
    sys.path.insert(0, os.path.join(phase0_dir, "..", "..", "target_fishing_v1_freeze", "code"))
    os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
        os.path.join(phase0_dir, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))
    import target_fishing as TF
    _, _, full_df = TF._load()
    needed = set(r["smiles"] for r in records)
    sub = full_df[full_df["smiles"].isin(needed)][["smiles", "murcko_scaffold"]].drop_duplicates(subset="smiles")
    scaffold_by_smiles = dict(zip(sub["smiles"], sub["murcko_scaffold"]))
    out = []
    for r in records:
        sc = scaffold_by_smiles.get(r["smiles"])
        if not isinstance(sc, str) or not sc:
            sc = f"__no_ring__:{r['smiles']}"
        out.append(sc)
    return out


def group_kfold(n_items, group_of_item, k, seed):
    """Scaffold-grouped K-fold: whole scaffold groups assigned to folds
       (greedy balancing by group size), so no scaffold straddles a
       train/test boundary within a split."""
    groups = {}
    for i, g in enumerate(group_of_item):
        groups.setdefault(g, []).append(i)
    group_keys = list(groups.keys())
    rng = random.Random(seed)
    rng.shuffle(group_keys)
    fold_sizes = [0] * k
    fold_of_group = {}
    # greedy: assign each group to the currently-smallest fold
    for g in sorted(group_keys, key=lambda gk: -len(groups[gk])):
        smallest = min(range(k), key=lambda f: fold_sizes[f])
        fold_of_group[g] = smallest
        fold_sizes[smallest] += len(groups[g])
    folds = [[] for _ in range(k)]
    for i, g in enumerate(group_of_item):
        folds[fold_of_group[g]].append(i)
    return folds


def mrr_of(records_subset, k, alpha, tier):
    total = 0.0
    n = 0
    for r in records_subset:
        scores = S.score_query(r["neighbours"], k=k, alpha=alpha)
        ranked = S.rank_from_scores(scores)
        rank = S.best_rank(ranked, r[tier])
        total += (1.0 / rank) if rank is not None else 0.0
        n += 1
    return total / n if n else 0.0


def select_hyperparams(train_records, train_group, seed):
    inner_folds_idx = group_kfold(len(train_records), train_group, N_INNER, seed)
    best = None
    for k in K_GRID:
        for alpha in ALPHA_GRID:
            scores = []
            for fi in range(N_INNER):
                val_idx = set(inner_folds_idx[fi])
                val_records = [train_records[i] for i in range(len(train_records)) if i in val_idx]
                if not val_records:
                    continue
                scores.append(mrr_of(val_records, k, alpha, "annotated_targets"))
            mean_score = sum(scores) / len(scores) if scores else 0.0
            if best is None or mean_score > best[0]:
                best = (mean_score, k, alpha)
    return best[1], best[2], best[0]


def nested_cv(records, scaffolds, seed=SEED):
    outer_folds_idx = group_kfold(len(records), scaffolds, N_OUTER, seed)
    out_of_sample = {}  # query index -> {"k":..,"alpha":.., "rank_any":.., "rank_primary":..}
    fold_choices = []
    for fo in range(N_OUTER):
        test_idx = set(outer_folds_idx[fo])
        train_idx = [i for i in range(len(records)) if i not in test_idx]
        train_records = [records[i] for i in train_idx]
        train_group = [scaffolds[i] for i in train_idx]
        k, alpha, inner_mrr = select_hyperparams(train_records, train_group, seed + fo)
        fold_choices.append({"outer_fold": fo, "k": k, "alpha": alpha, "inner_val_mrr": round(inner_mrr, 4),
                              "n_train": len(train_idx), "n_test": len(test_idx)})
        for i in test_idx:
            r = records[i]
            scores = S.score_query(r["neighbours"], k=k, alpha=alpha)
            ranked = S.rank_from_scores(scores)
            rank_any = S.best_rank(ranked, r["annotated_targets"])
            rank_prim = S.best_rank(ranked, r["primary_targets"]) if r["primary_targets"] else None
            out_of_sample[i] = {"k": k, "alpha": alpha, "rank_any": rank_any, "rank_primary": rank_prim}
    return out_of_sample, fold_choices


def fixed_baseline_ranks(records, k, alpha):
    out = {}
    for i, r in enumerate(records):
        scores = S.score_query(r["neighbours"], k=k, alpha=alpha)
        ranked = S.rank_from_scores(scores)
        rank_any = S.best_rank(ranked, r["annotated_targets"])
        rank_prim = S.best_rank(ranked, r["primary_targets"]) if r["primary_targets"] else None
        out[i] = {"rank_any": rank_any, "rank_primary": rank_prim}
    return out


def bestsim_ranks(records):
    out = {}
    for i, r in enumerate(records):
        ranked = sorted(r["bestsim_results"], key=lambda ts: -ts[1])
        rank_map = {t: j + 1 for j, (t, s) in enumerate(ranked)}
        rank_any = min((rank_map[t] for t in r["annotated_targets"] if t in rank_map), default=None)
        rank_prim = (min((rank_map[t] for t in r["primary_targets"] if t in rank_map), default=None)
                     if r["primary_targets"] else None)
        out[i] = {"rank_any": rank_any, "rank_primary": rank_prim}
    return out


def topk_hit(rank, k):
    return 1.0 if (rank is not None and rank <= k) else 0.0


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def summarize_tier(ranks_dict, key, indices):
    n = len(indices)
    if n == 0:
        return {"n": 0}
    top1 = sum(topk_hit(ranks_dict[i][key], 1) for i in indices) / n
    top10 = sum(topk_hit(ranks_dict[i][key], 10) for i in indices) / n
    mrr = sum(rr(ranks_dict[i][key]) for i in indices) / n
    return {"n": n, "top_1": round(top1, 4), "top_10": round(top10, 4), "mrr": round(mrr, 4)}


def paired_bca(ranks_a, ranks_b, key, metric_fn, indices, scaffolds):
    diffs = [metric_fn(ranks_a[i][key]) - metric_fn(ranks_b[i][key]) for i in indices]
    groups = [scaffolds[i] for i in indices]
    return BCA.scaffold_clustered_bca(diffs, groups)


def run_for_set(setname):
    records = load_capture(setname)
    scaffolds = load_scaffolds(records)
    print(f"Set {setname}: {len(records)} queries, {len(set(scaffolds))} distinct scaffold groups", flush=True)

    print("  Running nested CV (weighted k-NN, tuned)...", flush=True)
    oos, fold_choices = nested_cv(records, scaffolds)
    weighted_ranks = {i: {"rank_any": oos[i]["rank_any"], "rank_primary": oos[i]["rank_primary"]} for i in oos}

    print("  Computing fixed baselines (10-NN unweighted, best_similarity)...", flush=True)
    tenn_ranks = fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = bestsim_ranks(records)

    all_idx = list(range(len(records)))
    primary_idx = [i for i, r in enumerate(records) if r["primary_targets"]]

    out = {
        "n_queries": len(records),
        "n_scaffold_groups": len(set(scaffolds)),
        "n_queries_with_primary_target": len(primary_idx),
        "fold_choices": fold_choices,
        "any_annotated": {
            "weighted_knn_nested_cv": summarize_tier(weighted_ranks, "rank_any", all_idx),
            "unweighted_10nn": summarize_tier(tenn_ranks, "rank_any", all_idx),
            "best_similarity": summarize_tier(bs_ranks, "rank_any", all_idx),
        },
        "primary_target": {
            "weighted_knn_nested_cv": summarize_tier(weighted_ranks, "rank_primary", primary_idx),
            "unweighted_10nn": summarize_tier(tenn_ranks, "rank_primary", primary_idx),
            "best_similarity": summarize_tier(bs_ranks, "rank_primary", primary_idx),
        },
        "paired_bca_top1_vs_best_similarity": {
            "any_annotated": paired_bca(weighted_ranks, bs_ranks, "rank_any", lambda r: topk_hit(r, 1), all_idx, scaffolds),
            "primary_target": (paired_bca(weighted_ranks, bs_ranks, "rank_primary", lambda r: topk_hit(r, 1), primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_mrr_vs_best_similarity": {
            "any_annotated": paired_bca(weighted_ranks, bs_ranks, "rank_any", rr, all_idx, scaffolds),
            "primary_target": (paired_bca(weighted_ranks, bs_ranks, "rank_primary", rr, primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_top1_vs_10nn": {
            "any_annotated": paired_bca(weighted_ranks, tenn_ranks, "rank_any", lambda r: topk_hit(r, 1), all_idx, scaffolds),
            "primary_target": (paired_bca(weighted_ranks, tenn_ranks, "rank_primary", lambda r: topk_hit(r, 1), primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_mrr_vs_10nn": {
            "any_annotated": paired_bca(weighted_ranks, tenn_ranks, "rank_any", rr, all_idx, scaffolds),
            "primary_target": (paired_bca(weighted_ranks, tenn_ranks, "rank_primary", rr, primary_idx, scaffolds)
                                if primary_idx else None),
        },
    }
    return out


def main():
    report = {}
    for setname in ["A", "B"]:
        report[setname] = run_for_set(setname)
        out_path = os.path.join(RESULTS_DIR, "phase0b_h4_report.json")
        with open(out_path, "w") as f:
            json.dump(report, f, indent=2, default=str)
        print(f"  saved -> {out_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
