"""
Widened-grid rerun, requested directly ("widen the k and alpha grids and
rerun") after the Rev 6 response found k hitting the old grid's edge
(k=100, in 3 of 5 outer folds for Set A under scaffold-strict leakage) --
the alpha widening alone (analyze_h4_rev6.py) wasn't the full story;
k needed the same treatment, which in turn required recapturing with a
larger NEIGHBOR_CAP (capture.py / capture_scaffold_strict.py, 150 -> 300)
since k cannot usefully exceed how many neighbours were actually stored
per query.

Grids:
  k:     {5, 10, 15, 25, 50, 75, 100, 150, 200, 300} -- was capped at 100
  alpha: {0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24}  -- extends the
         rev6 grid's max (12) further, even though rev6 already found an
         interior optimum around 3-8, to leave a larger safety margin
         before calling that settled.

Reads the SAME capture files as analyze_h4_rev6.py (capture_{A,B}.jsonl,
capture_scaffold_strict_{A,B}.jsonl), regenerated with NEIGHBOR_CAP=300 so
k=300 is actually meaningful rather than silently truncating to whatever
was stored. Writes to a NEW output file (phase0b_h4_widegrid_report.json)
rather than overwriting the rev6 report, preserving that run's provenance.

Usage: python3 analyze_h4_widegrid.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402
import analyze_h4 as H4  # noqa: E402
import analyze_h4_rev6 as H4R6  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

K_GRID = [5, 10, 15, 25, 50, 75, 100, 150, 200, 300]
ALPHA_GRID = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
N_OUTER = H4.N_OUTER
N_INNER = H4.N_INNER
SEED = H4.SEED


def select_hyperparams_wide(train_records, train_group, seed):
    inner_folds_idx = H4.group_kfold(len(train_records), train_group, N_INNER, seed)
    best = None
    for k in K_GRID:
        for alpha in ALPHA_GRID:
            scores = []
            for fi in range(N_INNER):
                val_idx = set(inner_folds_idx[fi])
                val_records = [train_records[i] for i in range(len(train_records)) if i in val_idx]
                if not val_records:
                    continue
                scores.append(H4.mrr_of(val_records, k, alpha, "annotated_targets"))
            mean_score = sum(scores) / len(scores) if scores else 0.0
            if best is None or mean_score > best[0]:
                best = (mean_score, k, alpha)
    return best[1], best[2], best[0]


def nested_cv_wide(records, scaffolds, seed=SEED):
    outer_folds_idx = H4.group_kfold(len(records), scaffolds, N_OUTER, seed)
    out_of_sample = {}
    fold_choices = []
    for fo in range(N_OUTER):
        test_idx = set(outer_folds_idx[fo])
        train_idx = [i for i in range(len(records)) if i not in test_idx]
        train_records = [records[i] for i in train_idx]
        train_group = [scaffolds[i] for i in train_idx]
        k, alpha, inner_mrr = select_hyperparams_wide(train_records, train_group, seed + fo)
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


def run_variant(setname, variant):
    records = H4R6.load_capture_variant(setname, variant)
    scaffolds = H4.load_scaffolds(records)
    max_neighbours_seen = max(len(r["neighbours"]) for r in records)
    print(f"Set {setname} [{variant}]: {len(records)} queries, {len(set(scaffolds))} scaffold groups, "
          f"max neighbours stored per query = {max_neighbours_seen}", flush=True)

    oos, fold_choices = nested_cv_wide(records, scaffolds)
    weighted_ranks = {i: {"rank_any": oos[i]["rank_any"], "rank_primary": oos[i]["rank_primary"]} for i in oos}
    tenn_ranks = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = H4.bestsim_ranks(records)

    all_idx = list(range(len(records)))
    primary_idx = [i for i, r in enumerate(records) if r["primary_targets"]]

    out = {
        "variant": variant, "n_queries": len(records), "n_scaffold_groups": len(set(scaffolds)),
        "n_queries_with_primary_target": len(primary_idx),
        "max_neighbours_stored_per_query": max_neighbours_seen,
        "k_grid": K_GRID, "alpha_grid": ALPHA_GRID,
        "fold_choices": fold_choices,
        "any_annotated": {
            "weighted_knn_nested_cv": H4.summarize_tier(weighted_ranks, "rank_any", all_idx),
            "unweighted_10nn": H4.summarize_tier(tenn_ranks, "rank_any", all_idx),
            "best_similarity": H4.summarize_tier(bs_ranks, "rank_any", all_idx),
        },
        "primary_target": {
            "weighted_knn_nested_cv": H4.summarize_tier(weighted_ranks, "rank_primary", primary_idx),
            "unweighted_10nn": H4.summarize_tier(tenn_ranks, "rank_primary", primary_idx),
            "best_similarity": H4.summarize_tier(bs_ranks, "rank_primary", primary_idx),
        },
        "paired_bca_top1_vs_best_similarity": {
            "any_annotated": H4.paired_bca(weighted_ranks, bs_ranks, "rank_any", lambda r: H4.topk_hit(r, 1), all_idx, scaffolds),
            "primary_target": (H4.paired_bca(weighted_ranks, bs_ranks, "rank_primary", lambda r: H4.topk_hit(r, 1), primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_mrr_vs_best_similarity": {
            "any_annotated": H4.paired_bca(weighted_ranks, bs_ranks, "rank_any", H4.rr, all_idx, scaffolds),
            "primary_target": (H4.paired_bca(weighted_ranks, bs_ranks, "rank_primary", H4.rr, primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_top1_vs_10nn": {
            "any_annotated": H4.paired_bca(weighted_ranks, tenn_ranks, "rank_any", lambda r: H4.topk_hit(r, 1), all_idx, scaffolds),
            "primary_target": (H4.paired_bca(weighted_ranks, tenn_ranks, "rank_primary", lambda r: H4.topk_hit(r, 1), primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_mrr_vs_10nn": {
            "any_annotated": H4.paired_bca(weighted_ranks, tenn_ranks, "rank_any", H4.rr, all_idx, scaffolds),
            "primary_target": (H4.paired_bca(weighted_ranks, tenn_ranks, "rank_primary", H4.rr, primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_10nn_vs_best_similarity_top1": {
            "any_annotated": H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", lambda r: H4.topk_hit(r, 1), all_idx, scaffolds),
            "primary_target": (H4.paired_bca(tenn_ranks, bs_ranks, "rank_primary", lambda r: H4.topk_hit(r, 1), primary_idx, scaffolds)
                                if primary_idx else None),
        },
        "paired_bca_10nn_vs_best_similarity_mrr": {
            "any_annotated": H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", H4.rr, all_idx, scaffolds),
            "primary_target": (H4.paired_bca(tenn_ranks, bs_ranks, "rank_primary", H4.rr, primary_idx, scaffolds)
                                if primary_idx else None),
        },
    }
    return out


def main():
    report = {}
    for setname in ["A", "B"]:
        report[setname] = {}
        for variant in ["near_dup_only", "scaffold_strict"]:
            report[setname][variant] = run_variant(setname, variant)
            out_path = os.path.join(RESULTS_DIR, "phase0b_h4_widegrid_report.json")
            with open(out_path, "w") as f:
                json.dump(report, f, indent=2, default=str)
            print(f"  saved -> {out_path}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
