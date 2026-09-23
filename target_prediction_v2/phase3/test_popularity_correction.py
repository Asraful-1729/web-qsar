"""
Phase 3A: actually TEST score.py's popularity_correct() -- BUILD_PLAN.md
lists this as an "untested lever, still planned" (the function exists,
was never run against real data). This applies it to the same n=1600
pooled v2 captures used for the density-adaptive re-fit, at k=10
unweighted pooling (alpha=0, matching the density-adaptive rule's own
comparator), and measures whether it actually improves retrieval --
not assumed, tested.

popularity_correct(scores, popularity, eps) divides each target's raw
k-NN vote-mass score by (global_popularity + eps) -- rev5 Section 3's
"same logic as SEA's background model." Target popularity = distinct
compound count per target in v2_index (real range: 1 to 195,809,
confirmed a severe skew before running this test).

Tested at multiple eps values (not just the function's default eps=1.0,
which was never validated) -- and stratified by the JUST-refit density
threshold (density>=8, phase3/PHASE3A_DENSITY_REFIT.md), since
popularity correction is scoped to layer on top of the pooling branch,
which only fires in the dense regime.

Usage: python3 test_popularity_correction.py
"""
import json
import os
import sys

PHASE0B_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0b")
sys.path.insert(0, PHASE0B_DIR)
import score as S  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

EPS_GRID = [0.1, 1.0, 5.0, 20.0, 100.0]
K = 10


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_records():
    records = []
    for letter in ("A", "B"):
        with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_v2_{letter}.jsonl")) as f:
            for line in f:
                r = json.loads(line)
                r["_set"] = letter
                records.append(r)
    return records


def main():
    with open(os.path.join(RESULTS_DIR, "target_popularity_v2.json")) as f:
        popularity = json.load(f)
    print(f"{len(popularity)} targets with known popularity", flush=True)

    records = load_records()
    print(f"{len(records)} pooled queries loaded", flush=True)

    scaffolds = [f"{r['_set']}:{r['molecule_chembl_id']}" for r in records]  # per-query grouping fallback
    # reuse the same scaffold grouping the density refit used, for a proper clustered CI
    sys.path.insert(0, os.path.join(HERE))
    from fit_density_adaptive_rule_v2 import load_scaffolds_v2
    real_scaffolds = load_scaffolds_v2(records)
    scaffold_groups = [f"{r['_set']}:{sc}" for r, sc in zip(records, real_scaffolds)]

    densities = [density_of(r) for r in records]

    raw_scores_by_query = []
    for r in records:
        raw_scores_by_query.append(S.score_query(r["neighbours"], k=K, alpha=0))

    results = {}
    for eps in EPS_GRID:
        diffs_any_all, diffs_prim_all, groups_prim_all = [], [], []
        diffs_any_dense, diffs_prim_dense, groups_prim_dense = [], [], []
        groups_dense = []
        for i, r in enumerate(records):
            raw = raw_scores_by_query[i]
            corrected = S.popularity_correct(raw, popularity, eps=eps)
            raw_ranked = S.rank_from_scores(raw)
            corr_ranked = S.rank_from_scores(corrected)

            raw_rank_any = S.best_rank(raw_ranked, r["annotated_targets"])
            corr_rank_any = S.best_rank(corr_ranked, r["annotated_targets"])
            diffs_any_all.append(rr(corr_rank_any) - rr(raw_rank_any))

            if r["primary_targets"]:
                raw_rank_prim = S.best_rank(raw_ranked, r["primary_targets"])
                corr_rank_prim = S.best_rank(corr_ranked, r["primary_targets"])
                diffs_prim_all.append(rr(corr_rank_prim) - rr(raw_rank_prim))
                groups_prim_all.append(scaffold_groups[i])

            if densities[i] >= 8:  # the just-refit pooling threshold
                diffs_any_dense.append(rr(corr_rank_any) - rr(raw_rank_any))
                groups_dense.append(scaffold_groups[i])
                if r["primary_targets"]:
                    diffs_prim_dense.append(rr(corr_rank_prim) - rr(raw_rank_prim))
                    groups_prim_dense.append(scaffold_groups[i])

        ci_any_all = BCA.scaffold_clustered_bca(diffs_any_all, scaffold_groups)
        ci_any_dense = BCA.scaffold_clustered_bca(diffs_any_dense, groups_dense) if diffs_any_dense else None
        ci_prim_all = BCA.scaffold_clustered_bca(diffs_prim_all, groups_prim_all) if diffs_prim_all else None
        ci_prim_dense = BCA.scaffold_clustered_bca(diffs_prim_dense, groups_prim_dense) if diffs_prim_dense else None

        results[eps] = {
            "n_any_all": len(diffs_any_all),
            "mean_diff_any_all": round(sum(diffs_any_all) / len(diffs_any_all), 5),
            "ci_any_all": ci_any_all,
            "n_any_dense": len(diffs_any_dense),
            "mean_diff_any_dense": round(sum(diffs_any_dense) / len(diffs_any_dense), 5) if diffs_any_dense else None,
            "ci_any_dense": ci_any_dense,
            "n_prim_all": len(diffs_prim_all),
            "mean_diff_prim_all": round(sum(diffs_prim_all) / len(diffs_prim_all), 5) if diffs_prim_all else None,
            "ci_prim_all": ci_prim_all,
            "n_prim_dense": len(diffs_prim_dense),
            "mean_diff_prim_dense": round(sum(diffs_prim_dense) / len(diffs_prim_dense), 5) if diffs_prim_dense else None,
            "ci_prim_dense": ci_prim_dense,
        }
        print(f"eps={eps}: any-tier all n={len(diffs_any_all)} mean_diff={results[eps]['mean_diff_any_all']} "
              f"CI={ci_any_all}", flush=True)
        print(f"         any-tier dense (density>=8) n={len(diffs_any_dense)} "
              f"mean_diff={results[eps]['mean_diff_any_dense']} CI={ci_any_dense}", flush=True)
        print(f"         primary-tier all n={len(diffs_prim_all)} mean_diff={results[eps]['mean_diff_prim_all']} "
              f"CI={ci_prim_all}", flush=True)
        print(f"         primary-tier dense (density>=8) n={len(diffs_prim_dense)} "
              f"mean_diff={results[eps]['mean_diff_prim_dense']} CI={ci_prim_dense}", flush=True)

    out_path = os.path.join(RESULTS_DIR, "popularity_correction_test.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
