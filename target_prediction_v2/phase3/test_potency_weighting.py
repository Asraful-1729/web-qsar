"""
Phase 3A: test score.py's potency_weight() (logistic ramp on pchembl,
NOT the s^alpha similarity weighting H4b already closed as a genuine
null -- a different lever, operating on potency not similarity, still
open per BUILD_PLAN.md's architecture list). Tests potency-weighted k=10
pooling against plain unweighted k=10 pooling, same n=1600 v2 captures
and scaffold-clustered BCa methodology as the density refit and the
popularity-correction test.

Grid: rev5 Section 3 flags the functional form (threshold, floor,
steepness) as "unvalidated -- tune it" -- swept here rather than assumed
at the function's default (threshold=6.0, floor=0.2, steepness=2.0).

Usage: python3 test_potency_weighting.py
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

# (threshold, floor, steepness)
GRID = [
    (6.0, 0.2, 2.0),   # function default
    (5.0, 0.2, 2.0),
    (7.0, 0.2, 2.0),
    (6.0, 0.0, 2.0),   # true floor (hard-ish at the low end)
    (6.0, 0.5, 2.0),   # gentler floor
    (6.0, 0.2, 0.5),   # gentler ramp
    (6.0, 0.2, 5.0),   # steeper ramp
]


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
    records = load_records()
    print(f"{len(records)} pooled queries loaded", flush=True)

    sys.path.insert(0, HERE)
    from fit_density_adaptive_rule_v2 import load_scaffolds_v2
    real_scaffolds = load_scaffolds_v2(records)
    scaffold_groups = [f"{r['_set']}:{sc}" for r, sc in zip(records, real_scaffolds)]
    densities = [density_of(r) for r in records]

    raw_scores_by_query = [S.score_query(r["neighbours"], k=10, alpha=0) for r in records]

    results = {}
    for (thr, floor, steep) in GRID:
        key = f"thr={thr}_floor={floor}_steep={steep}"

        def pfn(pchembl, thr=thr, floor=floor, steep=steep):
            return S.potency_weight(pchembl, threshold=thr, floor=floor, steepness=steep)

        diffs_any_all, diffs_any_dense, groups_dense = [], [], []
        diffs_prim_all, groups_prim_all = [], []
        for i, r in enumerate(records):
            weighted = S.score_query(r["neighbours"], k=10, alpha=0, potency_fn=pfn)
            raw = raw_scores_by_query[i]
            raw_ranked = S.rank_from_scores(raw)
            w_ranked = S.rank_from_scores(weighted)

            raw_rank_any = S.best_rank(raw_ranked, r["annotated_targets"])
            w_rank_any = S.best_rank(w_ranked, r["annotated_targets"])
            diffs_any_all.append(rr(w_rank_any) - rr(raw_rank_any))

            if r["primary_targets"]:
                raw_rank_prim = S.best_rank(raw_ranked, r["primary_targets"])
                w_rank_prim = S.best_rank(w_ranked, r["primary_targets"])
                diffs_prim_all.append(rr(w_rank_prim) - rr(raw_rank_prim))
                groups_prim_all.append(scaffold_groups[i])

            if densities[i] >= 8:
                diffs_any_dense.append(rr(w_rank_any) - rr(raw_rank_any))
                groups_dense.append(scaffold_groups[i])

        ci_any_all = BCA.scaffold_clustered_bca(diffs_any_all, scaffold_groups)
        ci_any_dense = BCA.scaffold_clustered_bca(diffs_any_dense, groups_dense) if diffs_any_dense else None
        ci_prim_all = BCA.scaffold_clustered_bca(diffs_prim_all, groups_prim_all) if diffs_prim_all else None

        results[key] = {"ci_any_all": ci_any_all, "ci_any_dense": ci_any_dense, "ci_prim_all": ci_prim_all}
        print(f"{key}: any_all={ci_any_all['point_estimate']:.4f} "
              f"[{ci_any_all['ci_low']:.4f},{ci_any_all['ci_high']:.4f}] excl0={ci_any_all['ci_excludes_zero']} | "
              f"any_dense={ci_any_dense['point_estimate']:.4f} excl0={ci_any_dense['ci_excludes_zero']} | "
              f"prim_all={ci_prim_all['point_estimate']:.4f} excl0={ci_prim_all['ci_excludes_zero']}", flush=True)

    out_path = os.path.join(RESULTS_DIR, "potency_weighting_test.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
