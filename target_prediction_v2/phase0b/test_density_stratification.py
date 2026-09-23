"""
Direct test of the density-stratification hypothesis from
PHASE0B_DIAGNOSIS_A_VS_B.md: "Set A vs Set B" may be a proxy for
neighbourhood density, not the real causal variable behind the H4a
(pooling-beats-best-similarity) effect being significant for Set B and not
for Set A. If density is the real variable, splitting SET A's OWN 200
queries by density should reproduce the same split in outcome: pooling
should show a real, positive, significant effect on Set A's dense subset
(matching Set B's pattern) and no effect (or a null one) on the sparse
subset (matching Set A's overall pattern) -- WITHOUT needing to invoke
"approved vs clinical" as an explanation at all.

Uses `capture_scaffold_strict_A.jsonl` (300-neighbour, scaffold-strict
leakage control -- the stricter, Rev-6-established-correct control, so
this test is run under the same conditions the H4a/H4b conclusions
actually rest on, not the looser near-duplicate-only control).

Density metric: n_neighbours at Tanimoto >= 0.5 (the single strongest,
most mechanistically direct within-set predictor found in the diagnosis --
"how many decent analogs are actually available to pool votes over" is
exactly the mechanism the density hypothesis proposes). max_similarity_to_
index is reported alongside as a robustness cross-check (it should split
the same way if the effect is genuinely about density, not an artifact of
which specific metric was used to define the strata).

Two views:
  1. A clean median split (dense half vs sparse half) -- the primary
     significance test, adequately powered (100 queries/scaffold-clustered
     BCa per half).
  2. A quartile dose-response view -- weaker per-quartile power (~50
     queries each) but shows whether the effect is a step function or a
     gradient, which a single median split cannot distinguish.

Usage: python3 test_density_stratification.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S  # noqa: E402
import bca as BCA  # noqa: E402
import analyze_h4 as H4  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")


def load_capture_scaffold_strict_A():
    with open(os.path.join(RESULTS_DIR, "capture_scaffold_strict_A.jsonl")) as f:
        return [json.loads(line) for line in f]


def density_metric(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rank_sets(records):
    tenn_ranks = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = H4.bestsim_ranks(records)
    return tenn_ranks, bs_ranks


def stratum_report(records, scaffolds, label):
    tenn_ranks, bs_ranks = rank_sets(records)
    all_idx = list(range(len(records)))
    primary_idx = [i for i, r in enumerate(records) if r["primary_targets"]]

    def topk_hit(r, k):
        return 1.0 if (r is not None and r <= k) else 0.0

    def rr(r):
        return (1.0 / r) if r is not None else 0.0

    any_10nn = H4.summarize_tier(tenn_ranks, "rank_any", all_idx)
    any_bs = H4.summarize_tier(bs_ranks, "rank_any", all_idx)

    bca_top1 = H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", lambda r: topk_hit(r, 1), all_idx, scaffolds)
    bca_mrr = H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", rr, all_idx, scaffolds)

    prim_bca_top1 = (H4.paired_bca(tenn_ranks, bs_ranks, "rank_primary", lambda r: topk_hit(r, 1), primary_idx, scaffolds)
                      if len(primary_idx) >= 10 else {"note": f"too few primary-target queries (n={len(primary_idx)}) for a meaningful bootstrap"})

    return {
        "label": label,
        "n_queries": len(records),
        "n_queries_with_primary_target": len(primary_idx),
        "mean_density_n_ge_0.5": round(sum(density_metric(r) for r in records) / len(records), 2),
        "mean_max_similarity": round(sum(r["max_similarity_to_index"] for r in records) / len(records), 4),
        "any_annotated_unweighted_10nn": any_10nn,
        "any_annotated_best_similarity": any_bs,
        "paired_bca_10nn_vs_bestsim_top1_any": bca_top1,
        "paired_bca_10nn_vs_bestsim_mrr_any": bca_mrr,
        "paired_bca_10nn_vs_bestsim_top1_primary": prim_bca_top1,
    }


def main():
    records = load_capture_scaffold_strict_A()
    scaffolds_all = H4.load_scaffolds(records)

    # attach density + scaffold to each record for sorting/splitting
    indexed = list(enumerate(records))
    by_density = sorted(indexed, key=lambda kv: density_metric(kv[1]))

    n = len(records)
    report = {}

    # --- Median split ---
    sparse_idx = [i for i, r in by_density[: n // 2]]
    dense_idx = [i for i, r in by_density[n // 2:]]
    sparse_records = [records[i] for i in sparse_idx]
    dense_records = [records[i] for i in dense_idx]
    sparse_scaffolds = [scaffolds_all[i] for i in sparse_idx]
    dense_scaffolds = [scaffolds_all[i] for i in dense_idx]

    print("Computing median-split strata (sparse half vs dense half)...", flush=True)
    report["median_split"] = {
        "sparse_half": stratum_report(sparse_records, sparse_scaffolds, "Set A sparse half (n_neighbours>=0.5 below median)"),
        "dense_half": stratum_report(dense_records, dense_scaffolds, "Set A dense half (n_neighbours>=0.5 at/above median)"),
    }

    # --- Quartiles ---
    print("Computing quartile strata (dose-response)...", flush=True)
    q_size = n // 4
    quartiles = [by_density[0:q_size], by_density[q_size:2*q_size], by_density[2*q_size:3*q_size], by_density[3*q_size:]]
    report["quartiles"] = {}
    for qi, q in enumerate(quartiles):
        q_idx = [i for i, r in q]
        q_records = [records[i] for i in q_idx]
        q_scaffolds = [scaffolds_all[i] for i in q_idx]
        report["quartiles"][f"Q{qi+1}"] = stratum_report(q_records, q_scaffolds, f"Set A density quartile {qi+1} (Q1=sparsest)")

    out_path = os.path.join(RESULTS_DIR, "density_stratification_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
