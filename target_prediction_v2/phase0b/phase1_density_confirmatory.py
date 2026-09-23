"""
Phase 1 entry requirement (rev7 Section 2.2, items 1-3 / BUILD_PLAN.md
Phase 1 items 8-10): confirmatory re-test of the density finding at a
properly-powered scale, replacing Phase 0b's 200-query quartile split.

Three things, all on the SAME n=800/set scaffold-strict captures:

1. DECILE (not quartile) dose-response of pooling benefit (unweighted
   10-NN vs. best-similarity) vs. neighbourhood density -- resolves
   whether the Q2 "valley" found at quartile resolution on 200 queries is
   a real local minimum or a small-sample artifact (rev7 Section 1.3).

2. Similarity-SPREAD diagnostic, run at the same scale, both overall and
   specifically within whichever decile(s) replace the old Q2 band --
   properly powered version of the n=37/half preliminary check in
   PHASE0B_REV7_RESPONSE.md.

3. Joint density + target-reference-depth check, done the RIGHT way this
   time (rev7 flag-3): NOT a flat linear model (Phase 0b's attempt was
   the wrong tool for a non-monotonic relationship and settled nothing).
   Instead: within EACH density decile, check whether target-reference-
   depth still separates recovered/not-recovered queries. If depth's
   within-decile predictive gap is small and roughly constant across
   deciles, that supports "depth doesn't add much beyond density." If a
   later decile shows a large depth gap, that decile's queries are still
   bottlenecked on target-side evidence even though neighbourhood density
   is high, meaning the two variables are NOT interchangeable there.

Usage: python3 phase1_density_confirmatory.py --set A
       python3 phase1_density_confirmatory.py --set B
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))

import target_fishing as TF  # noqa: E402
import analyze_h4 as H4  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
N_DECILES = 10


def load_capture(setname):
    with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_{setname}.jsonl")) as f:
        return [json.loads(line) for line in f]


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def spread_gap_of(record):
    sims = sorted((sim for _, sim, _ in record["neighbours"]), reverse=True)
    return (sims[0] - sims[9]) if len(sims) >= 10 else None


def topk_hit(r, k):
    return 1.0 if (r is not None and r <= k) else 0.0


def rr(r):
    return (1.0 / r) if r is not None else 0.0


def decile_report(records, scaffolds, label):
    tenn_ranks = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = H4.bestsim_ranks(records)
    all_idx = list(range(len(records)))
    primary_idx = [i for i, r in enumerate(records) if r["primary_targets"]]

    top1_bca = H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", lambda r: topk_hit(r, 1), all_idx, scaffolds)
    mrr_bca = H4.paired_bca(tenn_ranks, bs_ranks, "rank_any", rr, all_idx, scaffolds)
    prim_top1_bca = (H4.paired_bca(tenn_ranks, bs_ranks, "rank_primary", lambda r: topk_hit(r, 1), primary_idx, scaffolds)
                      if len(primary_idx) >= 15 else None)
    return {
        "label": label, "n_queries": len(records), "n_with_primary": len(primary_idx),
        "mean_density": round(float(np.mean([density_of(r) for r in records])), 2),
        "paired_bca_top1_any": top1_bca, "paired_bca_mrr_any": mrr_bca,
        "paired_bca_top1_primary": prim_top1_bca,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["A", "B"])
    args = ap.parse_args()

    records = load_capture(args.set)
    scaffolds_all = H4.load_scaffolds(records)
    print(f"Set {args.set}: {len(records)} queries loaded", flush=True)

    _, _, full_df = TF._load()
    depth_by_target = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size().to_dict()

    indexed = list(enumerate(records))
    by_density = sorted(indexed, key=lambda kv: density_of(kv[1]))
    n = len(records)
    bin_size = n // N_DECILES

    print(f"Computing {N_DECILES}-decile dose-response...", flush=True)
    deciles = {}
    for d in range(N_DECILES):
        lo = d * bin_size
        hi = (d + 1) * bin_size if d < N_DECILES - 1 else n
        bucket = by_density[lo:hi]
        idx = [i for i, r in bucket]
        recs = [records[i] for i in idx]
        scafs = [scaffolds_all[i] for i in idx]
        deciles[f"D{d+1}"] = decile_report(recs, scafs, f"Set {args.set} density decile {d+1}/10 (D1=sparsest)")
        print(f"  D{d+1}: n={len(recs)}, mean_density={deciles[f'D{d+1}']['mean_density']}, "
              f"top1_delta={deciles[f'D{d+1}']['paired_bca_top1_any']['point_estimate']}, "
              f"sig={deciles[f'D{d+1}']['paired_bca_top1_any']['ci_excludes_zero']}", flush=True)

    print("Computing similarity-spread diagnostic (overall + within lowest-density deciles)...", flush=True)
    tenn_ranks = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = H4.bestsim_ranks(records)
    rows = []
    for i, r in enumerate(records):
        sg = spread_gap_of(r)
        if sg is None:
            continue
        outcome = rr(tenn_ranks[i]["rank_any"]) - rr(bs_ranks[i]["rank_any"])
        rows.append({"idx": i, "density": density_of(r), "spread_gap": sg, "outcome": outcome})

    # find whichever deciles are still in the "valley" region (any-annotated Top1 delta < 0)
    valley_deciles = [d for d, rep in deciles.items() if rep["paired_bca_top1_any"]["point_estimate"] < 0]
    valley_idx = set()
    for d in valley_deciles:
        dnum = int(d[1:]) - 1
        lo = dnum * bin_size
        hi = (dnum + 1) * bin_size if dnum < N_DECILES - 1 else n
        valley_idx.update(i for i, r in by_density[lo:hi])

    valley_rows = [row for row in rows if row["idx"] in valley_idx]

    def group_compare(subset, label):
        if len(subset) < 10:
            return {"note": f"n={len(subset)} too small", "label": label}
        s_sorted = sorted(subset, key=lambda r: r["spread_gap"])
        half = len(s_sorted) // 2
        narrow, wide = s_sorted[:half], s_sorted[half:]
        return {
            "label": label, "n": len(subset),
            "narrow_half": {"n": len(narrow), "mean_spread": round(float(np.mean([r["spread_gap"] for r in narrow])), 3),
                             "mean_outcome": round(float(np.mean([r["outcome"] for r in narrow])), 4)},
            "wide_half": {"n": len(wide), "mean_spread": round(float(np.mean([r["spread_gap"] for r in wide])), 3),
                          "mean_outcome": round(float(np.mean([r["outcome"] for r in wide])), 4)},
        }

    spread_report = {
        "valley_deciles_identified": valley_deciles,
        "overall": group_compare(rows, "all queries"),
        "within_valley_deciles": group_compare(valley_rows, "valley deciles only"),
    }

    print("Computing within-decile target-depth predictiveness (joint check, non-linear via stratification)...", flush=True)
    depth_by_decile = {}
    for d in range(N_DECILES):
        lo = d * bin_size
        hi = (d + 1) * bin_size if d < N_DECILES - 1 else n
        bucket_idx = [i for i, r in by_density[lo:hi]]
        bucket_records = [records[i] for i in bucket_idx]
        ranked_bs = [sorted(r["bestsim_results"], key=lambda ts: -ts[1]) for r in bucket_records]
        recovered, not_recovered = [], []
        for r, ranked in zip(bucket_records, ranked_bs):
            top10 = set(t for t, s in ranked[:10])
            annotated = set(r["annotated_targets"])
            depths = [depth_by_target.get(t, 0) for t in annotated]
            mean_depth = float(np.mean(depths)) if depths else None
            if mean_depth is None:
                continue
            (recovered if (annotated & top10) else not_recovered).append(mean_depth)
        depth_by_decile[f"D{d+1}"] = {
            "n_recovered": len(recovered), "n_not_recovered": len(not_recovered),
            "mean_depth_recovered": round(float(np.mean(recovered)), 1) if recovered else None,
            "mean_depth_not_recovered": round(float(np.mean(not_recovered)), 1) if not_recovered else None,
            "gap": (round(float(np.mean(recovered)) - float(np.mean(not_recovered)), 1)
                    if recovered and not_recovered else None),
        }

    out = {
        "set": args.set, "n_queries": n,
        "decile_dose_response": deciles,
        "spread_diagnostic": spread_report,
        "within_decile_depth_predictiveness": depth_by_decile,
    }
    out_path = os.path.join(RESULTS_DIR, f"phase1_density_confirmatory_{args.set}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
