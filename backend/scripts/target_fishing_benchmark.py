"""
Scaffold-split validation benchmark for target_fishing.py.

The question this answers: when the reference pool has NEVER seen the query
compound, or anything sharing its Murcko scaffold, can target-fishing still
recover a target independently known to interact with that compound? This
is the only benchmark design that actually tests scaffold-hopping
generalization rather than database lookup — a random compound split would
let near-identical analogs leak between "reference" and "test," inflating
apparent performance (see module docstring precedent: target_fishing.py's
whole 2026 redesign already treats "is this really independent of what we
already had" as the central question).

METHODOLOGY (recorded here AND in the report's own "methodology" block —
this needs to survive independent of anyone remembering this file):

1. Scaffold groups: every DISTINCT compound in the index (857,232 as of
   the 2026 ChEMBL pull) is grouped by its precomputed Murcko scaffold.
   A compound with no ring system (no scaffold — RDKit gives it an empty
   string, ~0.36% of the pool) is treated as its own singleton scaffold
   group (its own SMILES as the group key) — it shares no real structural
   core with anything else, so it can't leak via scaffold-sharing either.

2. Eligible holdout scaffolds: only scaffold groups with 1-50 members.
   This is NOT an arbitrary cutoff — the real group-size distribution has
   catastrophic skew (median group size is 1, but the single largest is
   *benzene*, c1ccccc1, with 14,381 members — 1.6% of the entire pool —
   and several other generic fragments have 1,000+). Holding out a
   mega-scaffold would (a) gut the reference pool for every target that
   happens to have a benzene-containing analog, and (b) produce a test
   set with no coherent "one chemical series" identity, defeating the
   point of a SCAFFOLD split specifically. Verified empirically before
   writing this script (see conversation record), not assumed.

3. Random sample (fixed seed for reproducibility) of eligible scaffolds,
   accumulated until reaching TARGET_N_HELDOUT_COMPOUNDS held-out
   compounds (a whole scaffold group is always taken together, never
   split — that's the actual leakage boundary).

4. Leakage removal: every row (across ALL targets, not just the one being
   evaluated) whose compound is in a held-out scaffold group is removed
   from the reference index used for searching. Ground truth itself still
   comes from the ORIGINAL, pre-removal data.

5. Ground truth: for each held-out compound, every DISTINCT target it was
   originally associated with is a separate evaluation instance (a
   5-target compound = 5 chances to recover, each scored independently —
   not collapsed into one pass/fail). Capped at MAX_TARGETS_PER_COMPOUND
   per compound (randomly sampled if it has more) so one extremely
   promiscuous compound (the real max in this pool is 1,810 distinct
   targets for a single compound) can't dominate the aggregate stats.

6. Search: each held-out compound is queried (threshold=0.0, so a true
   target's rank is always defined, never silently excluded by a
   threshold) against the REDUCED reference via target_fishing.py's own
   _search_against() — the exact production aggregation/ranking code,
   not a reimplementation.

7. Two rankings are scored for comparison: by evidence_score (what real
   users see) and by best_similarity alone (a naive-similarity baseline)
   — this directly answers whether the evidence_score formula improves
   ranking quality over raw Tanimoto, or is just cosmetic.

8. Near-exact split: if a held-out compound's own top hit (after removing
   its whole scaffold group) still comes back at best_similarity >= 0.99,
   that compound's instances are routed to the SECONDARY report bucket
   (near-exact retrieval) instead of PRIMARY (novel-scaffold
   generalization) — a folded 2048-bit Morgan fingerprint can collide for
   distinct structures, and a stereoisomer/tautomer/salt-form variant can
   get a technically-different Murcko scaffold despite being chemically
   near-identical to something still in the reference; lumping those into
   the primary "can it generalize to a new scaffold" number would inflate
   it with what's still effectively lookup.

9. Metrics (primary bucket is the headline; secondary reported alongside
   for comparison, never merged into it): Top-1/5/10/20 recovery, Mean
   Reciprocal Rank, random-baseline expected recovery (K / n_targets in
   the reduced reference) and enrichment (observed / random), each
   computed for both rankings (evidence_score, best_similarity).
   Stratified by (a) the target's ORIGINAL reference depth (how many
   distinct compounds supported it before any removal: 1-2, 3-10, 11-50,
   >50) and (b) the achieved similarity band for recovered instances
   (0.0-0.5, 0.5-0.6, ..., 0.9-0.99 — 0.99+ is the secondary bucket).

Usage (run from the repo root):
    python -m scripts.target_fishing_benchmark [--n-heldout 5000] [--seed 42] \\
        [--out ligand_audit... ] (see --out default below)
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RDLogger.DisableLog("rdApp.*")

import target_fishing as TF

MIN_SCAFFOLD_GROUP_SIZE = 1
MAX_SCAFFOLD_GROUP_SIZE = 50    # see module docstring point 2 — empirically justified, not a guess
MAX_TARGETS_PER_COMPOUND = 20   # see point 5 — caps one promiscuous compound's influence
NEAR_EXACT_THRESHOLD = 0.99     # see point 8
TOP_KS = (1, 5, 10, 20)
REF_DEPTH_BUCKETS = [(1, 2), (3, 10), (11, 50), (51, float("inf"))]
SIM_BUCKETS = [(0.0, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.99)]


def _scaffold_key(smiles, scaffold):
    """A compound with no Murcko scaffold (no ring system) gets its own
       SMILES as a synthetic singleton scaffold key — see docstring
       point 1. NaN-safe: pandas gives NaN for an empty/missing scaffold
       string on CSV round-trip."""
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring_system__:{smiles}"


def build_scaffold_groups(compound_df):
    compound_df = compound_df.copy()
    compound_df["scaffold_key"] = [
        _scaffold_key(s, sc) for s, sc in zip(compound_df["smiles"], compound_df["murcko_scaffold"])
    ]
    groups = compound_df.groupby("scaffold_key")["smiles"].apply(list)
    return groups  # Series: scaffold_key -> [smiles, ...]


def select_heldout_scaffolds(groups, target_n_compounds, seed):
    sizes = groups.apply(len)
    eligible = sizes[(sizes >= MIN_SCAFFOLD_GROUP_SIZE) & (sizes <= MAX_SCAFFOLD_GROUP_SIZE)]
    keys = list(eligible.index)
    rng = random.Random(seed)
    rng.shuffle(keys)

    held_out_keys = []
    held_out_smiles = set()
    for k in keys:
        held_out_keys.append(k)
        held_out_smiles.update(groups[k])
        if len(held_out_smiles) >= target_n_compounds:
            break
    return held_out_keys, held_out_smiles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-heldout", type=int, default=5000, help="Target number of held-out compounds")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--index-dir", default=os.environ.get("TARGET_FISHING_INDEX_DIR",
                     os.path.join(os.path.dirname(__file__), "..", "target_fishing_index")))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "..", "target_fishing_benchmark_report.json"))
    args = ap.parse_args()

    t0 = time.time()
    os.environ["TARGET_FISHING_INDEX_DIR"] = args.index_dir
    print("Loading full index...", flush=True)
    full_fps, full_pop_counts, full_df = TF._load()
    n_total_targets_full = int(full_df["target_chembl"].nunique())
    print(f"  {len(full_df)} compound-target pairs, {full_df['smiles'].nunique()} compounds, "
          f"{n_total_targets_full} targets ({time.time()-t0:.0f}s)", flush=True)

    print("Building scaffold groups...", flush=True)
    compound_df = full_df.drop_duplicates(subset="smiles")[["smiles", "murcko_scaffold"]].reset_index(drop=True)
    groups = build_scaffold_groups(compound_df)
    sizes = groups.apply(len)
    print(f"  {len(groups)} distinct scaffold groups (incl. singleton no-ring-system compounds); "
          f"size range {sizes.min()}-{sizes.max()}, median {int(sizes.median())} ({time.time()-t0:.0f}s)", flush=True)

    held_out_keys, held_out_smiles = select_heldout_scaffolds(groups, args.n_heldout, args.seed)
    print(f"  selected {len(held_out_keys)} scaffold groups -> {len(held_out_smiles)} held-out compounds "
          f"(target was {args.n_heldout})", flush=True)

    # ---- ground truth: from the ORIGINAL, pre-removal data ----
    print("Building ground truth + original reference-depth-per-target (pre-removal)...", flush=True)
    original_depth_by_target = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size()
    heldout_mask = full_df["smiles"].isin(held_out_smiles).values
    ground_truth_df = full_df.loc[heldout_mask, ["smiles", "target_chembl"]].drop_duplicates()
    rng = random.Random(args.seed)
    ground_truth = {}   # smiles -> [target_chembl, ...] (capped)
    for smi, sub in ground_truth_df.groupby("smiles"):
        targets = sub["target_chembl"].tolist()
        if len(targets) > MAX_TARGETS_PER_COMPOUND:
            targets = rng.sample(targets, MAX_TARGETS_PER_COMPOUND)
        ground_truth[smi] = targets
    n_eval_instances = sum(len(v) for v in ground_truth.values())
    print(f"  {len(ground_truth)} held-out compounds with >=1 known target, "
          f"{n_eval_instances} total (compound, true-target) evaluation instances", flush=True)

    # ---- build the leakage-free reduced reference ----
    reference_mask = ~heldout_mask
    reference_df = full_df.loc[reference_mask].reset_index(drop=True)
    reference_fps = full_fps[reference_mask]
    reference_pop_counts = full_pop_counts[reference_mask]
    n_targets_in_reference = int(reference_df["target_chembl"].nunique())
    print(f"  reduced reference: {len(reference_df)} pairs, {reference_df['smiles'].nunique()} compounds, "
          f"{n_targets_in_reference} targets ({time.time()-t0:.0f}s)", flush=True)

    # fingerprint lookup for held-out compounds, keyed by smiles (need
    # ONE representative row's fingerprint per compound from the FULL
    # arrays, since a held-out compound may have appeared multiple times
    # -- once per target it was originally tested against)
    first_row_per_smiles = full_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]

    # ---- run the benchmark ----
    print(f"\nRunning {len(ground_truth)} queries against the reduced reference...", flush=True)
    instances = []
    t_search0 = time.time()
    for i, (smi, true_targets) in enumerate(ground_truth.items()):
        row_idx = int(first_row_per_smiles[smi])
        q_packed = full_fps[row_idx]
        q_pop = int(full_pop_counts[row_idx])
        results = TF._search_against(reference_fps, reference_pop_counts, reference_df, q_packed, q_pop,
                                      threshold=0.0, max_compounds_per_target=1)
        by_evidence = {r["target_chembl"]: (rank + 1, r["best_similarity"]) for rank, r in enumerate(results)}
        by_similarity_order = sorted(results, key=lambda r: -r["best_similarity"])
        by_similarity = {r["target_chembl"]: (rank + 1) for rank, r in enumerate(by_similarity_order)}
        near_exact = bool(results) and results[0]["best_similarity"] >= NEAR_EXACT_THRESHOLD

        for true_target in true_targets:
            rank_evidence, achieved_sim = by_evidence.get(true_target, (None, 0.0))
            rank_similarity = by_similarity.get(true_target, None)
            instances.append({
                "smiles": smi,
                "true_target": true_target,
                "rank_evidence": rank_evidence,
                "rank_similarity": rank_similarity,
                "achieved_similarity": achieved_sim,
                "original_reference_depth": int(original_depth_by_target.get(true_target, 0)),
                "near_exact": near_exact,
            })

        if (i + 1) % 500 == 0 or (i + 1) == len(ground_truth):
            elapsed = time.time() - t_search0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(ground_truth) - (i + 1)) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(ground_truth)}] {elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining", flush=True)

    print(f"\nSearch phase done in {time.time()-t_search0:.0f}s. Computing metrics...", flush=True)

    # ---- metrics ----
    def sim_bucket(s):
        for lo, hi in SIM_BUCKETS:
            if lo <= s < hi:
                return f"{lo:.1f}-{hi:.2f}"
        return None

    def depth_bucket(d):
        for lo, hi in REF_DEPTH_BUCKETS:
            if lo <= d <= hi:
                return f"{lo}-{hi if hi != float('inf') else '+'}"
        return None

    def compute_metrics(rows, rank_key):
        n = len(rows)
        if n == 0:
            return {"n_instances": 0}
        ranks = [r[rank_key] for r in rows]
        out = {"n_instances": n}
        for k in TOP_KS:
            n_recovered = sum(1 for rk in ranks if rk is not None and rk <= k)
            recovery = n_recovered / n
            random_baseline = min(k / n_targets_in_reference, 1.0)
            out[f"top_{k}_recovery"] = round(recovery, 4)
            out[f"top_{k}_random_baseline"] = round(random_baseline, 4)
            out[f"top_{k}_enrichment"] = round(recovery / random_baseline, 2) if random_baseline > 0 else None
        mrr = sum((1.0 / rk) if rk is not None else 0.0 for rk in ranks) / n
        out["mrr"] = round(mrr, 4)
        return out

    primary = [r for r in instances if not r["near_exact"]]
    secondary = [r for r in instances if r["near_exact"]]

    report = {
        "methodology": {
            "n_indexed_compound_target_pairs_full": len(full_df),
            "n_indexed_compounds_full": int(full_df["smiles"].nunique()),
            "n_indexed_targets_full": n_total_targets_full,
            "n_distinct_scaffold_groups": len(groups),
            "scaffold_group_size_range_eligible_for_holdout": [MIN_SCAFFOLD_GROUP_SIZE, MAX_SCAFFOLD_GROUP_SIZE],
            "n_scaffold_groups_held_out": len(held_out_keys),
            "n_heldout_compounds": len(held_out_smiles),
            "n_heldout_compounds_with_ground_truth": len(ground_truth),
            "n_evaluation_instances": len(instances),
            "n_evaluation_instances_primary_novel_scaffold": len(primary),
            "n_evaluation_instances_secondary_near_exact": len(secondary),
            "max_targets_per_compound_cap": MAX_TARGETS_PER_COMPOUND,
            "near_exact_similarity_threshold": NEAR_EXACT_THRESHOLD,
            "n_targets_in_reduced_reference": n_targets_in_reference,
            "search_threshold": 0.0,
            "seed": args.seed,
            "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
        "primary_novel_scaffold": {
            "by_evidence_score": compute_metrics(primary, "rank_evidence"),
            "by_best_similarity": compute_metrics(primary, "rank_similarity"),
        },
        "secondary_near_exact_retrieval": {
            "by_evidence_score": compute_metrics(secondary, "rank_evidence"),
            "by_best_similarity": compute_metrics(secondary, "rank_similarity"),
        },
        "stratification_by_reference_depth_primary_only": {},
        "stratification_by_similarity_band_primary_only": {},
    }

    for lo, hi in REF_DEPTH_BUCKETS:
        label = f"{lo}-{hi if hi != float('inf') else '+'}"
        subset = [r for r in primary if lo <= r["original_reference_depth"] <= hi]
        report["stratification_by_reference_depth_primary_only"][label] = compute_metrics(subset, "rank_evidence")

    for lo, hi in SIM_BUCKETS:
        label = f"{lo:.1f}-{hi:.2f}"
        subset = [r for r in primary if lo <= r["achieved_similarity"] < hi]
        report["stratification_by_similarity_band_primary_only"][label] = compute_metrics(subset, "rank_evidence")

    out_path = os.path.abspath(args.out)
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    instances_path = out_path.replace(".json", "_instances.csv")
    pd.DataFrame(instances).to_csv(instances_path, index=False)

    print(f"\nDone in {time.time()-t0:.0f}s total.")
    print(f"Report: {out_path}")
    print(f"Raw per-instance results: {instances_path}")
    print(json.dumps(report["methodology"], indent=2))
    print("\nPRIMARY (novel scaffold) — by evidence_score:")
    print(json.dumps(report["primary_novel_scaffold"]["by_evidence_score"], indent=2))
    print("\nPRIMARY (novel scaffold) — by best_similarity alone:")
    print(json.dumps(report["primary_novel_scaffold"]["by_best_similarity"], indent=2))
    print("\nSECONDARY (near-exact retrieval) — by evidence_score:")
    print(json.dumps(report["secondary_near_exact_retrieval"]["by_evidence_score"], indent=2))


if __name__ == "__main__":
    main()
