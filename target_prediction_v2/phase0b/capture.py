"""
Phase 0b (target_prediction_v2_procedure_rev5.md, Section 5) -- the shared
per-query capture backbone. Rev 5's Block 1/2 tasks (P0b-1, 3, 4, 5, 6, 7, 8,
9, 11, 14, 15, 16) all need the SAME expensive per-query computation (leakage
removal + Tanimoto against the frozen v1 reference), just re-scored or
re-aggregated differently. Phase 0 (rev3) ran a separate full search per
experiment; this version computes each query's raw neighbourhood ONCE and
saves enough detail (per-neighbour tanimoto + per-neighbour target/pchembl
annotations, up to NEIGHBOR_CAP compounds) that every downstream analysis --
best_similarity ranking, unweighted 10-NN, weighted k-NN at any (k, alpha),
potency-weighted variants, nested-CV grid search, bootstrap resampling --
is a fast in-memory re-score, not a re-query.

Re-uses phase0/leakage.py and phase0/metrics.py unmodified (same frozen v1
reference, same leakage-control logic) -- no duplication of that logic here.

Output: results/capture_<setname>.jsonl, one JSON object per query:
  {
    "set": "A"|"B"|"D", "smiles": ..., "molecule_chembl_id": ...,
    "annotated_targets": [...], "mechanistic_targets": [...] (see note),
    "primary_targets": [...],
    "max_similarity_to_index": float,   -- top hit's own tanimoto (P0b's
                                            "standard stratification")
    "bestsim_results": [[target_chembl, best_similarity], ...]  -- ALL
        targets with best_similarity >= 0.15 (bounds storage; every sweep
        threshold in P0b-14, 0.2-0.6, is >= 0.15 so nothing needed by the
        sweep is truncated)
    "neighbours": [[smiles, tanimoto, [[target_chembl, pchembl_or_null], ...]], ...]
        -- top NEIGHBOR_CAP distinct reference COMPOUNDS by tanimoto
        (post-leakage-removal), each with its full per-target annotation
        list restricted to the reduced reference. This is what lets any
        (k, alpha, potency-weight-function) be re-scored later without
        re-running the search.
  }

NOTE on "mechanistic_targets": eval_sets.json's mechanistic_targets and
primary_targets fields are IDENTICAL for every query (verified: all 6,984
rows in mechanisms_raw.csv have molecular_mechanism=1 -- ChEMBL's own
mechanism.json endpoint apparently only returns validated
molecular-mechanism rows in bulk, so "any mechanism.json row" and
"molecular_mechanism=1" never actually differ in our pull). So
eval_sets.json's "mechanistic_targets" is really Rev 5's PRIMARY/INTENDED
tier (drug_mechanism presence, §4.1), not the broader two-of-three
"mechanistically supported" tier Rev 5 defines -- that broader tier needs
per-document data (see mechanism_support.py) and is computed separately,
not carried in this capture file.

Usage: python3 capture.py --set A --n 200 --seed 42
       python3 capture.py --set B --n 200 --seed 42
       python3 capture.py --set D --n -1 --seed 42   (n=-1 = all)
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index"))

import target_fishing as TF  # noqa: E402
import leakage as LK  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

NEIGHBOR_CAP = 300          # widened from 150 -- rev6 rerun found k hitting the old grid's
                            # edge (100) in several folds, meaning the old cap of 150 stored
                            # neighbours wasn't leaving headroom for a genuinely wider k grid
BESTSIM_MIN_STORE = 0.15    # below the lowest P0b-14 sweep point (0.2) -- safety margin
LEAKAGE_CUTOFF = 0.95       # primary cutoff; leakage sensitivity handled separately (see leakage_sensitivity.py)


def load_index():
    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE0_DIR, "data", "leakage_fp2.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])
    with open(os.path.join(PHASE0_DIR, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)
    return leak_idx, eval_sets


def capture_one(row, leak_idx, cutoff=LEAKAGE_CUTOFF):
    smi = row["smiles"]
    reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, cutoff)
    q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)

    results = TF._search_against(reduced_fps, reduced_pop, reduced_df, q_packed, q_pop,
                                  threshold=0.0, max_compounds_per_target=1)
    max_sim = max((r["best_similarity"] for r in results), default=0.0)
    bestsim_results = [[r["target_chembl"], r["best_similarity"]] for r in results
                        if r["best_similarity"] >= BESTSIM_MIN_STORE]

    tanimoto = LK.tanimoto_to_all(q_packed, reduced_fps, reduced_pop, q_pop)
    tmp = reduced_df.copy()
    tmp["tanimoto"] = tanimoto
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    top_n = by_compound.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)

    neighbours = []
    if not top_n.empty:
        neighbor_smiles = set(top_n.index)
        neighbor_rows = tmp[tmp["smiles"].isin(neighbor_smiles)][["smiles", "target_chembl", "pchembl_value"]]
        by_smi = {}
        for r in neighbor_rows.itertuples(index=False):
            pv = None if r.pchembl_value != r.pchembl_value else float(r.pchembl_value)  # NaN-safe
            by_smi.setdefault(r.smiles, []).append([r.target_chembl, pv])
        for nsmi, sim in top_n["tanimoto"].items():
            neighbours.append([nsmi, round(float(sim), 4), by_smi.get(nsmi, [])])

    return {
        "smiles": smi,
        "molecule_chembl_id": row["molecule_chembl_id"],
        "max_phase": row["max_phase"],
        "annotated_targets": row["annotated_targets"],
        "mechanistic_targets": row["mechanistic_targets"],
        "primary_targets": row["primary_targets"],
        "max_similarity_to_index": round(float(max_sim), 4),
        "bestsim_results": bestsim_results,
        "neighbours": neighbours,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["A", "B", "D"])
    ap.add_argument("--n", type=int, default=200, help="-1 = use all matched queries")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    leak_idx, eval_sets = load_index()
    key = {"A": "set_A_approved", "B": "set_B_clinical", "D": "set_D_sparse_target"}[args.set]
    rows = list(eval_sets[key])

    rng = random.Random(args.seed)
    rng.shuffle(rows)
    if args.n > 0:
        rows = rows[:args.n]

    out_path = os.path.join(RESULTS_DIR, f"capture_{args.set}.jsonl")
    print(f"Capturing Set {args.set}: {len(rows)} queries (seed={args.seed}) -> {out_path}", flush=True)
    t0 = time.time()
    with open(out_path, "w") as f:
        for i, row in enumerate(rows):
            rec = capture_one(row, leak_idx)
            rec["set"] = args.set
            f.write(json.dumps(rec) + "\n")
            f.flush()
            if (i + 1) % 25 == 0 or (i + 1) == len(rows):
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                remaining = (len(rows) - (i + 1)) / rate if rate > 0 else 0
                print(f"  [{args.set}] {i+1}/{len(rows)} ({elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining)", flush=True)

    print(f"Done: {out_path} ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
