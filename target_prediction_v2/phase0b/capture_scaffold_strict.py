"""
Rev 6 flag-3 (blocking): the original H4 decomposition (capture.py ->
analyze_h4.py) ran under NEAR-DUPLICATE-ONLY leakage control (Tanimoto
>=0.95 under two fingerprints), the same control Phase 0 always used for
compound-level searches. Separately, scaffold_leakage.py showed that
removing the query's WHOLE Bemis-Murcko scaffold group (not just near-
duplicates) drops recovery by a further 4-7 points on Set A. Those two
facts were never combined: the H4a/H4b comparison has not yet been run
under the stricter control, so it isn't known whether "pooling beats
best-similarity" (H4a) or "weighting doesn't beat unweighted pooling"
(H4b) survive it.

This script is capture.py's neighbour-capture logic, with ONE change: after
the standard near-duplicate leakage removal, ALSO remove every reference
compound sharing the query's own Murcko scaffold (same logic as
scaffold_leakage.py's reference reduction), before computing best_similarity
results and the neighbour pool. Everything downstream (score.py, bca.py,
the nested-CV/BCa machinery in analyze_h4.py) is reused unchanged -- this
produces a second capture file in the exact same schema as capture.py's, so
analyze_h4.py can be pointed at it with no code changes, only a different
input path.

Usage: python3 capture_scaffold_strict.py --set A --n 200 --seed 42
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

NEIGHBOR_CAP = 300          # widened from 150 -- matches capture.py's widening, see its comment
BESTSIM_MIN_STORE = 0.15
NEAR_DUP_CUTOFF = 0.95


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def load_index():
    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE0_DIR, "data", "leakage_fp2.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])
    with open(os.path.join(PHASE0_DIR, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))
    scaffold_key_by_smiles = {s: scaffold_key(s, sc) for s, sc in scaffold_by_smiles.items()}
    return leak_idx, eval_sets, full_df, scaffold_key_by_smiles


def capture_one(row, leak_idx, full_df, scaffold_key_by_smiles, cutoff=NEAR_DUP_CUTOFF):
    smi = row["smiles"]
    # Step 1: standard near-duplicate reduction (same as capture.py)
    reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, cutoff)
    q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)

    # Step 2: ADDITIONALLY strip every compound sharing the query's own scaffold
    q_scaffold = scaffold_key_by_smiles.get(smi, f"__no_ring__:{smi}")
    same_scaffold_mask = reduced_df["smiles"].map(lambda s: scaffold_key_by_smiles.get(s, f"__no_ring__:{s}")) == q_scaffold
    keep_mask = ~same_scaffold_mask.values
    scaffold_group_size = int(same_scaffold_mask.sum())

    sc_fps = reduced_fps[keep_mask]
    sc_pop = reduced_pop[keep_mask]
    sc_df = reduced_df.loc[keep_mask]

    results = TF._search_against(sc_fps, sc_pop, sc_df, q_packed, q_pop,
                                  threshold=0.0, max_compounds_per_target=1)
    max_sim = max((r["best_similarity"] for r in results), default=0.0)
    bestsim_results = [[r["target_chembl"], r["best_similarity"]] for r in results
                        if r["best_similarity"] >= BESTSIM_MIN_STORE]

    tanimoto = LK.tanimoto_to_all(q_packed, sc_fps, sc_pop, q_pop)
    tmp = sc_df.copy()
    tmp["tanimoto"] = tanimoto
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    top_n = by_compound.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)

    neighbours = []
    if not top_n.empty:
        neighbor_smiles = set(top_n.index)
        neighbor_rows = tmp[tmp["smiles"].isin(neighbor_smiles)][["smiles", "target_chembl", "pchembl_value"]]
        by_smi = {}
        for r in neighbor_rows.itertuples(index=False):
            pv = None if r.pchembl_value != r.pchembl_value else float(r.pchembl_value)
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
        "scaffold_group_size_removed": scaffold_group_size,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["A", "B", "D"])
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    leak_idx, eval_sets, full_df, scaffold_key_by_smiles = load_index()
    key = {"A": "set_A_approved", "B": "set_B_clinical", "D": "set_D_sparse_target"}[args.set]
    rows = list(eval_sets[key])

    rng = random.Random(args.seed)
    rng.shuffle(rows)
    if args.n > 0:
        rows = rows[:args.n]

    out_path = os.path.join(RESULTS_DIR, f"capture_scaffold_strict_{args.set}.jsonl")
    print(f"Capturing Set {args.set} (scaffold-strict): {len(rows)} queries (seed={args.seed}) -> {out_path}", flush=True)
    t0 = time.time()
    with open(out_path, "w") as f:
        for i, row in enumerate(rows):
            rec = capture_one(row, leak_idx, full_df, scaffold_key_by_smiles)
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
