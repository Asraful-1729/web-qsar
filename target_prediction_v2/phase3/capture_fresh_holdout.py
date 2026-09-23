"""
Fresh-sample sanity check, Phase A step 1 of making v2 usable: every prior
Phase 3 test (potency weighting, orthologue ablation, stacker, calibration)
reused the SAME n=1600 queries (capture_scaffold_strict_v2_{A,B}.jsonl)
repeatedly across multiple rounds of tuning -- a real risk of overfitting
to that specific set that has never been checked. This draws a genuinely
NEW sample (different seed, explicit exclusion of every molecule_chembl_id
already used) with the exact same capture methodology (near-dup leakage
removal under two independent fingerprints, then whole-scaffold-group
removal), so the resulting numbers are comparable to everything already
measured.

This is NOT the locked scaffold test (that stays sealed) -- a cheap,
reversible intermediate check, run before deciding whether opening the
real locked test is even worth doing.

Usage: python3 capture_fresh_holdout.py --set A --n 400 --seed 123
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
PHASE2_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase2")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(os.path.join(PHASE2_DIR, "v2_index"))

import target_fishing as TF  # noqa: E402
import leakage as LK  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

NEIGHBOR_CAP = 300
BESTSIM_MIN_STORE = 0.15
NEAR_DUP_CUTOFF = 0.95


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def load_already_used_ids():
    used = set()
    for letter in ("A", "B"):
        path = os.path.join(RESULTS_DIR, f"capture_scaffold_strict_v2_{letter}.jsonl")
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                used.add(r["molecule_chembl_id"])
    return used


def load_index():
    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE2_DIR, "data", "v2_leakage_fp2_filtered.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])
    with open(os.path.join(PHASE0_DIR, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)
    with open(os.path.join(PHASE2_DIR, "data", "stage6_std_cache.json")) as f:
        std_cache = json.load(f)
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))
    scaffold_key_by_smiles = {s: scaffold_key(s, sc) for s, sc in scaffold_by_smiles.items()}
    return leak_idx, eval_sets, full_df, scaffold_key_by_smiles, std_cache


def capture_one(row, leak_idx, full_df, scaffold_key_by_smiles, cutoff=NEAR_DUP_CUTOFF):
    smi = row["smiles"]
    reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, cutoff)
    q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)

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
    ap.add_argument("--set", required=True, choices=["A", "B"])
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    already_used = load_already_used_ids()
    print(f"{len(already_used)} molecule_chembl_ids already used across all prior tuning -- excluding", flush=True)

    leak_idx, eval_sets, full_df, scaffold_key_by_smiles, std_cache = load_index()
    key = {"A": "set_A_approved", "B": "set_B_clinical"}[args.set]
    all_rows = list(eval_sets[key])

    v2_index_smiles = set(full_df["smiles"].unique())
    rows = []
    n_dropped_used = 0
    n_dropped_no_v2 = 0
    for r in all_rows:
        if r["molecule_chembl_id"] in already_used:
            n_dropped_used += 1
            continue
        v2_smi = std_cache.get(r["molecule_chembl_id"])
        if not v2_smi or v2_smi not in v2_index_smiles:
            n_dropped_no_v2 += 1
            continue
        r2 = dict(r)
        r2["smiles"] = v2_smi
        rows.append(r2)
    print(f"Set {args.set}: {len(all_rows)} total, {n_dropped_used} already used (excluded), "
          f"{n_dropped_no_v2} not in v2 index, {len(rows)} eligible fresh candidates", flush=True)

    rng = random.Random(args.seed)
    rng.shuffle(rows)
    if args.n > 0:
        rows = rows[:args.n]

    out_path = os.path.join(RESULTS_DIR, f"capture_fresh_holdout_{args.set}.jsonl")
    print(f"Capturing FRESH Set {args.set} holdout: {len(rows)} queries "
          f"(seed={args.seed}, never used in any prior tuning) -> {out_path}", flush=True)
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
