"""
THE LOCKED SCAFFOLD TEST -- opened once, per explicit user authorization,
after the fresh-holdout sanity check (test_fresh_holdout_validation.py)
showed real positive signal. Sealed since Phase 1 (phase1/build_holdouts.py,
795 compounds / 598 scaffold groups).

CRITICAL, disclosed finding checked BEFORE running this: 476 of 795 locked
compounds (60%) were already used somewhere in prior Phase 3 tuning
(capture_scaffold_strict_v2_{A,B} or the fresh holdout) -- the tuning
scripts sampled directly from eval_sets.json without excluding the locked
test's scaffold groups. This capture pulls all resolvable rows (909, some
duplicate smiles across Set A/B); the scoring script (test_locked_test_g2.py)
reports the CLEAN 433-row never-before-touched subset as the real,
unbiased G2 answer, with the full/contaminated set shown only for context.

Same rigorous methodology as every other capture this session (near-dup
leakage removal under two independent fingerprints, then whole-scaffold-
group removal) -- not relaxed for this one.

Usage: python3 capture_locked_test.py
"""
import json
import os
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


def load_index():
    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE2_DIR, "data", "v2_leakage_fp2_filtered.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])
    with open(os.path.join(PHASE2_DIR, "data", "stage6_std_cache.json")) as f:
        std_cache = json.load(f)
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))
    scaffold_key_by_smiles = {s: scaffold_key(s, sc) for s, sc in scaffold_by_smiles.items()}
    return leak_idx, full_df, scaffold_key_by_smiles, std_cache


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
    os.makedirs(RESULTS_DIR, exist_ok=True)
    leak_idx, full_df, scaffold_key_by_smiles, std_cache = load_index()

    with open(os.path.join(RESULTS_DIR, "locked_test_rows.json")) as f:
        locked_data = json.load(f)
    all_rows = locked_data["all"]
    clean_ids = set(locked_data["clean_molecule_chembl_ids"])

    v2_index_smiles = set(full_df["smiles"].unique())
    rows = []
    n_dropped = 0
    for r in all_rows:
        v2_smi = std_cache.get(r["molecule_chembl_id"])
        if not v2_smi or v2_smi not in v2_index_smiles:
            n_dropped += 1
            continue
        r2 = dict(r)
        r2["smiles"] = v2_smi
        r2["is_clean"] = r["molecule_chembl_id"] in clean_ids
        rows.append(r2)
    print(f"{len(all_rows)} locked-test rows, {n_dropped} not in v2 index, "
          f"{len(rows)} to capture ({sum(1 for r in rows if r['is_clean'])} clean)", flush=True)

    out_path = os.path.join(RESULTS_DIR, "capture_locked_test.jsonl")
    print(f"Capturing THE LOCKED SCAFFOLD TEST: {len(rows)} queries -> {out_path}", flush=True)
    t0 = time.time()
    with open(out_path, "w") as f:
        for i, row in enumerate(rows):
            rec = capture_one(row, leak_idx, full_df, scaffold_key_by_smiles)
            rec["is_clean"] = row["is_clean"]
            f.write(json.dumps(rec) + "\n")
            f.flush()
            if (i + 1) % 25 == 0 or (i + 1) == len(rows):
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                remaining = (len(rows) - (i + 1)) / rate if rate > 0 else 0
                print(f"  {i+1}/{len(rows)} ({elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining)", flush=True)

    print(f"Done: {out_path} ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
