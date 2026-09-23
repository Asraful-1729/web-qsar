"""
Phase 3A: re-fit prerequisite. This is phase0b/capture_scaffold_strict.py's
exact logic (near-duplicate leakage removal under two independent
fingerprints, THEN whole-scaffold-group removal, then neighbour capture),
pointed at the Phase 2 rebuilt index instead of v1's frozen one --
required because Phase 2 item 3's Potency-record inclusion roughly
triples native-human pair volume (1,312,849 -> 3,973,438), which
materially shifts the density distribution the original density-adaptive
rule (density>=12) was fit against (PHASE1_DENSITY_ADAPTIVE_RULE.md).
That threshold cannot be reused unchanged -- BUILD_PLAN.md Phase 3A's own
text requires a re-fit before using it for anything.

Query smiles resolution: eval_sets.json's own `smiles` field was matched
against v1's (buggy) SaltRemover-based standardization. Phase 2 uses the
fixed rdMolStandardize.FragmentParent() (item 5's fix) instead, so a
query drug's standardized SMILES can differ from the old index's version.
Confirmed empirically before writing this: 88.9%/93.4% of Set A/B query
drugs are present in Phase 2's stage6_std_cache.json at all (the rest
simply have no activity data in the Phase 2 pull's scope); of those
present, 98%+ have IDENTICAL smiles to the old value (the fix mostly
preserves existing matches). This script resolves each query's smiles
via stage6_std_cache.json (the authoritative current mapping), not
eval_sets.json's stale field, and skips any drug absent from the v2
index -- a real, disclosed coverage change from n=1965/2750 to ~88.9%/
93.4% of that, not silently patched over.

Usage: python3 capture_scaffold_strict_v2.py --set A --n 800 --seed 42
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
    smi = row["smiles"]  # already resolved to the current v2 standardized smiles by caller
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
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    leak_idx, eval_sets, full_df, scaffold_key_by_smiles, std_cache = load_index()
    key = {"A": "set_A_approved", "B": "set_B_clinical"}[args.set]
    all_rows = list(eval_sets[key])

    # resolve each drug's CURRENT (v2) standardized smiles; drop drugs absent
    # from the v2 index (disclosed coverage change, see module docstring).
    # Checking std_cache alone is NOT enough -- std_cache covers every
    # distinct compound across ALL Phase 2 stages, but a compound can have
    # zero surviving native-human pairs after Stage 7's filtering (e.g. all
    # its pairs referenced a non-single-protein target) and so be absent
    # from v2_index/compounds.csv.gz specifically. Confirmed as a real bug,
    # not hypothetical: Set B's first full run crashed with a KeyError on
    # exactly this case before this check was added.
    v2_index_smiles = set(full_df["smiles"].unique())
    rows = []
    n_dropped_no_v2 = 0
    for r in all_rows:
        v2_smi = std_cache.get(r["molecule_chembl_id"])
        if not v2_smi or v2_smi not in v2_index_smiles:
            n_dropped_no_v2 += 1
            continue
        r2 = dict(r)
        r2["smiles"] = v2_smi
        rows.append(r2)
    print(f"Set {args.set}: {len(all_rows)} total drugs, {len(rows)} present in v2 index "
          f"({100*len(rows)/len(all_rows):.1f}%), {n_dropped_no_v2} dropped", flush=True)

    rng = random.Random(args.seed)
    rng.shuffle(rows)
    if args.n > 0:
        rows = rows[:args.n]

    out_path = os.path.join(RESULTS_DIR, f"capture_scaffold_strict_v2_{args.set}.jsonl")
    print(f"Capturing Set {args.set} (scaffold-strict, v2 index): {len(rows)} queries "
          f"(seed={args.seed}) -> {out_path}", flush=True)
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
