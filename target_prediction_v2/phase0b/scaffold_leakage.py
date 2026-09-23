"""
Phase 0b, P0b-10 (scaffold-level half only -- see module docstring's caveat
on document-level leakage, which stays untested by design).

Phase 0's leakage control (and capture.py's, reused here) removes only
near-duplicate COMPOUNDS (Tanimoto >= cutoff under two fingerprints). It
does NOT remove the query's whole analogue series -- other compounds
sharing the same Bemis-Murcko scaffold but similar only in the 0.5-0.9
Tanimoto range, well below the near-duplicate cutoff, stay in the
reference. Rev 5 (Section 2's leakage row, Section 5's P0b-10) flags this
as the real, still-untested leakage risk: a query's whole chemical series
being available to the search, not just its near-exact duplicates.

This script re-runs a SAMPLE of queries with a strictly harsher reference
reduction -- remove EVERY compound sharing the query's own Murcko scaffold,
in addition to the standard near-duplicate removal -- and compares recovery
against the same queries' near-duplicate-only numbers (already captured).
A large drop under scaffold removal means analogue-series leakage was
doing real work in the near-duplicate-only setting; a small drop means it
was not.

Document-level leakage (same ChEMBL publication/patent reporting many
near-identical compounds against one target) is NOT tested here -- it
needs document_chembl_id for the full ~1.3M-pair reference, which v1's
aggregated index does not carry and which is out of scope to pull for the
whole index in Phase 0b's "cheap counting" spirit (rev5 explicitly
anticipates this staying open: "Scaffold- and document-level holdout
remain untested and remain the real risk"). Left as an explicit Phase 1
benchmark-card item (four split types already includes document/assay-
campaign).

Usage: python3 scaffold_leakage.py --set A --n 100
"""
import argparse
import json
import os
import random
import sys

import numpy as np

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index"))

import target_fishing as TF  # noqa: E402
import leakage as LK  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["A", "B", "D"])
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE0_DIR, "data", "leakage_fp2.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])

    with open(os.path.join(PHASE0_DIR, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)
    key = {"A": "set_A_approved", "B": "set_B_clinical", "D": "set_D_sparse_target"}[args.set]
    rows = list(eval_sets[key])
    rng = random.Random(args.seed)
    rng.shuffle(rows)
    rows = rows[:args.n]

    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))

    results = []
    for i, row in enumerate(rows):
        smi = row["smiles"]
        annotated = set(row["annotated_targets"])

        # standard (near-duplicate-only) leakage removal, as in capture.py / Phase 0
        reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, 0.95)
        q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)
        results_nd = TF._search_against(reduced_fps, reduced_pop, reduced_df, q_packed, q_pop,
                                         threshold=0.0, max_compounds_per_target=1)
        rank_nd = None
        for j, r in enumerate(results_nd):
            if r["target_chembl"] in annotated:
                rank_nd = j + 1
                break

        # additionally strip the query's WHOLE scaffold group
        q_scaffold = scaffold_key(smi, scaffold_by_smiles.get(smi))
        same_scaffold_smiles = set(full_df.loc[full_df["smiles"].map(
            lambda s: scaffold_key(s, scaffold_by_smiles.get(s))) == q_scaffold, "smiles"])
        extra_mask = ~reduced_df["smiles"].isin(same_scaffold_smiles).values
        sc_fps = reduced_fps[extra_mask]
        sc_pop = reduced_pop[extra_mask]
        sc_df = reduced_df.loc[extra_mask]
        results_sc = TF._search_against(sc_fps, sc_pop, sc_df, q_packed, q_pop,
                                         threshold=0.0, max_compounds_per_target=1)
        rank_sc = None
        for j, r in enumerate(results_sc):
            if r["target_chembl"] in annotated:
                rank_sc = j + 1
                break

        results.append({
            "smiles": smi, "scaffold_group_size": len(same_scaffold_smiles),
            "rank_near_dup_only": rank_nd, "rank_plus_scaffold_removed": rank_sc,
        })
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(rows)}", flush=True)

    n = len(results)

    def topk(key_, k):
        return sum(1 for r in results if r[key_] is not None and r[key_] <= k) / n

    out = {
        "set": args.set, "n": n, "seed": args.seed,
        "near_dup_only": {"top_1": round(topk("rank_near_dup_only", 1), 4), "top_10": round(topk("rank_near_dup_only", 10), 4)},
        "plus_scaffold_removed": {"top_1": round(topk("rank_plus_scaffold_removed", 1), 4), "top_10": round(topk("rank_plus_scaffold_removed", 10), 4)},
        "mean_scaffold_group_size": round(sum(r["scaffold_group_size"] for r in results) / n, 2),
        "per_query": results,
    }
    out_path = os.path.join(HERE, "results", f"scaffold_leakage_{args.set}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "per_query"}, indent=2), flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
