"""
Phase 0b -- supplementary, cheap capture: for each sampled query drug, its
OWN (target_chembl -> pchembl_value) pairs as recorded directly in v1's
index (i.e., the potency of the query compound itself against each of its
annotated targets). Needed for P0b-5's ground-truth-filtering variant
("only count a target as truly annotated if the QUERY DRUG's own activity
against it clears the threshold") -- distinct from index-filtering (removing
weak-potency NEIGHBOUR evidence from voting, already covered by capture.py's
per-neighbour pchembl + score.py's potency_fn).

No leakage removal, no Tanimoto -- just a direct lookup against the full
v1 index (loaded once, in-memory groupby), so this is fast even for
hundreds of queries.

Output: data/query_own_potency.json -- {smiles: {target_chembl: pchembl_or_null}}
"""
import json
import os
import random
import sys

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index"))

import target_fishing as TF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    _, _, full_df = TF._load()

    with open(os.path.join(PHASE0_DIR, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)

    all_rows = []
    for key, n, seed in [("set_A_approved", 200, 42), ("set_B_clinical", 200, 42), ("set_D_sparse_target", -1, 42)]:
        rows = list(eval_sets[key])
        rng = random.Random(seed)
        rng.shuffle(rows)
        if n > 0:
            rows = rows[:n]
        all_rows.extend(rows)

    smiles_needed = set(r["smiles"] for r in all_rows)
    sub = full_df[full_df["smiles"].isin(smiles_needed)][["smiles", "target_chembl", "pchembl_value"]]

    out = {}
    for smi, grp in sub.groupby("smiles"):
        out[smi] = {
            row.target_chembl: (None if row.pchembl_value != row.pchembl_value else float(row.pchembl_value))
            for row in grp.itertuples(index=False)
        }

    out_path = os.path.join(HERE, "data", "query_own_potency.json")
    with open(out_path, "w") as f:
        json.dump(out, f)
    print(f"{len(out)} queries -> {out_path}")


if __name__ == "__main__":
    main()
