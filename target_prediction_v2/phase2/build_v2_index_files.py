"""
Phase 3A prep: convert Phase 2's rebuilt data (stage7_native_human_pairs.jsonl
+ stage6_fingerprints.npz/stage6_fp_meta.json) into v1's exact index file
format (compounds.csv.gz + fingerprints.npz, one row per (compound,target)
pair, fingerprint duplicated per row) -- so the already-validated
target_fishing.py / capture_scaffold_strict.py / leakage.py pipeline can be
reused UNCHANGED against the new index, just by pointing
TARGET_FISHING_INDEX_DIR at this new directory. Rebuilding that pipeline
from scratch would risk silently drifting from the exact, already-proven
search/aggregation logic (_search_against's vectorized groupby, the
evidence_score formula, etc.) -- reuse is strictly safer here.

target_pref_name / target_id are written as empty (None) -- confirmed by
reading target_fishing.py's _search_against(): both are carried through
via groupby(...).agg(..., "first") purely for DISPLAY, never used in the
actual Tanimoto/aggregation math. Not fetched in any Phase 2 stage (an
acceptable, disclosed gap for this refit -- the density-adaptive rule
doesn't need target display names).

Output: v2_index/compounds.csv.gz, v2_index/fingerprints.npz
(row-aligned, same shape/dtype convention as v1's).

Usage: python3 build_v2_index_files.py
"""
import json
import os
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
OUT_DIR = os.path.join(HERE, "v2_index")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading Stage 6 fingerprint metadata + fingerprints...", flush=True)
    with open(os.path.join(DATA_DIR, "stage6_fp_meta.json")) as f:
        fp_meta = json.load(f)
    fps_loaded = np.load(os.path.join(DATA_DIR, "stage6_fingerprints.npz"))["fps"]
    print(f"{len(fp_meta)} distinct structures, {fps_loaded.shape} fingerprint array", flush=True)

    t0 = time.time()
    rows = []
    fp_rows = []
    n_no_fp = 0
    with open(os.path.join(DATA_DIR, "stage7_native_human_pairs.jsonl")) as f:
        for i, line in enumerate(f):
            r = json.loads(line)
            smi = r["smiles"]
            meta = fp_meta.get(smi)
            if meta is None or meta.get("fp_index") is None:
                n_no_fp += 1
                continue
            rows.append({
                "target_chembl": r["target_chembl_id"],
                "target_pref_name": None,
                "target_id": None,
                "smiles": smi,
                "pchembl_value": r["pchembl_value"],
                "murcko_scaffold": r["murcko_scaffold"],
            })
            fp_rows.append(meta["fp_index"])
            if (i + 1) % 500000 == 0:
                print(f"  {i+1} pairs processed ({time.time()-t0:.0f}s)", flush=True)

    print(f"{len(rows)} pairs with a valid fingerprint, {n_no_fp} dropped (no fingerprint)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "compounds.csv.gz"), index=False, compression="gzip")
    print(f"Wrote compounds.csv.gz ({time.time()-t0:.0f}s)", flush=True)

    fps_out = fps_loaded[np.array(fp_rows)]
    np.savez_compressed(os.path.join(OUT_DIR, "fingerprints.npz"), fps=fps_out)
    print(f"Wrote fingerprints.npz, shape {fps_out.shape} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
