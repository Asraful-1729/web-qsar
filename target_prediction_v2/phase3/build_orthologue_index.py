"""
Phase 3A orthologue ablation prep: v1-format index (compounds.csv.gz +
fingerprints.npz) built from stage7_orthologue_pairs.jsonl instead of the
native-human pairs -- same conversion pattern as
phase2/build_v2_index_files.py, reused for the much smaller orthologue
tier (34,080 pairs) so the already-validated target_fishing.py search
logic can run against it too.

Usage: python3 build_orthologue_index.py
"""
import json
import os
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2_DIR = os.path.join(HERE, "..", "phase2")
DATA_DIR = os.path.join(PHASE2_DIR, "data")
OUT_DIR = os.path.join(HERE, "orthologue_index")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    with open(os.path.join(DATA_DIR, "stage6_fp_meta.json")) as f:
        fp_meta = json.load(f)
    fps_loaded = np.load(os.path.join(DATA_DIR, "stage6_fingerprints.npz"))["fps"]

    t0 = time.time()
    rows, fp_rows = [], []
    n_no_fp = 0
    with open(os.path.join(DATA_DIR, "stage7_orthologue_pairs.jsonl")) as f:
        for line in f:
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
                "species_provenance": r["species_provenance"],
            })
            fp_rows.append(meta["fp_index"])

    print(f"{len(rows)} orthologue pairs with a valid fingerprint, {n_no_fp} dropped", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "compounds.csv.gz"), index=False, compression="gzip")
    fps_out = fps_loaded[np.array(fp_rows)]
    np.savez_compressed(os.path.join(OUT_DIR, "fingerprints.npz"), fps=fps_out)
    print(f"Wrote index, shape {fps_out.shape} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
