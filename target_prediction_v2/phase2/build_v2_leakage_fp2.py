"""
Phase 3A prep: second, independent fingerprint (RDKit topological,
2048-bit) for every distinct compound in Phase 2's rebuilt index --
mirrors phase0/build_leakage_fp2.py exactly, against the new, larger
compound set (1,590,210 vs. v1's 857,232), for the same reason: a near-
duplicate that folds/collides similarly under Morgan/ECFP4 (the main
search fingerprint) is unlikely to do so under a genuinely different
algorithm family too, per Phase 0's leakage-sensitivity finding (L5).

Requires RDKit -- run with the vsdock conda env's python3.

Output: data/v2_leakage_fp2.npz -- {"smiles": ..., "fps": (n,256) uint8}

Usage (from the vsdock env): python3 build_v2_leakage_fp2.py
"""
import json
import os
import time

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")


def packed_topo_fp(mol):
    fp = Chem.RDKFingerprint(mol, fpSize=2048)
    bits = np.zeros(2048, dtype=np.uint8)
    on = list(fp.GetOnBits())
    if on:
        bits[on] = 1
    return np.packbits(bits)


def main():
    with open(os.path.join(DATA_DIR, "stage6_fp_meta.json")) as f:
        fp_meta = json.load(f)
    distinct = list(fp_meta.keys())
    print(f"{len(distinct)} distinct compounds to fingerprint (topological, 2048-bit)", flush=True)

    out_path = os.path.join(DATA_DIR, "v2_leakage_fp2.npz")
    ckpt_path = os.path.join(DATA_DIR, "v2_leakage_fp2_progress.json")

    start_i = 0
    fps, kept_smiles = [], []
    if os.path.exists(ckpt_path) and os.path.exists(out_path):
        with open(ckpt_path) as f:
            start_i = json.load(f)["n_done"]
        loaded = np.load(out_path, allow_pickle=True)
        fps = list(loaded["fps"])
        kept_smiles = list(loaded["smiles"])
        print(f"Resuming from {start_i}", flush=True)

    t0 = time.time()
    n_fail = 0
    for i in range(start_i, len(distinct)):
        smi = distinct[i]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_fail += 1
            continue
        fps.append(packed_topo_fp(mol))
        kept_smiles.append(smi)
        if (i + 1) % 50000 == 0:
            fps_arr = np.vstack(fps)
            np.savez_compressed(out_path, smiles=np.array(kept_smiles, dtype=object), fps=fps_arr)
            with open(ckpt_path, "w") as f:
                json.dump({"n_done": i + 1}, f)
            print(f"  {i+1}/{len(distinct)} ({time.time()-t0:.0f}s)", flush=True)

    fps_arr = np.vstack(fps)
    np.savez_compressed(out_path, smiles=np.array(kept_smiles, dtype=object), fps=fps_arr)
    with open(ckpt_path, "w") as f:
        json.dump({"n_done": len(distinct)}, f)
    print(f"Done: {len(kept_smiles)} fingerprinted, {n_fail} failed -> {out_path} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
