"""
Phase 0 (Section 4, point 2) -- precomputes a SECOND, independent fingerprint
(RDKit's classic path-based topological fingerprint, 2048 bits) for every
DISTINCT compound in the frozen v1 index. v1's own search fingerprint is
Morgan/ECFP4 (circular); this one is path-based, a genuinely different
algorithm family, so a near-duplicate that happens to fold/collide
similarly under Morgan is unlikely to do so under this one too -- that is
exactly the point of the review's leakage-sensitivity request ("a filter
using the same fingerprint as the model can miss close analogues").

This is a one-time, read-only precompute against the FROZEN snapshot
(target_fishing_v1_freeze/) -- output goes into THIS phase0 workspace, never
back into the freeze.

Output: data/leakage_fp2.npz -- {"smiles": array of standardized SMILES
(distinct compounds, same order as "fps"), "fps": (n, 256) uint8 packed bits}
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "target_fishing_v1_freeze", "target_fishing_index"))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def packed_topo_fp(mol):
    fp = Chem.RDKFingerprint(mol, fpSize=2048)
    bits = np.zeros(2048, dtype=np.uint8)
    on = list(fp.GetOnBits())
    if on:
        bits[on] = 1
    return np.packbits(bits)


def main():
    import target_fishing as TF
    _, _, df = TF._load()
    distinct = df.drop_duplicates(subset="smiles")["smiles"].tolist()
    print(f"{len(distinct)} distinct compounds to fingerprint (topological, 2048-bit)", flush=True)

    t0 = time.time()
    fps = []
    kept_smiles = []
    n_fail = 0
    for i, smi in enumerate(distinct):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_fail += 1
            continue
        fps.append(packed_topo_fp(mol))
        kept_smiles.append(smi)
        if (i + 1) % 50000 == 0:
            print(f"  {i+1}/{len(distinct)} ({time.time()-t0:.0f}s)", flush=True)

    fps_arr = np.vstack(fps)
    out_path = os.path.join(DATA_DIR, "leakage_fp2.npz")
    np.savez_compressed(out_path, smiles=np.array(kept_smiles, dtype=object), fps=fps_arr)
    print(f"Done: {len(kept_smiles)} fingerprinted, {n_fail} failed -> {out_path} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
