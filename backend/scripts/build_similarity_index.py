"""
A3 — builds the similarity-search index from a raw COCONUT database export
(https://coconut.naturalproducts.net/download, CC0-licensed — free to use,
modify, and redistribute, no attribution required).

Input: the COCONUT "Lite" CSV (~550MB extracted, ~738K compounds as of the
September 2026 release) — too large to process at query time or bundle
directly, so this is a one-off, offline build step (matching this project's
existing pattern: panel_results_v2.csv, docking_registry.json, and the
QSAR/docking model buckets are all pre-built offline, not computed live).

Output (both written to --out-dir, default backend/similarity_index/):
  - fingerprints.npy   : (n_compounds, 2048) uint8, PACKED Morgan/ECFP4
                         bits (radius=2, 2048 bits — same convention as
                         serving/featurize.py and scripts/generate_decoys.py)
  - compounds.csv      : id, smiles, name, molecular_weight, murcko_framework,
                         np_classifier_pathway/superclass/class,
                         chemical_super_class — everything the search UI
                         needs, without re-parsing the full ~550MB raw file.
  - manifest.json       : row count, build timestamp, source info.

Both output files are row-aligned (row i in compounds.csv <-> row i in
fingerprints.npy) — similarity.py's search loads them together and relies
on this.

Usage:
    python -m scripts.build_similarity_index --input /path/to/coconut_csv_lite-09-2026.csv
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")  # this dataset has ~738K rows; don't spam stderr per-molecule parse warning

KEEP_COLUMNS = [
    "identifier", "canonical_smiles", "name", "molecular_weight", "murcko_framework",
    "np_classifier_pathway", "np_classifier_superclass", "np_classifier_class", "chemical_super_class",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to the raw COCONUT CSV (e.g. coconut_csv_lite-09-2026.csv)")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "..", "similarity_index"))
    ap.add_argument("--chunksize", type=int, default=20000)
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    fps_out = []
    rows_out = []
    n_total = 0
    n_parsed = 0
    t0 = time.time()

    for chunk in pd.read_csv(args.input, usecols=KEEP_COLUMNS, chunksize=args.chunksize, low_memory=False):
        for row in chunk.itertuples(index=False):
            n_total += 1
            smi = row.canonical_smiles
            if not isinstance(smi, str) or not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = morgan.GetFingerprint(mol)
            arr = np.zeros((2048,), dtype=np.uint8)
            for bit in fp.GetOnBits():
                arr[bit] = 1
            fps_out.append(np.packbits(arr))
            rows_out.append({
                "id": row.identifier, "smiles": Chem.MolToSmiles(mol), "name": row.name,
                "molecular_weight": row.molecular_weight, "murcko_framework": row.murcko_framework,
                "np_classifier_pathway": row.np_classifier_pathway,
                "np_classifier_superclass": row.np_classifier_superclass,
                "np_classifier_class": row.np_classifier_class,
                "chemical_super_class": row.chemical_super_class,
            })
            n_parsed += 1
        print(f"  {n_total} read, {n_parsed} parsed OK ({time.time()-t0:.0f}s elapsed)", flush=True)

    fps_arr = np.vstack(fps_out)
    np.save(os.path.join(out_dir, "fingerprints.npy"), fps_arr)
    pd.DataFrame(rows_out).to_csv(os.path.join(out_dir, "compounds.csv"), index=False)
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump({
            "n_compounds": n_parsed, "n_input_rows": n_total,
            "source": "COCONUT (coconut.naturalproducts.net), CC0 license",
            "fingerprint": "Morgan/ECFP4 radius=2 fpSize=2048, packed bits",
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }, f, indent=2)

    print(f"Done: {n_parsed}/{n_total} compounds indexed in {time.time()-t0:.0f}s -> {out_dir}")


if __name__ == "__main__":
    main()
