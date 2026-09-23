"""
Phase 2 rebuild, Stage 6 (PHASE2_REBUILD_EXECUTION_PLAN.md): standardize
every distinct compound across all Stage 1-5b outputs and fingerprint
once per unique standardized structure.

Reuses v1's exact conventions (build_target_fishing_index.py): Morgan
radius=2, fpSize=2048, packed via np.packbits; Murcko scaffold via
MurckoScaffold.GetScaffoldForMol -- same fingerprint/scaffold code, only
the standardization step changes.

THE FIX (item 5's decision): SaltRemover.StripMol() -> rdMolStandardize.
FragmentParent(), confirmed live to correctly collapse n:1 stoichiometry
salts (the ladostigil tartrate case) that StripMol leaves as disconnected
duplicate fragments -- 2,167 rows (~0.17%) of v1's own shipped index were
found to carry this exact bug.

Checkpointed (a real lesson from this session -- Stage 5a lost 434s of
progress to a single crash near the end with no resume support; this
stage processes ~1.6M distinct compounds, too large to risk the same
mistake): periodic saves of the standardization + fingerprint cache,
resumable.

Requires RDKit -- run with the vsdock conda env's python3, not the
default interpreter (confirmed via `conda env list` that vsdock has
RDKit 2026.03.6; no other available env does).

Usage (from the vsdock env):
    python3 stage6_standardize.py
"""
import json
import os
import time

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
CHECKPOINT_EVERY = 25000

SOURCE_FILES = [
    "stage1_activities.jsonl", "stage2_censored.jsonl", "stage3_potency.jsonl",
    "stage4_inactives.jsonl", "stage5b_orthologue_activities.jsonl",
]


def build_distinct_compounds():
    out_path = os.path.join(DATA_DIR, "stage6_distinct_compounds.json")
    if os.path.exists(out_path):
        with open(out_path) as f:
            return json.load(f)
    distinct = {}
    for fn in SOURCE_FILES:
        path = os.path.join(DATA_DIR, fn)
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                mcid, smi = r.get("molecule_chembl_id"), r.get("canonical_smiles")
                if mcid and smi and mcid not in distinct:
                    distinct[mcid] = smi
        print(f"  after {fn}: {len(distinct)} cumulative distinct compounds", flush=True)
    with open(out_path, "w") as f:
        json.dump(distinct, f)
    return distinct


def standardize(smiles):
    """FragmentParent (item 5's fix, replaces SaltRemover.StripMol) -> canonicalize."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.FragmentParent(mol)
    except Exception:
        return None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(mol)


def main():
    print("Building/loading distinct-compound map...", flush=True)
    distinct = build_distinct_compounds()
    print(f"{len(distinct)} distinct compounds total", flush=True)

    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    std_cache_path = os.path.join(DATA_DIR, "stage6_std_cache.json")
    std_cache = {}
    if os.path.exists(std_cache_path):
        with open(std_cache_path) as f:
            std_cache = json.load(f)
        print(f"Resuming standardization: {len(std_cache)} already done", flush=True)

    mol_items = list(distinct.items())
    remaining = [(m, s) for m, s in mol_items if m not in std_cache]
    print(f"{len(remaining)} compounds remaining to standardize", flush=True)

    t0 = time.time()
    for i, (mcid, smi) in enumerate(remaining):
        std_cache[mcid] = standardize(smi)
        if (i + 1) % CHECKPOINT_EVERY == 0:
            with open(std_cache_path, "w") as f:
                json.dump(std_cache, f)
            print(f"  {i+1}/{len(remaining)} standardized this run ({time.time()-t0:.0f}s)", flush=True)
    with open(std_cache_path, "w") as f:
        json.dump(std_cache, f)
    print(f"Standardization done: {len(std_cache)} total ({time.time()-t0:.0f}s)", flush=True)

    print("Grouping by standardized SMILES (this is where salt-form collapsing happens)...", flush=True)
    by_std = {}
    for mcid, std in std_cache.items():
        if std is None:
            continue
        by_std.setdefault(std, []).append(mcid)
    print(f"{len(by_std)} distinct standardized structures "
          f"(collapsed from {len(std_cache)} distinct molecule_chembl_ids, "
          f"{sum(1 for v in by_std.values() if len(v) > 1)} structures have >1 contributing salt-form id)", flush=True)

    fp_cache_meta_path = os.path.join(DATA_DIR, "stage6_fp_meta.json")
    fp_npz_path = os.path.join(DATA_DIR, "stage6_fingerprints.npz")

    fp_meta = {}
    fps = []
    if os.path.exists(fp_cache_meta_path) and os.path.exists(fp_npz_path):
        with open(fp_cache_meta_path) as f:
            fp_meta = json.load(f)
        loaded = np.load(fp_npz_path)
        fps = list(loaded["fps"])
        print(f"Resuming fingerprinting: {len(fp_meta)} already done", flush=True)

    std_list = list(by_std.keys())
    remaining_std = [s for s in std_list if s not in fp_meta]
    print(f"{len(remaining_std)} structures remaining to fingerprint", flush=True)

    t0 = time.time()
    for i, std in enumerate(remaining_std):
        mol = Chem.MolFromSmiles(std)
        if mol is None:
            fp_meta[std] = {"scaffold": None, "fp_index": None,
                             "contributing_molecule_chembl_ids": by_std[std]}
            continue
        fp = morgan.GetFingerprint(mol)
        bits = np.zeros(2048, dtype=np.uint8)
        on = list(fp.GetOnBits())
        if on:
            bits[on] = 1
        packed = np.packbits(bits)
        try:
            scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        except Exception:
            scaffold = None
        fps.append(packed)
        fp_meta[std] = {"scaffold": scaffold, "fp_index": len(fps) - 1,
                         "contributing_molecule_chembl_ids": by_std[std]}

        if (i + 1) % CHECKPOINT_EVERY == 0:
            np.savez_compressed(fp_npz_path, fps=np.array(fps))
            with open(fp_cache_meta_path, "w") as f:
                json.dump(fp_meta, f)
            print(f"  {i+1}/{len(remaining_std)} fingerprinted this run ({time.time()-t0:.0f}s)", flush=True)

    np.savez_compressed(fp_npz_path, fps=np.array(fps))
    with open(fp_cache_meta_path, "w") as f:
        json.dump(fp_meta, f)
    print(f"Fingerprinting done: {len(fp_meta)} total structures, {len(fps)} fingerprints "
          f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
