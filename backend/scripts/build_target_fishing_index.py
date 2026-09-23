"""
A2 — builds the target-fishing index from a broad ChEMBL bioactivity pull,
independent of this app's own ~64 QSAR-modeled targets (see
scratchpad/fetch_chembl_activities.py for the fetch step and its documented
filter: human targets, IC50/Ki/Kd/EC50, assay_confidence_score>=8, pchembl_
value present).

This REPLACES the earlier version of this script, which only fingerprinted
models/curated/*.csv (this app's own 64 QSAR-panel targets, reused as a
target-fishing pool) — that meant target fishing could never return a target
outside those same 64, no matter what compound was queried. See the ~Sept
2026 audit that found this (177,875 compounds, 64 targets, exactly this
app's own QSAR panel, not a target-fishing-purposed dataset).

IMPORTANT correction from an earlier draft of this docstring:
assay_confidence_score>=8 does NOT reliably restrict activity rows to real
single-protein targets on its own — spot-checked on a sample pull and found
target_pref_name values of "Homo sapiens" (organism-level, no specific
protein) and "HaCaT" (a cell line) passing that filter. ChEMBL's confidence
score describes the ASSAY's target-assignment confidence, which can still
be high for a cell-based or organism-level assay. The activity endpoint's
own `target_type` query param has no effect (confirmed: identical
total_count with or without it — it's silently ignored on this endpoint).
The real fix is a separate cross-reference: `single_protein_target_ids.json`
(built once via `https://www.ebi.ac.uk/chembl/api/data/target.json?
organism=Homo+sapiens&target_type=SINGLE+PROTEIN`, ~5,869 ids as of the
2026 ChEMBL release) — every activity row is checked against this set
before being kept.

Input: activities.jsonl (one raw ChEMBL activity record per line: fields
molecule_chembl_id, canonical_smiles, target_chembl_id, target_pref_name,
pchembl_value, standard_type — see the fetch script for how this was pulled).

Pipeline:
  0. Drop any row whose target_chembl_id isn't in the verified single-
     protein human target id set (see the correction note above).
  1. Standardize each DISTINCT raw SMILES once (salt-strip to the largest
     fragment, canonicalize, drop unparseable structures) — many activity
     rows share the same molecule_chembl_id (tested against several targets,
     or several assays against the same target), so standardizing per
     unique molecule_chembl_id instead of per row avoids redoing this
     millions of times.
  2. Fingerprint (Morgan/ECFP4, radius=2, 2048 bits — same convention as
     similarity.py/generate_decoys.py) and compute the Murcko scaffold once
     per unique standardized SMILES.
  3. Aggregate raw activity rows to one row per (compound, target) pair —
     duplicate assay measurements for the same pair are averaged (mean
     pchembl_value), matching the convention the earlier per-target curated
     CSVs already used (pchembl_value/pchembl_std/n_records).
  4. Map target_chembl_id -> this app's own target_id where one happens to
     exist (docking_registry.json) — most targets in this broader pool
     won't have one; that's expected and fine (see target_fishing.py's
     search(), which now works with target_pref_name alone when there's no
     internal target_id to route into TargetBrowser).

Output (backend/target_fishing_index/):
  - fingerprints.npz : (n_pairs, 256) uint8, packed Morgan bits — ONE row
                        per (compound, target) pair (a compound tested
                        against several targets gets its fingerprint
                        duplicated across several rows; simpler and still
                        cheap — 256 bytes/row even at 2M rows is ~512MB).
  - compounds.csv.gz : target_chembl, target_pref_name, target_id, smiles,
                        pchembl_value, murcko_scaffold
  - manifest.json    : row-aligned with fetch_meta.json's filter/source
                        info, plus this build's own counts/timestamp.

Usage (run from the repo root):
    python -m scripts.build_target_fishing_index --input /path/to/activities.jsonl \\
        --single-protein-ids /path/to/single_protein_target_ids.json [--fetch-meta /path/to/fetch_meta.json]
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, SaltRemover
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

REGISTRY = os.environ.get("DOCKING_REGISTRY", "docking_registry.json")
_salt_remover = SaltRemover.SaltRemover()


def _target_id_map():
    """CHEMBL203 -> CHEMBL203_EGFR, from this app's own docking_registry.json
       — lets a target-fishing hit that HAPPENS to also be one of this app's
       real docking/QSAR targets route straight into TargetBrowser's
       existing selection flow. Most targets in this broader pool won't be
       in the registry at all; those just get target_id=None (see
       target_fishing.py's search(), which handles that)."""
    out = {}
    if os.path.exists(REGISTRY):
        try:
            with open(REGISTRY) as f:
                data = json.load(f)
            targets = data.get("targets", data)
            ids = [t["target_id"] for t in targets] if isinstance(targets, list) else list(targets.keys())
            for tid in ids:
                out[tid.split("_", 1)[0]] = tid
        except Exception:
            pass
    return out


def _standardize(smiles):
    """Parse -> strip salts/counterions (keep the largest organic fragment)
       -> canonicalize. Returns None for anything RDKit can't parse or that
       ends up with zero atoms after stripping (e.g. a bare counterion
       logged as the 'compound' by mistake)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _salt_remover.StripMol(mol, dontRemoveEverything=True)
    except Exception:
        pass
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(mol)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to activities.jsonl")
    ap.add_argument("--single-protein-ids", required=True,
                     help="Path to single_protein_target_ids.json — a JSON list of verified "
                          "human single-protein target_chembl_ids (see module docstring for how "
                          "this was built); every row is checked against this set before being kept.")
    ap.add_argument("--fetch-meta", default=None, help="Path to fetch_meta.json (source/filter provenance)")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "..", "target_fishing_index"))
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    with open(args.single_protein_ids) as f:
        single_protein_ids = set(json.load(f))
    print(f"{len(single_protein_ids)} verified single-protein human target ids loaded", flush=True)

    t0 = time.time()
    print("Pass 1: reading raw activity rows...", flush=True)
    raw_rows = []
    n_lines = 0
    n_non_single_protein = 0
    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n_lines += 1
            try:
                r = json.loads(line)
            except Exception:
                continue
            if not r.get("canonical_smiles") or not r.get("molecule_chembl_id") or not r.get("target_chembl_id"):
                continue
            if r["target_chembl_id"] not in single_protein_ids:
                n_non_single_protein += 1
                continue
            raw_rows.append(r)
    print(f"  {n_lines} lines read, {n_non_single_protein} dropped (not a verified single-protein target), "
          f"{len(raw_rows)} usable rows ({time.time()-t0:.0f}s)", flush=True)

    print("Pass 2: standardizing + fingerprinting distinct compounds...", flush=True)
    std_cache = {}   # molecule_chembl_id -> standardized SMILES or None
    fp_cache = {}     # standardized SMILES -> (packed fingerprint, murcko scaffold)
    n_std_done = 0
    for r in raw_rows:
        mcid = r["molecule_chembl_id"]
        if mcid not in std_cache:
            std_cache[mcid] = _standardize(r["canonical_smiles"])
        std = std_cache[mcid]
        if std is None or std in fp_cache:
            continue
        mol = Chem.MolFromSmiles(std)
        if mol is None:
            fp_cache[std] = None
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
        fp_cache[std] = (packed, scaffold)
        n_std_done += 1
        if n_std_done % 50000 == 0:
            print(f"  {n_std_done} distinct compounds fingerprinted ({time.time()-t0:.0f}s)", flush=True)
    print(f"  {n_std_done} distinct compounds total ({time.time()-t0:.0f}s)", flush=True)

    print("Pass 3: aggregating to (compound, target) pairs...", flush=True)
    # key: (standardized_smiles, target_chembl_id) -> {target_pref_name, pchembl values[]}
    pairs = {}
    for r in raw_rows:
        std = std_cache.get(r["molecule_chembl_id"])
        if not std or fp_cache.get(std) is None:
            continue
        key = (std, r["target_chembl_id"])
        entry = pairs.setdefault(key, {"target_pref_name": r.get("target_pref_name"), "pchembl": []})
        pv = r.get("pchembl_value")
        if pv is not None:
            try:
                entry["pchembl"].append(float(pv))
            except (TypeError, ValueError):
                pass
        if not entry["target_pref_name"] and r.get("target_pref_name"):
            entry["target_pref_name"] = r["target_pref_name"]
    print(f"  {len(pairs)} distinct (compound, target) pairs ({time.time()-t0:.0f}s)", flush=True)

    print("Pass 4: writing index...", flush=True)
    tid_map = _target_id_map()
    fps_out = []
    rows_out = []
    for (std, target_chembl), entry in pairs.items():
        packed, scaffold = fp_cache[std]
        pchembl_vals = entry["pchembl"]
        mean_pchembl = round(sum(pchembl_vals) / len(pchembl_vals), 2) if pchembl_vals else None
        fps_out.append(packed)
        rows_out.append({
            "target_chembl": target_chembl,
            "target_pref_name": entry["target_pref_name"],
            "target_id": tid_map.get(target_chembl),
            "smiles": std,
            "pchembl_value": mean_pchembl,
            "murcko_scaffold": scaffold,
        })

    fps_arr = np.vstack(fps_out)
    np.savez_compressed(os.path.join(out_dir, "fingerprints.npz"), fps=fps_arr)
    pd.DataFrame(rows_out).to_csv(os.path.join(out_dir, "compounds.csv.gz"), index=False, compression="gzip")

    manifest = {
        "n_compound_target_pairs": len(rows_out),
        "n_distinct_compounds": n_std_done,
        "n_distinct_targets": len(set(r["target_chembl"] for r in rows_out)),
        "n_targets_mapped_to_internal_id": sum(1 for r in rows_out if r["target_id"]),
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if args.fetch_meta and os.path.exists(args.fetch_meta):
        with open(args.fetch_meta) as f:
            manifest["fetch"] = json.load(f)
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Done in {time.time()-t0:.0f}s: {json.dumps(manifest, indent=2)}", flush=True)


if __name__ == "__main__":
    main()
