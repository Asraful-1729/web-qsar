"""
Phase 0 (target_prediction_v2_procedure_rev3.md, Sections 3 and 4) -- turns
the raw ChEMBL pull (fetch_chembl_reference.py's output) and the FROZEN v1
index (target_fishing_v1_freeze/, untouched -- see that directory's own
V1_FREEZE.md) into the four evaluation sets and the three-level ground truth
Section 3.1 asks for.

Why matching is done by STANDARDIZED SMILES, not molecule_chembl_id:
v1's own compounds.csv.gz never stored molecule_chembl_id (see
build_target_fishing_index.py's docstring: it aggregates ACROSS
molecule_chembl_ids sharing a standardized structure before writing output,
and the original per-row activities.jsonl this was built from is no longer
on disk -- confirmed absent from the repo). So a drug is "in the v1 index"
iff its standardized SMILES matches a standardized SMILES already indexed.
This is the SAME standardization used at index-build time (salt-strip +
canonicalize, reused verbatim from build_target_fishing_index.py) so the
match is exact, not a fuzzy fallback.

Consequence worth flagging up front (Phase 0 measures this, doesn't assume
it): a drug can fail to appear in v1's OWN evaluable set for two entirely
different reasons that must not be conflated --
  (a) it is a genuine coverage gap (ChEMBL has no qualifying bioactivity
      row for it at all under v1's activity filter), or
  (b) it exists in ChEMBL bioactivity data but its specific salt/tautomer
      form standardizes differently than however v1 indexed it.
This script reports the match rate; build_eval_sets.py does not try to
paper over unmatched drugs with fuzzy matching, because that would hide
exactly the kind of "how much of the drug-realistic evaluation is even
attemptable" number Phase 0 needs to see honestly.

Three-level ground truth (Section 3.1), per matched drug:
  - annotated_targets:    every target_chembl already linked to this
                           compound's SMILES in the v1 index (i.e. what v1
                           itself would call "known" for this compound,
                           under its own activity filter).
  - mechanistic_targets:  the subset of annotated_targets that ALSO appear
                           in ChEMBL's own curated mechanism.json table for
                           that molecule_chembl_id (any row, not just
                           molecular_mechanism=1) -- operationalises
                           "stronger evidence" as "ChEMBL curators reviewed
                           and recorded a specific mechanism," which is a
                           real, defensible, but coarser proxy than
                           "at least two independent assays" (that count is
                           not recoverable from v1's index -- see module
                           docstring above -- so this script does NOT claim
                           to implement the assay-count criterion; only the
                           curated-mechanism criterion).
  - primary_targets:      the subset of mechanistic_targets where ChEMBL
                           flags molecular_mechanism=1 (its own "this is the
                           molecular mechanism of action," as opposed to an
                           ancillary/off-target pharmacology entry).
Multi-target drugs keep ALL primary targets recorded by ChEMBL (rare but
real, e.g. dual-mechanism drugs) -- Section 3.1's "count a hit if any
annotated mechanism target is retrieved, use the best rank" is applied at
SCORING time (metrics.py), not by collapsing here.

Sets:
  A: max_phase == 4 (approved) drugs matched into the v1 index.
  B: max_phase in {1, 2, 3} (clinical) drugs matched into the v1 index.
     NOTE -- this is NOT "clinical AND preclinical" as Section 4.1 names
     Set B. ChEMBL's max_phase field has no defined value for genuine
     preclinical-only compounds (max_phase is null/0 for the vast majority
     of the database, which is NOT the same population as "preclinical
     drug candidates" -- it is mostly ordinary screening compounds with no
     development-phase information at all). A real preclinical set would
     need a curated source (e.g. a specific preclinical pipeline database)
     that is out of scope for this pull. Set B here is clinical-only
     (phases 1-3); the gap is recorded in the Phase 0 report rather than
     silently relabeling ordinary ChEMBL compounds as "preclinical."
  C: scaffold-held-out -- NOT built here. Reuses the existing, already-
     validated scaffold-holdout machinery in
     target_fishing_v1_freeze/code/target_fishing_benchmark.py verbatim
     (see run_phase0.py), because that script already IS a correct
     novel-chemistry generalization test; duplicating it here would risk a
     silent behavioural drift between "the v1 benchmark" and "Phase 0's
     idea of the v1 benchmark."
  D: sparse-target -- built here as a query set: matched drugs (from A + B)
     whose annotated_targets include at least one target with ORIGINAL
     reference depth (n distinct compounds in the FULL v1 index, before any
     leakage removal) of 10 or fewer. This directly operationalises H8/H5's
     "does the method still work when there isn't much ChEMBL evidence,"
     using drug-realistic queries rather than sampling sparse-target
     compounds from the general pool (which the existing scaffold-holdout
     benchmark's own stratification already covers separately, see
     run_phase0.py).

Output: data/eval_sets.json --
  {
    "set_A_approved": [ {smiles, molecule_chembl_id, max_phase,
                          annotated_targets, mechanistic_targets,
                          primary_targets}, ... ],
    "set_B_clinical": [ ... same shape ... ],
    "set_D_sparse_target": [ ... same shape, subset of A+B ... ],
    "match_stats": {...}
  }

Usage: python3 build_eval_sets.py
"""
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "..", "target_fishing_v1_freeze", "code"))

INDEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "target_fishing_v1_freeze", "target_fishing_index")
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(INDEX_DIR)

import pandas as pd  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import SaltRemover  # noqa: E402

RDLogger.DisableLog("rdApp.*")

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_salt_remover = SaltRemover.SaltRemover()

SPARSE_TARGET_DEPTH_MAX = 10  # matches Section 4's "sparse-target compounds (targets with few known actives)"


def _standardize(smiles):
    """Verbatim copy of build_target_fishing_index.py's _standardize() --
       must match exactly, or drug SMILES will fail to line up with how v1
       canonicalized the same structures at index-build time."""
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


def load_v1_index():
    import target_fishing as TF
    _, _, df = TF._load()
    return df


def main():
    t0 = time.time()
    print("Loading v1 frozen index (read-only)...", flush=True)
    df = load_v1_index()
    print(f"  {len(df)} pairs, {df['smiles'].nunique()} compounds, {df['target_chembl'].nunique()} targets", flush=True)

    original_depth_by_target = df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size()

    annotated_by_smiles = {}
    for smi, sub in df.groupby("smiles"):
        annotated_by_smiles[smi] = sorted(sub["target_chembl"].unique().tolist())

    print("Loading ChEMBL drug pull + mechanism table...", flush=True)
    drugs = pd.read_csv(os.path.join(DATA_DIR, "drugs_raw.csv"))
    mechs = pd.read_csv(os.path.join(DATA_DIR, "mechanisms_raw.csv"))
    print(f"  {len(drugs)} drug rows (max_phase 1-4), {len(mechs)} mechanism rows", flush=True)

    mech_targets_by_mol = {}       # molecule_chembl_id -> set(target_chembl) [any mechanism row]
    primary_targets_by_mol = {}    # molecule_chembl_id -> set(target_chembl) [molecular_mechanism == 1]
    for row in mechs.itertuples(index=False):
        mech_targets_by_mol.setdefault(row.molecule_chembl_id, set()).add(row.target_chembl_id)
        if row.molecular_mechanism == 1:
            primary_targets_by_mol.setdefault(row.molecule_chembl_id, set()).add(row.target_chembl_id)

    print("Standardizing + matching drug SMILES against the v1 index...", flush=True)
    n_unparseable = 0
    n_no_match = 0
    n_matched = 0
    matched_rows = []
    for i, row in enumerate(drugs.itertuples(index=False)):
        std = _standardize(row.canonical_smiles)
        if std is None:
            n_unparseable += 1
            continue
        annotated = annotated_by_smiles.get(std)
        if not annotated:
            n_no_match += 1
            continue
        n_matched += 1
        mcid = row.molecule_chembl_id
        mech_t = sorted(mech_targets_by_mol.get(mcid, set()) & set(annotated))
        prim_t = sorted(primary_targets_by_mol.get(mcid, set()) & set(annotated))
        min_depth = min((original_depth_by_target.get(t, 0) for t in annotated), default=None)
        matched_rows.append({
            "molecule_chembl_id": mcid,
            "smiles": std,
            "max_phase": int(row.max_phase),
            "annotated_targets": annotated,
            "mechanistic_targets": mech_t,
            "primary_targets": prim_t,
            "min_annotated_target_reference_depth": int(min_depth) if min_depth is not None else None,
        })
        if (i + 1) % 3000 == 0:
            print(f"  {i+1}/{len(drugs)} processed ({time.time()-t0:.0f}s)", flush=True)

    print(f"  unparseable: {n_unparseable}, no match in v1 index: {n_no_match}, matched: {n_matched}", flush=True)

    by_mcid = {r["molecule_chembl_id"]: r for r in matched_rows}
    set_a = [r for r in matched_rows if r["max_phase"] == 4]
    set_b = [r for r in matched_rows if r["max_phase"] in (1, 2, 3)]
    set_d = [r for r in matched_rows
             if r["min_annotated_target_reference_depth"] is not None
             and r["min_annotated_target_reference_depth"] <= SPARSE_TARGET_DEPTH_MAX]

    n_with_mechanistic = sum(1 for r in matched_rows if r["mechanistic_targets"])
    n_with_primary = sum(1 for r in matched_rows if r["primary_targets"])

    out = {
        "match_stats": {
            "n_drug_rows_pulled": len(drugs),
            "n_unparseable": n_unparseable,
            "n_no_match_in_v1_index": n_no_match,
            "n_matched_total": n_matched,
            "n_matched_with_any_mechanistic_target": n_with_mechanistic,
            "n_matched_with_any_primary_target": n_with_primary,
            "pct_matched_with_primary_target": round(100 * n_with_primary / n_matched, 1) if n_matched else None,
            "sparse_target_depth_threshold": SPARSE_TARGET_DEPTH_MAX,
            "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
        "set_A_approved": set_a,
        "set_B_clinical": set_b,
        "set_D_sparse_target": set_d,
    }
    out_path = os.path.join(DATA_DIR, "eval_sets.json")
    with open(out_path, "w") as f:
        json.dump(out, f)
    print(f"\nSet A (approved, matched): {len(set_a)}", flush=True)
    print(f"Set B (clinical phase 1-3, matched): {len(set_b)}", flush=True)
    print(f"Set D (sparse-target subset of A+B): {len(set_d)}", flush=True)
    print(f"With >=1 primary/mechanism-of-action target: {n_with_primary} ({out['match_stats']['pct_matched_with_primary_target']}%)", flush=True)
    print(f"\nWrote {out_path} ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
