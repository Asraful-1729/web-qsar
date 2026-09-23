"""
Phase 2 rebuild, Stage 2 (PHASE2_REBUILD_EXECUTION_PLAN.md): censored
('<'/'<=') relation records -- item 3's decision. ~174K records (measured
in item 3's scoping pass; restricting to standard_units=nM loses only
~0.26% of them -- 173,779 of 174,234 -- so filtered server-side here
rather than handling arbitrary units).

pchembl_value is NEVER populated by ChEMBL for non-'=' relations (the L4
finding) -- computed locally here from standard_value (nM): pChEMBL =
9 - log10(value_nM), the standard nM->pChEMBL conversion (pChEMBL =
-log10(molar); nM = 1e-9 M, so -log10(value_nM * 1e-9) = 9 - log10(value_nM)).

Tags every row is_censored=True (a LOWER BOUND on potency, not exact --
item 1/3's decision: Phase 3A's weight-fitting must not treat this as an
exact measurement).

Usage:
    python3 stage2_fetch_censored.py --pilot 20000
    python3 stage2_fetch_censored.py
"""
import argparse
import math
import os

from fetch_common import run_concurrent_fetch

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")

FIELDS = ("molecule_chembl_id,target_chembl_id,standard_type,standard_relation,"
          "standard_value,standard_units,confidence_score,assay_type,"
          "assay_chembl_id,document_chembl_id,canonical_smiles")

QUERY = (f"{BASE}/activity.json?target_organism=Homo+sapiens"
         f"&standard_type__in=IC50,Ki,Kd,EC50&standard_relation__in=%3C,%3C%3D"
         f"&confidence_score__gte=8&standard_value__isnull=false&standard_units=nM"
         f"&only={FIELDS}")


def transform(r):
    try:
        value_nm = float(r["standard_value"])
        if value_nm <= 0:
            return None
    except (TypeError, ValueError):
        return None
    r["pchembl_value_computed"] = round(9 - math.log10(value_nm), 3)
    r["is_censored"] = True
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    suffix = "_pilot" if args.pilot else ""
    out_path = os.path.join(DATA_DIR, f"stage2_censored{suffix}.jsonl")
    ckpt_path = os.path.join(DATA_DIR, f"stage2_completed_offsets{suffix}.json")

    run_concurrent_fetch(QUERY, out_path, ckpt_path, pilot=args.pilot, row_transform=transform)


if __name__ == "__main__":
    main()
