"""
Phase 2 rebuild, Stage 4 (PHASE2_REBUILD_EXECUTION_PLAN.md): measured-
inactive evidence -- item 6's decision. Genuinely tested-and-inactive
compounds (standard_relation '>'/'>=', standard_value >=10uM -- symmetric
with the existing active/ground-truth boundary), NOT "never tested"
compounds misused as presumed-negatives (the bias PIDGIN [19] avoids).
Measured in item 6's scoping pass: 490,821 records (487,802 '>' + 3,019
'>=').

pchembl_value is not populated for non-'=' relations (L4) -- computed
locally here for completeness/consistency with Stages 2-3, though the
measured-inactive flag itself (not the exact magnitude) is what Phase 3's
functional-form fitting primarily needs.

Tags is_measured_inactive=True.

Usage:
    python3 stage4_fetch_inactives.py --pilot 20000
    python3 stage4_fetch_inactives.py
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
         f"&standard_type__in=IC50,Ki,Kd,EC50&standard_relation__in=%3E,%3E%3D"
         f"&standard_value__gte=10000&standard_units=nM&confidence_score__gte=8"
         f"&only={FIELDS}")


def transform(r):
    try:
        value_nm = float(r["standard_value"])
        if value_nm <= 0:
            return None
    except (TypeError, ValueError):
        return None
    r["pchembl_value_computed"] = round(9 - math.log10(value_nm), 3)
    r["is_measured_inactive"] = True
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    suffix = "_pilot" if args.pilot else ""
    out_path = os.path.join(DATA_DIR, f"stage4_inactives{suffix}.jsonl")
    ckpt_path = os.path.join(DATA_DIR, f"stage4_completed_offsets{suffix}.json")

    run_concurrent_fetch(QUERY, out_path, ckpt_path, pilot=args.pilot, row_transform=transform)


if __name__ == "__main__":
    main()
