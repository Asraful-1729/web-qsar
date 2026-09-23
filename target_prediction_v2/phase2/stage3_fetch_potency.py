"""
Phase 2 rebuild, Stage 3 (PHASE2_REBUILD_EXECUTION_PLAN.md): ChEMBL
'Potency'-type records -- item 3's decision, closes L3. Confirmed NOT
narrowly PubChem-sourced (predominantly src_id=1 Scientific Literature) --
that mischaracterization in the original plan is corrected in BUILD_PLAN.md
Section 4 already; not repeated as a decision here.

Quality filters applied SERVER-SIDE (cheaper than pull-then-discard,
confirmed both filters work and keep the large majority of the raw pool):
  - assay_type__in=B,F  (drops non-binding/functional codes, e.g. 'T') --
    measured: 2,961,425 of ~2,985,896 raw records survive (98.8%).
  - standard_units=nM   (drops '%'-unit and other non-molar records) --
    measured: 2,985,566 of ~2,985,896 already use nM (99.99%).
pchembl_value is not populated for Potency-type records (confirmed by
inspection during item 3's scoping) -- computed locally here, same nM
formula as Stage 2.

Tags source_type='Potency' for provenance (kept distinct from the four
traditional bioactivity types).

Usage:
    python3 stage3_fetch_potency.py --pilot 20000
    python3 stage3_fetch_potency.py
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
          "assay_chembl_id,document_chembl_id,canonical_smiles,src_id")

QUERY = (f"{BASE}/activity.json?target_organism=Homo+sapiens"
         f"&standard_type=Potency&confidence_score__gte=8"
         f"&standard_value__isnull=false&standard_units=nM&assay_type__in=B,F"
         f"&only={FIELDS}")


def transform(r):
    try:
        value_nm = float(r["standard_value"])
        if value_nm <= 0:
            return None
    except (TypeError, ValueError):
        return None
    r["pchembl_value_computed"] = round(9 - math.log10(value_nm), 3)
    r["source_type"] = "Potency"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    suffix = "_pilot" if args.pilot else ""
    out_path = os.path.join(DATA_DIR, f"stage3_potency{suffix}.jsonl")
    ckpt_path = os.path.join(DATA_DIR, f"stage3_completed_offsets{suffix}.json")

    run_concurrent_fetch(QUERY, out_path, ckpt_path, pilot=args.pilot, row_transform=transform)


if __name__ == "__main__":
    main()
