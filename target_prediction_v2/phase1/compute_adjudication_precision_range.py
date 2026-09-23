"""
Closes the loop on the adjudication study (PHASE1_ADJUDICATION_STUDY_STATUS.md,
BUILD_PLAN.md gate G1): once a qualified reviewer fills in the `final_verdict`
column of `data/adjudication_pilot_output.csv`, this script turns that into
the actual precision-range number G1 needs -- "X% of the model's top
false positives are actually real, just unannotated" -- instead of
requiring another engineering pass after the review is done.

Expected `final_verdict` values (free text is fine, but these exact strings
are recognized and counted without further normalization -- pick one per
row, matching `ai_suggested_verdict`'s categories or overriding it):
  - "real_but_unannotated (salt-form gap)"
  - "real_but_unannotated (orthologue/homologue gap)"
  - "real_but_unannotated (other)"
  - "genuine_miss"
  - "uncertain"  -- reviewed but genuinely could not be resolved either way

Rows with a blank `final_verdict` are reported as NOT YET REVIEWED and
excluded from the computed range -- this script never fabricates a verdict
for an unreviewed row.

Usage: python3 compute_adjudication_precision_range.py
"""
import csv
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = os.path.join(HERE, "data", "adjudication_pilot_output.csv")

REAL_BUT_UNANNOTATED_PREFIXES = ("real_but_unannotated",)
GENUINE_MISS_VALUES = {"genuine_miss", "genuine_miss (no mechanical evidence found)"}
UNCERTAIN_VALUES = {"uncertain", "uncertain (literature co-occurrence found, not verified)"}


def classify(verdict):
    v = (verdict or "").strip()
    if not v:
        return "not_reviewed"
    if v.startswith(REAL_BUT_UNANNOTATED_PREFIXES):
        return "real_but_unannotated"
    if v in GENUINE_MISS_VALUES:
        return "genuine_miss"
    if v in UNCERTAIN_VALUES:
        return "uncertain"
    return "other_unrecognized"


def main():
    if not os.path.exists(INPUT_PATH):
        print(f"Not found: {INPUT_PATH} -- run adjudication_pilot.py first.")
        return

    with open(INPUT_PATH) as f:
        rows = list(csv.DictReader(f))

    counts = Counter(classify(r.get("final_verdict")) for r in rows)
    n_total = len(rows)
    n_reviewed = n_total - counts["not_reviewed"]

    print(f"{n_total} total candidate pairs; {n_reviewed} have a final_verdict filled in; "
          f"{counts['not_reviewed']} not yet reviewed.\n")

    if counts["other_unrecognized"]:
        print(f"WARNING: {counts['other_unrecognized']} row(s) have a final_verdict value "
              f"that doesn't match any recognized category -- check spelling against the "
              f"values listed in this script's docstring before trusting the range below.\n")

    if n_reviewed == 0:
        floor = sum(1 for r in rows if r.get("ai_salt_form_gap") == "Yes"
                    or r.get("ai_orthologue_gap") == "Yes")
        print(f"No rows reviewed yet -- nothing to report. As an interim, disclosed-as-not-final "
              f"reference point: {floor}/{n_total} ({100*floor/n_total:.1f}%) are already "
              f"database-CONFIRMED real-but-unannotated by the mechanical checks alone "
              f"(PHASE1_ADJUDICATION_STUDY_STATUS.md) -- this is a floor, not the reviewed answer.")
        return

    real = counts["real_but_unannotated"]
    genuine = counts["genuine_miss"]
    uncertain = counts["uncertain"]

    lo = real / n_reviewed
    hi = (real + uncertain) / n_reviewed  # uncertain rows treated as upper-bound-only, not counted as real

    print("Reviewed-row breakdown:")
    print(f"  real_but_unannotated : {real:3d} ({100*real/n_reviewed:.1f}%)")
    print(f"  genuine_miss         : {genuine:3d} ({100*genuine/n_reviewed:.1f}%)")
    print(f"  uncertain            : {uncertain:3d} ({100*uncertain/n_reviewed:.1f}%)")
    print()
    print(f"Precision range for G1 (among REVIEWED rows, n={n_reviewed}):")
    print(f"  {100*lo:.1f}% - {100*hi:.1f}% of the model's top 'false positives' are "
          f"actually real, just unannotated in the current ground truth.")
    print(f"  (low = confirmed real_but_unannotated only; high = also counting still-uncertain "
          f"rows as if they resolve real -- an upper bound, not a point estimate.)")

    if counts["not_reviewed"] > 0:
        print(f"\nNote: {counts['not_reviewed']} row(s) remain unreviewed -- this range covers "
              f"only the {n_reviewed} rows reviewed so far, not the full 153-row worksheet.")


if __name__ == "__main__":
    main()
