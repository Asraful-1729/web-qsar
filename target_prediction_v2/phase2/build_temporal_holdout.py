"""
Phase 2/Phase 1 follow-up: the temporal holdout, previously fully blocked
(no per-pair date existed anywhere in this program -- PHASE7_RELEASE_READINESS.md).
Unblocked by stage_document_year_backfill.py + Stage 7's first_seen_year
field (earliest resolvable publication year among a pair's contributing
documents).

REAL COVERAGE, disclosed precisely, not glossed over: 1,567,975 of
3,973,653 native-human pairs (39.5%) have a resolvable first_seen_year.
Root cause of the gap: Stage 3 (Potency records, ~half of all raw activity
rows) is 99.997% dominated by a single ChEMBL document (CHEMBL1201862, a
bulk-deposited dataset with no publication year in ChEMBL's own database)
-- confirmed, not a bug in the backfill (Stages 1/2/4 resolve at 90-99.7%
individually). A pair gets a year here if ANY of its contributing rows
(across all 4 native-human stages) resolves one -- so a pair with BOTH
Potency-only evidence AND e.g. Stage 1 evidence still gets a real year.

Partition: same scaffold-group-safe logic as phase1/build_holdouts.py,
but temporally ordered instead of randomly shuffled -- the most RECENT
scaffold groups (by their earliest known evidence year) become the locked
temporal test set, mirroring a genuine "does this generalize to newer
chemistry" holdout, not a random resample.

Scope: only compounds/pairs with a resolvable year participate in the
temporal split (the other 60.5% remain usable for random/scaffold splits,
just not this one) -- disclosed, not silently dropped elsewhere.

Usage: python3 build_temporal_holdout.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
TEST_FRACTION = 0.20


def main():
    print("Loading native-human pairs with a resolvable first_seen_year...", flush=True)
    # compound (smiles) -> earliest year across ALL its pairs (a compound's
    # "first appearance" in this dataset, not per-target)
    compound_year = {}
    compound_scaffold = {}
    n_pairs = 0
    n_pairs_with_year = 0
    with open(os.path.join(DATA_DIR, "stage7_native_human_pairs.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            n_pairs += 1
            smi = r["smiles"]
            compound_scaffold[smi] = r.get("murcko_scaffold") or f"__no_ring__:{smi}"
            yr = r.get("first_seen_year")
            if yr is not None:
                n_pairs_with_year += 1
                if smi not in compound_year or yr < compound_year[smi]:
                    compound_year[smi] = yr

    n_compounds = len(compound_scaffold)
    n_compounds_with_year = len(compound_year)
    print(f"{n_pairs} pairs total, {n_pairs_with_year} with a resolvable year "
          f"({100*n_pairs_with_year/n_pairs:.1f}%)", flush=True)
    print(f"{n_compounds} distinct compounds, {n_compounds_with_year} with a resolvable year "
          f"({100*n_compounds_with_year/n_compounds:.1f}%)", flush=True)

    # group eligible (year-resolved) compounds by scaffold; a scaffold
    # group's year = its EARLIEST compound's year (the group "first
    # appeared" when its first member did)
    scaffold_groups = {}
    for smi, yr in compound_year.items():
        sk = compound_scaffold[smi]
        scaffold_groups.setdefault(sk, {"smiles": [], "year": None})
        scaffold_groups[sk]["smiles"].append(smi)
        if scaffold_groups[sk]["year"] is None or yr < scaffold_groups[sk]["year"]:
            scaffold_groups[sk]["year"] = yr

    print(f"{len(scaffold_groups)} distinct scaffold groups among year-eligible compounds", flush=True)

    # sort scaffold groups OLDEST to NEWEST; the most recent
    # TEST_FRACTION of compounds (by scaffold-group year) become the
    # locked temporal test set
    ordered = sorted(scaffold_groups.items(), key=lambda kv: kv[1]["year"])
    total_eligible = sum(len(v["smiles"]) for _, v in ordered)
    target_test_n = int(round(total_eligible * TEST_FRACTION))

    test_groups, train_groups = [], []
    running = 0
    # walk from the END (most recent) backward to build the test set
    for sk, v in reversed(ordered):
        if running < target_test_n:
            test_groups.append((sk, v))
            running += len(v["smiles"])
        else:
            train_groups.append((sk, v))
    train_groups = list(reversed(train_groups))  # restore chronological order for reporting

    test_smiles = [s for _, v in test_groups for s in v["smiles"]]
    train_smiles = [s for _, v in train_groups for s in v["smiles"]]
    test_years = [v["year"] for _, v in test_groups]
    train_years = [v["year"] for _, v in train_groups]

    out = {
        "method": "scaffold-group-safe, temporally ordered (most recent scaffold groups -> test)",
        "n_native_human_pairs_total": n_pairs,
        "n_pairs_with_resolvable_year": n_pairs_with_year,
        "pct_pairs_with_year": round(100 * n_pairs_with_year / n_pairs, 1),
        "n_compounds_total": n_compounds,
        "n_compounds_with_resolvable_year": n_compounds_with_year,
        "known_coverage_gap_cause": "Stage 3 (Potency) is 99.997% dominated by one undated bulk-deposit "
                                     "document (CHEMBL1201862) -- confirmed root cause, not a backfill bug. "
                                     "Stages 1/2/4 individually resolve at 90-99.7%.",
        "test_fraction_target": TEST_FRACTION,
        "temporal_test": {
            "n_compounds": len(test_smiles),
            "n_scaffold_groups": len(test_groups),
            "year_range": [min(test_years), max(test_years)] if test_years else None,
            "smiles": test_smiles,
        },
        "temporal_train": {
            "n_compounds": len(train_smiles),
            "n_scaffold_groups": len(train_groups),
            "year_range": [min(train_years), max(train_years)] if train_years else None,
        },
        "excluded_no_year": {
            "n_compounds": n_compounds - n_compounds_with_year,
            "note": "these compounds have no resolvable first_seen_year and do not participate in "
                    "the temporal split -- still usable for random/scaffold splits (phase1/data/holdouts.json)",
        },
    }

    out_path = os.path.join(DATA_DIR, "temporal_holdout.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nTemporal test: {len(test_smiles)} compounds, {len(test_groups)} scaffold groups, "
          f"years {min(test_years)}-{max(test_years)}", flush=True)
    print(f"Temporal train: {len(train_smiles)} compounds, {len(train_groups)} scaffold groups, "
          f"years {min(train_years)}-{max(train_years)}", flush=True)
    print(f"Excluded (no resolvable year): {n_compounds - n_compounds_with_year} compounds", flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
