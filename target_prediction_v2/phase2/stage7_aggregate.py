"""
Phase 2 rebuild, Stage 7 (PHASE2_REBUILD_EXECUTION_PLAN.md): aggregate raw
activity rows (Stages 1-5b) to one row per (compound, target) pair,
reusing v1's Pass-3/Pass-4 aggregation pattern with the decisions from
items 1-5:

  - MAX pchembl per pair, not mean (item 3's correction to v1's actual
    code, which used mean -- confirmed by reading build_target_fishing_index.py).
    The max is taken only among genuine positive/censored/potency evidence
    (is_measured_inactive rows are tracked as a flag, never allowed to
    drag the max down -- their own pchembl is definitionally low, <=5, so
    this matters for correctness not just tidiness).
  - is_censored / is_measured_inactive / confidence_score / best_relation
    all taken from the record that produced the winning max-pchembl value
    (or, for is_measured_inactive, True if ANY row for the pair was
    inactive-tagged, regardless of which record won the max).
  - n_documents / assay_types_seen aggregated across ALL rows for the pair
    (every evidence type contributes to these, not just the winner).
  - contributing_molecule_chembl_ids taken directly from Stage 6's
    fp_meta (already computed once per standardized structure).
  - Native-human tier and orthologue tier written to SEPARATE files, never
    merged (item 4's decision) -- orthologue tier aggregated per
    (compound, human_target, species) triple, not collapsed across
    species, since Phase 2's job is making granular data available, not
    deciding how Phase 3A should combine multi-species evidence.

Requires RDKit only insofar as it reads Stage 6's outputs (no RDKit calls
of its own) -- runs fine under the default python3, no special env needed.

Usage: python3 stage7_aggregate.py
"""
import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")

RELATION_PRIORITY = {"=": 0, "<": 1, "<=": 1, ">": 2, ">=": 2}  # for best_relation tie-breaking display only


def load_std_cache():
    with open(os.path.join(DATA_DIR, "stage6_std_cache.json")) as f:
        return json.load(f)


def load_document_year_map():
    """Backfilled separately (stage_document_year_backfill.py) -- same
    class of gap as confidence_score: activity.json doesn't serialize a
    date/year field itself, document.json does (confirmed live). Enables
    a real temporal split -- the one G1 blocker that was previously fully
    unaddressed, per PHASE7_RELEASE_READINESS.md."""
    path = os.path.join(DATA_DIR, "document_year_map.json")
    with open(path) as f:
        return json.load(f)


def load_single_protein_ids():
    """Cross-reference required because confidence_score__gte=8 alone does
    NOT reliably restrict to real single-protein targets -- confirmed both
    by v1's own build script docstring (found organism/cell-line targets
    passing that filter) and independently re-confirmed here: CHEMBL372
    (target_type=ORGANISM, pref_name='Homo sapiens') passed every Stage
    1-4 fetch's confidence_score__gte=8 filter despite its actual assay
    confidence_score being 1 ('Target assigned is non-molecular') --
    18.99% of the pre-fix native-human pairs (931,457 of 4,905,110)
    referenced a non-single-protein target. v1 fixed this by cross-
    checking every row against single_protein_target_ids.json before
    keeping it (Pass 0 of build_target_fishing_index.py); Stage 0 built
    that same file but this aggregation step had not been applying it
    until this fix."""
    with open(os.path.join(DATA_DIR, "single_protein_target_ids.json")) as f:
        return set(json.load(f))


def load_fp_meta():
    with open(os.path.join(DATA_DIR, "stage6_fp_meta.json")) as f:
        return json.load(f)


def load_confidence_map():
    """activity.json never serializes confidence_score itself (confirmed
    empirically -- filterable server-side, but absent from every returned
    record, even unrestricted) -- the real value lives on the Assay,
    reached via assay_chembl_id. stage_confidence_backfill.py resolves
    this separately; without it every row's confidence_score would be null."""
    path = os.path.join(DATA_DIR, "assay_confidence_map.json")
    with open(path) as f:
        return json.load(f)


def row_pchembl(r):
    if r.get("pchembl_value") is not None:
        try:
            return float(r["pchembl_value"])
        except (TypeError, ValueError):
            return None
    if r.get("pchembl_value_computed") is not None:
        return float(r["pchembl_value_computed"])
    return None


def aggregate_native_human(std_cache, confidence_map, single_protein_ids, year_map):
    pairs = {}  # (std_smiles, target_chembl_id) -> accumulator
    files = ["stage1_activities.jsonl", "stage2_censored.jsonl",
             "stage3_potency.jsonl", "stage4_inactives.jsonl"]
    t0 = time.time()
    n_rows = 0
    n_no_target = 0
    n_non_single_protein = 0
    for fn in files:
        with open(os.path.join(DATA_DIR, fn)) as f:
            for line in f:
                r = json.loads(line)
                n_rows += 1
                mcid, tcid = r.get("molecule_chembl_id"), r.get("target_chembl_id")
                if not tcid:
                    n_no_target += 1
                    continue
                if tcid not in single_protein_ids:
                    n_non_single_protein += 1
                    continue
                std = std_cache.get(mcid)
                if not std:
                    continue
                key = (std, tcid)
                acc = pairs.setdefault(key, {
                    "documents": set(), "assay_types": set(),
                    "best_pchembl": None, "best_relation": None,
                    "best_confidence": None, "best_is_censored": False,
                    "is_measured_inactive": False, "first_seen_year": None,
                })
                if r.get("document_chembl_id"):
                    acc["documents"].add(r["document_chembl_id"])
                    yr = year_map.get(r["document_chembl_id"])
                    if yr is not None and (acc["first_seen_year"] is None or yr < acc["first_seen_year"]):
                        acc["first_seen_year"] = yr
                if r.get("assay_type"):
                    acc["assay_types"].add(r["assay_type"])
                if r.get("is_measured_inactive"):
                    acc["is_measured_inactive"] = True
                    continue  # inactive rows never contribute to best_pchembl
                pc = row_pchembl(r)
                if pc is not None and (acc["best_pchembl"] is None or pc > acc["best_pchembl"]):
                    acc["best_pchembl"] = pc
                    acc["best_relation"] = r.get("standard_relation")
                    acc["best_confidence"] = confidence_map.get(r.get("assay_chembl_id"))
                    acc["best_is_censored"] = bool(r.get("is_censored"))
        print(f"  processed {fn} ({n_rows} rows so far, {len(pairs)} pairs so far, "
              f"{time.time()-t0:.0f}s)", flush=True)
    print(f"Native-human: {n_rows} rows, {n_no_target} dropped (no target_chembl_id), "
          f"{n_non_single_protein} dropped (not a verified single-protein target), "
          f"{len(pairs)} distinct (compound,target) pairs", flush=True)
    return pairs


def aggregate_orthologue(std_cache, confidence_map, year_map):
    pairs = {}  # (std_smiles, human_target_chembl_id, species) -> accumulator
    t0 = time.time()
    n_rows = 0
    with open(os.path.join(DATA_DIR, "stage5b_orthologue_activities.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            n_rows += 1
            mcid = r.get("molecule_chembl_id")
            std = std_cache.get(mcid)
            if not std:
                continue
            species = r.get("species_provenance")
            for human_tcid in r.get("human_target_chembl_ids", []):
                key = (std, human_tcid, species)
                acc = pairs.setdefault(key, {
                    "documents": set(), "assay_types": set(),
                    "best_pchembl": None, "best_relation": None,
                    "best_confidence": None, "first_seen_year": None,
                })
                if r.get("document_chembl_id"):
                    acc["documents"].add(r["document_chembl_id"])
                    yr = year_map.get(r["document_chembl_id"])
                    if yr is not None and (acc["first_seen_year"] is None or yr < acc["first_seen_year"]):
                        acc["first_seen_year"] = yr
                if r.get("assay_type"):
                    acc["assay_types"].add(r["assay_type"])
                pc = row_pchembl(r)
                if pc is not None and (acc["best_pchembl"] is None or pc > acc["best_pchembl"]):
                    acc["best_pchembl"] = pc
                    acc["best_relation"] = r.get("standard_relation")
                    acc["best_confidence"] = confidence_map.get(r.get("assay_chembl_id"))
    print(f"Orthologue: {n_rows} rows, {len(pairs)} distinct (compound,human_target,species) pairs "
          f"({time.time()-t0:.0f}s)", flush=True)
    return pairs


def write_output(pairs, fp_meta, out_path, orthologue=False):
    std_to_contrib = {std: v["contributing_molecule_chembl_ids"] for std, v in fp_meta.items()}
    std_to_scaffold = {std: v["scaffold"] for std, v in fp_meta.items()}

    with open(out_path, "w") as f:
        for key, acc in pairs.items():
            if orthologue:
                std, target_chembl_id, species = key
            else:
                std, target_chembl_id = key
                species = "Homo sapiens"
            row = {
                "smiles": std,
                "murcko_scaffold": std_to_scaffold.get(std),
                "contributing_molecule_chembl_ids": std_to_contrib.get(std, []),
                "target_chembl_id": target_chembl_id,
                "species_provenance": species,
                "pchembl_value": round(acc["best_pchembl"], 3) if acc["best_pchembl"] is not None else None,
                "best_relation": acc["best_relation"],
                "confidence_score": acc["best_confidence"],
                "is_censored": acc.get("best_is_censored", False),
                "is_measured_inactive": acc.get("is_measured_inactive", False),
                "n_documents": len(acc["documents"]),
                "document_chembl_ids": sorted(acc["documents"]),
                "assay_types_seen": sorted(acc["assay_types"]),
                "first_seen_year": acc.get("first_seen_year"),
            }
            f.write(json.dumps(row) + "\n")
    print(f"Wrote {out_path}", flush=True)


def main():
    print("Loading Stage 6 outputs + confidence backfill + single-protein allowlist + document years...", flush=True)
    std_cache = load_std_cache()
    fp_meta = load_fp_meta()
    confidence_map = load_confidence_map()
    single_protein_ids = load_single_protein_ids()
    year_map = load_document_year_map()

    print("\nAggregating native-human tier...", flush=True)
    native_pairs = aggregate_native_human(std_cache, confidence_map, single_protein_ids, year_map)
    write_output(native_pairs, fp_meta, os.path.join(DATA_DIR, "stage7_native_human_pairs.jsonl"))

    print("\nAggregating orthologue tier...", flush=True)
    ortho_pairs = aggregate_orthologue(std_cache, confidence_map, year_map)
    write_output(ortho_pairs, fp_meta, os.path.join(DATA_DIR, "stage7_orthologue_pairs.jsonl"), orthologue=True)


if __name__ == "__main__":
    main()
