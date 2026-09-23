"""
Ligand-selection audit for every target in docking_registry.json.

Answers the question "are the recorded reference ligands actually the
intended crystallographic ligands?" using ONLY already-committed metadata
(docking_registry.json + panel_results_v2.csv) — no downloads, no network
calls, safe to re-run any time.

For any target whose raw/cleaned structure files ALSO happen to already be
present under docking_targets/<id>/ (e.g. downloaded during this session),
the report additionally computes REAL heavy-atom count / MW / pocket-
residue count directly from the files instead of approximating from the
CSV text — those rows are marked measurement_source="file", everything
else measurement_source="csv" or "unavailable".

Usage:
    python -m scripts.audit_ligand_selection [--out PATH]

Writes a CSV report (default: ligand_audit_report.csv at the repo root)
with one row per registry target, and prints a summary of flagged rows.
"""
import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.batch_validate import gene_for_target
from scripts.select_receptor import ADDITIVE_BLACKLIST, MIN_LIGAND_MW
from scripts.panel_candidates import _panel_df, _parse_top5, _LIG_RE
from docking.profile import REGISTRY
from docking.receptor_prep import MAX_LIGAND_MW, MIN_HEAVY_ATOMS, MIN_POCKET_RESIDUES

DOCKING_TARGETS_DIR = os.environ.get("DOCKING_TARGETS_DIR", "docking_targets")


def _pdb_id_from_source(pdb_source):
    if not pdb_source:
        return None
    return pdb_source.split("_raw")[0].split(".")[0].upper()


def _csv_corroboration(gene, pdb_id, resname, resname_full=None):
    """Cross-references (pdb_id, resname) against the panel CSV's own
       ranking for this gene. Returns (method, csv_mw_or_None, csv_rank).
       method is one of:
         csv_top5_match       — exact (pdb_id, resname) pair found in the
                                 CSV's own top5-with-real-ligand list
         csv_pdb_different_ligand — this pdb_id IS in the CSV's top5, but
                                 with a DIFFERENT resname recorded there
         csv_ranked_only      — pdb_id appears in all_pdb_ids_ranked but
                                 not in the (additive/MW-filtered) top5
         not_in_csv           — pdb_id doesn't appear in this gene's
                                 ranked list at all
         no_csv_row           — gene isn't in the panel CSV at all

       resname_full: the registry's reference_ligand_full_id, if set (see
       pdb_fetch.py's resolve_ligand_resname / receptor_prep.py's
       build_receptor) — the CSV records the TRUE, untruncated CCD code
       (e.g. "A1AWR"), but `resname` itself is deliberately the truncated,
       file-consistent form (e.g. "A1A" — what a live lookup against the
       raw structure must search for). Checking resname_full too, when
       present, is what tells apart a genuinely different ligand from a
       cosmetic truncation-of-the-same-one (confirmed on GENE_LRRK2/9C76)."""
    df = _panel_df()
    if df is None or gene not in df.index:
        return "no_csv_row", None, None
    row = df.loc[gene]
    ranked_ids = [p for p in str(row.get("all_pdb_ids_ranked") or "").split(";") if p]
    try:
        csv_rank = ranked_ids.index(pdb_id) + 1 if pdb_id else None
    except ValueError:
        csv_rank = None

    top5 = _parse_top5(row.get("top5_pdb_summary"))
    top5_by_pdb = {c["pdb_id"]: c["resname"] for c in top5}

    # MW lookup: search the RAW top5 text (not just the pre-filtered "best"
    # pick) for this exact pdb_id's chunk, then this exact resname within
    # it — finds the MW even when the registry's ligand wasn't the
    # heaviest candidate in that structure.
    csv_mw = None
    summary = row.get("top5_pdb_summary")
    if isinstance(summary, str) and pdb_id:
        for chunk in summary.split(" | "):
            if not chunk.strip().startswith(pdb_id):
                continue
            for want in (resname, resname_full):
                if not want:
                    continue
                for m in _LIG_RE.finditer(chunk):
                    if m.group(1).upper() == want.upper():
                        csv_mw = float(m.group(2))
                        break
                if csv_mw is not None:
                    break

    if pdb_id in top5_by_pdb:
        csv_resname = top5_by_pdb[pdb_id]
        same_ligand = csv_resname == resname or (resname_full and csv_resname == resname_full.upper())
        method = "csv_top5_match" if same_ligand else "csv_pdb_different_ligand"
    elif csv_rank is not None:
        method = "csv_ranked_only"
    else:
        method = "not_in_csv"
    return method, csv_mw, csv_rank


def _measure_from_files(target_id, pdb_source, resname, chain=None, ligand_chain=None):
    """If this target's raw+cleaned structures are already on disk, compute
       REAL heavy-atom-count/MW/pocket-residue-count directly — returns
       None if either file is missing (never triggers a download).

       chain MUST be passed through to extract_reference_ligand(): many
       real depositions carry 2+ copies of the reference ligand across
       different chains (e.g. a homomeric assembly — confirmed on
       GENE_CAMK2A/7REC, 6 copies of ligand 7ZV across chains B-G), and
       receptor_prep.py's build_receptor() always restricts to the ONE
       chain actually used for the receptor (see detect_chain.py). Without
       this, extract_reference_ligand's own "largest matching-name group"
       tie-break can pick a DIFFERENT chain's copy than the one the
       receptor was built around, measuring that copy's distance to a
       receptor it was never near — producing a spurious
       too_few_pocket_contacts flag for an otherwise entirely correct,
       already-validated target (confirmed: GENE_CAMK2A's own build-time
       binding_site_residues has 21 contacts; the unrestricted re-measurement
       here found 0 before this fix).

       ligand_chain: the ligand's OWN HETATM chain label, when different
       from `chain` (which chain the receptor was stripped to) — see
       receptor_prep.build_receptor's ligand_chain docstring. Without this,
       a target like GENE_UBA2 (ligand VAY filed under chain C, receptor
       built from chain B) would search chain B for a ligand that isn't
       filed there at all and wrongly report no_reference_ligand_recorded /
       measurement_error for an already-correct, already-fixed target."""
    raw_path = os.path.join(DOCKING_TARGETS_DIR, target_id, pdb_source or "")
    clean_path = os.path.join(DOCKING_TARGETS_DIR, target_id, "receptor_clean.pdb")
    if not (pdb_source and os.path.exists(raw_path) and os.path.exists(clean_path)):
        return None
    try:
        from docking.receptor_prep import extract_reference_ligand, pocket_residues
        coords, found_name, n_atoms, mw = extract_reference_ligand(raw_path, ref_resname=resname, chain=ligand_chain or chain)
        n_pocket = len(pocket_residues(clean_path, coords, cutoff=5.0))
        return {"n_atoms": n_atoms, "mw": round(mw, 1), "n_pocket_residues": n_pocket}
    except Exception as e:
        return {"error": str(e)}


def audit():
    with open(REGISTRY) as f:
        data = json.load(f)
    targets = data.get("targets", data if isinstance(data, list) else [])

    rows = []
    for t in targets:
        target_id = t["target_id"]
        pdb_source = t.get("pdb_source")
        pdb_id = _pdb_id_from_source(pdb_source)
        resname = (t.get("reference_ligand_resname") or "").upper()
        resname_full = t.get("reference_ligand_full_id")
        gene = gene_for_target(target_id)

        # A target with no recorded resname (site_source == "none_validated"
        # — deliberately blind-only, see docking/receptor_prep.py's
        # no-ligand fallback) has nothing to corroborate against the CSV or
        # re-measure from files: resname="" would otherwise fall through
        # extract_reference_ligand's OWN "no ref_resname given -> pick the
        # largest non-additive HETATM group automatically" behavior and
        # silently measure some UNRELATED stray HETATM group near that
        # chain as if it were "the ligand", producing nonsense
        # mw_below_floor/too_few_pocket_contacts/csv_pdb_different_ligand
        # flags for a target that was already correctly handled. Skip
        # straight to no_reference_ligand_recorded instead.
        if resname:
            method, csv_mw, csv_rank = _csv_corroboration(gene, pdb_id, resname, resname_full=resname_full)
            measured = _measure_from_files(target_id, pdb_source, resname, chain=t.get("chain"), ligand_chain=t.get("ligand_chain"))
        else:
            method, csv_mw, csv_rank = "no_ligand_recorded", None, None
            measured = None

        n_atoms = measured.get("n_atoms") if measured and "n_atoms" in measured else None
        mw = measured.get("mw") if measured and "mw" in measured else csv_mw
        n_pocket = measured.get("n_pocket_residues") if measured and "n_pocket_residues" in measured else None
        measurement_source = "file" if (measured and "n_atoms" in measured) else ("csv" if csv_mw is not None else "unavailable")

        flags = []
        if resname in ADDITIVE_BLACKLIST:
            flags.append("resname_in_additive_blacklist")
        if not resname:
            flags.append("no_reference_ligand_recorded")
        else:
            if n_atoms is not None and n_atoms < MIN_HEAVY_ATOMS:
                flags.append(f"too_few_heavy_atoms({n_atoms})")
            if mw is not None and mw < MIN_LIGAND_MW:
                flags.append(f"mw_below_floor({mw:.0f}Da)")
            if mw is not None and mw > MAX_LIGAND_MW:
                flags.append(f"mw_above_ceiling({mw:.0f}Da)")
            if n_pocket is not None and n_pocket < MIN_POCKET_RESIDUES:
                flags.append(f"too_few_pocket_contacts({n_pocket})")
            if method in ("not_in_csv", "csv_pdb_different_ligand"):
                flags.append(method)
            if measured and "error" in measured:
                flags.append(f"measurement_error: {measured['error']}")

        rows.append({
            "target_id": target_id,
            "gene": gene,
            "pdb_id": pdb_id or "",
            "resname": resname,
            "ligand_heavy_atoms": n_atoms if n_atoms is not None else "",
            "ligand_mw_da": mw if mw is not None else "",
            "measurement_source": measurement_source,
            "n_pocket_residues_5A": n_pocket if n_pocket is not None else "",
            "box_center": " ".join(f"{v:.2f}" for v in (t.get("center") or [])),
            "box_size": " ".join(f"{v:.2f}" for v in (t.get("box_size") or [])),
            "selection_method": method,
            "csv_rank": csv_rank if csv_rank is not None else "",
            "flags": ";".join(flags),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="ligand_audit_report.csv")
    args = ap.parse_args()

    rows = audit()
    fieldnames = ["target_id", "gene", "pdb_id", "resname", "ligand_heavy_atoms", "ligand_mw_da",
                  "measurement_source", "n_pocket_residues_5A", "box_center", "box_size",
                  "selection_method", "csv_rank", "flags"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    flagged = [r for r in rows if r["flags"]]
    measured_from_file = sum(1 for r in rows if r["measurement_source"] == "file")
    print(f"{len(rows)} targets audited, report written to {args.out}")
    print(f"{measured_from_file} target(s) had files on disk — measured directly, not CSV-approximated")
    print(f"{len(flagged)} target(s) flagged:")
    for r in flagged:
        print(f"  {r['target_id']:35s} {r['flags']}")


if __name__ == "__main__":
    main()
