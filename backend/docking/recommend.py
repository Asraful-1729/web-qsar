"""Disease -> ranked target -> 'why this structure' evidence.

Backs three read-only API endpoints (see app.py): a disease list, the
targets associated with a chosen disease (ranked by disease-association
score), and the evidence bundle behind one target's recommended docking
structure.

Scoped to the targets in our QSAR panel (the CHEMBL-prefixed entries in
docking_registry.json, ~50 today, regardless of whether their model bucket
is downloaded onto this machine yet) — panel_results_v2.csv covers ~669
genes total, most of which have no trained potency model and no docking
receptor at all; extending this beyond our QSAR panel is future work, not
this phase.

"panel_evidence" is panel_results_v2.csv's crystallographic context
(resolution/RSCC/RSR/rank) for the matching gene — supporting context, not
proof; see scripts/panel_candidates.py for why the CSV's own top pick can't
be trusted blindly (e.g. PDGFRB's best-ranked structure turned out to be a
bare sugar, not an inhibitor).
"""
import functools
import os

import pandas as pd

from scripts.batch_validate import gene_for_target, usable_targets
from scripts.panel_candidates import panel_evidence, PANEL_CSV
from docking.profile import load_registry
from serving import model_adapter as MA


@functools.lru_cache(maxsize=1)
def _panel_df():
    if not os.path.exists(PANEL_CSV):
        return None
    return pd.read_csv(PANEL_CSV)


def _target_gene_map():
    """(target_id -> gene, gene -> target_id) for every QSAR-panel target,
       i.e. every CHEMBL-prefixed entry in docking_registry.json (~50
       today) — NOT limited to whichever model buckets happen to be
       downloaded onto this machine right now. docking_registry.json is
       committed, portable metadata (unlike models/<target_id>/, which is
       fetched on demand from R2 by the Downloads tab), so a target we
       genuinely have a QSAR model for must still resolve to its real
       target_id and show up as has_qsar_model even before its bucket is
       downloaded — see targets_for_disease()'s model_installed field for
       the actually-on-this-disk answer, which the download-gate UI uses
       to decide whether to fetch it first."""
    reg = load_registry()
    t2g = {tid: gene_for_target(tid) for tid in reg if tid.startswith("CHEMBL")}
    g2t = {g: t for t, g in t2g.items()}
    return t2g, g2t


def list_diseases():
    """Every disease in the panel CSV, sorted by name — NOT restricted to
       diseases tied to a target we happen to have a local QSAR model
       for right now. That restriction made sense for the old static web
       deployment (models/ never changed while the server ran, so "ours"
       was a fixed, complete set); it's wrong for the desktop app's
       download-on-demand workflow, where a disease not yet showing a
       locally-downloaded QSAR target is still real and worth browsing —
       targets_for_disease() already shows the full target landscape
       (QSAR-modeled or docking-only) for whatever disease gets picked
       here, has_qsar_model correctly reflecting what's actually
       downloaded."""
    df = _panel_df()
    if df is None:
        return []
    out = (df[["diseaseId", "diseaseName", "is_therapeutic_area"]]
           .drop_duplicates(subset=["diseaseId"])
           .sort_values("diseaseName"))
    return [
        {"disease_id": r.diseaseId, "name": r.diseaseName, "is_therapeutic_area": bool(r.is_therapeutic_area)}
        for r in out.itertuples()
    ]


def targets_for_disease(disease_id):
    """EVERY protein target associated with disease_id in the Version 2 CSV
       (targetSymbol column), ranked by disease-association score (desc) —
       not just the subset we happen to have a trained QSAR model for. Each
       row carries has_qsar_model so the UI can show the full disease-target
       landscape while making clear which ones can actually be predicted/
       docked against; QSAR fields are only populated when a model exists
       (never invented for the rest)."""
    df = _panel_df()
    if df is None:
        return []
    _, g2t = _target_gene_map()
    rows = df[df["diseaseId"] == disease_id].sort_values("score", ascending=False)
    meta_by_id = {m["target_id"]: m["metrics"] for m in MA.list_targets_meta()}
    installed_ids = set(usable_targets())

    out = []
    seen = set()
    for r in rows.itertuples():
        symbol = r.targetSymbol
        if symbol in seen:   # a gene can appear more than once per disease in the CSV; keep the first (highest-scored, already sorted)
            continue
        seen.add(symbol)
        tid = g2t.get(symbol)
        has_model = tid is not None
        metrics = (meta_by_id.get(tid) or {}) if has_model else {}
        out.append({
            "target_id": tid,
            "target_symbol": symbol,
            "has_qsar_model": has_model,
            # Whether the model bucket is actually downloaded onto THIS
            # machine right now — has_qsar_model above is the download-
            # independent "do we have a model for this gene at all"
            # answer; the download-gate UI (useDownloadGate.ts) uses this
            # one to decide whether picking this target needs to fetch it
            # first.
            "model_installed": has_model and tid in installed_ids,
            "disease_score": round(float(r.score), 4),
            "test_r2": metrics.get("R2_Test"),
            "test_rmse": metrics.get("RMSE_Test"),
        })
    return out


def recommendation(target_id):
    """The full 'Why this?' evidence bundle for one target's recommended
       docking structure, or None if target_id isn't one of ours."""
    t2g, _ = _target_gene_map()
    gene = t2g.get(target_id)
    if gene is None:
        return None

    reg = load_registry()
    reg_entry = reg.get(target_id) or {}
    panel = panel_evidence(gene) or {}

    our_pdb_id = None
    src = reg_entry.get("pdb_source")
    if src:
        our_pdb_id = src.split("_raw")[0].split(".")[0].upper()

    rank = None
    panel_best = panel.get("best_pdb_id")
    if our_pdb_id and isinstance(panel.get("n_qualifying_structures"), (int, float)):
        # rank position of the structure we use, among the panel's own
        # crystallographic ranking (1-indexed), if it's in there
        df = _panel_df()
        if df is not None:
            match = df[df["targetSymbol"] == gene]
            if not match.empty:
                ranked = str(match.iloc[0].get("all_pdb_ids_ranked") or "").split(";")
                if our_pdb_id in ranked:
                    rank = ranked.index(our_pdb_id) + 1

    headline = (f"Recommended structure: {our_pdb_id}" if our_pdb_id
               else "No receptor structure prepared yet for this target")

    return {
        "target_id": target_id,
        "gene": gene,
        "headline": headline,
        "structure": {
            "pdb_id": our_pdb_id,
            "ligand_resname": reg_entry.get("reference_ligand_resname"),
            "rank_in_panel_evidence": rank,
        },
        "panel_evidence": {
            "top_ranked_pdb_id": panel_best,
            "top_ranked_chain": panel.get("best_chain"),
            "top_ranked_ligand": panel.get("ligand"),
            "resolution": panel.get("resolution"),
            "resolution_tier": panel.get("resolution_tier"),
            "r_free": panel.get("r_free"),
            "modeled_residues": panel.get("modeled_residues"),
            "ligand_RSCC": panel.get("ligand_RSCC"),
            "ligand_RSR": panel.get("ligand_RSR"),
            "n_qualifying_structures": panel.get("n_qualifying_structures"),
            "chembl_activity_records": panel.get("activity_records"),
            "note": "Crystallographic quality context only, not independent proof of pose accuracy.",
        },
    }
