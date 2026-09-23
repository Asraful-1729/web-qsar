"""
A5 — "Research Story Generator": assembles the linear evidence chain the
upgrade proposal specifies (natural source -> chemical identity ->
reported activity -> target prediction -> QSAR -> docking -> interaction
analysis -> ADMET -> off-target analysis -> overall evidence summary)
for ONE already-docked compound, reusing every module built earlier in
this phase rather than re-implementing anything:
  - A1's plant_source (already threaded through run_metadata)
  - A2's target_fishing.search() for target prediction + off-target
  - A4's literature.search() for reported bioactivity (best-effort --
    this is the report's only network-dependent section, so a failure
    here degrades that ONE section rather than the whole report)
  - admet.py's admet_profile() for the ADMET section (already degrades
    gracefully to deterministic-only if the ADMET-AI worker isn't up)
  - serving.model_adapter's QSAR prediction, when the caller doesn't
    already have one (a Screen job computes QSAR itself; a plain
    Docking-tab job doesn't)
  - A6's run_metadata (captured at submit time) for the Methods section

No text in evidence_summary/methods_draft/results_draft is model-
generated -- every sentence is templated from real numbers already
present in `dock_row`/`run_metadata`/the sections above, matching this
app's no-fabrication stance throughout: this assembles a structured,
citable DRAFT a researcher edits, not a finished paper, exactly as the
proposal's "manuscript-ready computational evidence package" phrasing
implies (a package to build on, not to publish verbatim).
"""
import datetime

PIPELINE_STAGES = [
    "Natural source", "Chemical identity", "Reported activity (literature)",
    "Target prediction", "QSAR prediction", "Docking & binding",
    "Interaction analysis", "ADMET profile", "Off-target analysis",
    "Overall evidence summary",
]


def _chemical_identity(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"available": False, "smiles": smiles, "note": "Could not parse this SMILES."}
    return {
        "available": True,
        "smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "n_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "n_hba": rdMolDescriptors.CalcNumHBA(mol),
        "n_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
    }


def _reported_activity(query):
    if not query:
        return {"available": False, "note": "No plant source or search term given for this compound — literature search skipped."}
    import literature as LIT
    try:
        r = LIT.search(query, max_results=5)
        return {"available": True, "query": query, "n_results": r["n_results"], "papers": r["papers"]}
    except LIT.LiteratureError as e:
        return {"available": False, "query": query, "note": f"Literature search failed: {e}"}
    except ValueError as e:
        return {"available": False, "query": query, "note": str(e)}


def _target_and_off_target(smiles, target_id, threshold=0.4):
    import target_fishing as TF
    if not TF.available():
        na = {"available": False, "note": "Target-fishing index not available in this build."}
        return na, na
    target_chembl = target_id.split("_", 1)[0] if target_id else None
    try:
        r = TF.search(smiles, threshold=threshold)
    except ValueError as e:
        na = {"available": False, "note": str(e)}
        return na, na
    on_target = next((h for h in r["results"] if h["target_chembl"] == target_chembl), None)
    off_targets = [h for h in r["results"] if h is not on_target]
    return (
        {"available": True, "on_target_supported": on_target is not None, "hit": on_target,
         "n_targets_searched": r["n_targets_searched"]},
        {"available": True, "n_off_targets": len(off_targets), "hits": off_targets[:10]},
    )


def _qsar(target_id, smiles, precomputed=None):
    """precomputed lets a Screen job (which already ran QSAR for this
       exact compound) pass that result straight through instead of
       re-predicting — see build_report's docstring."""
    if precomputed is not None:
        return {"available": True, **precomputed}
    from serving import model_adapter as MA
    try:
        target = MA.load_target(target_id)
    except Exception:
        return {"available": False, "note": "No QSAR model available for this target."}
    try:
        df = target.predict_smiles([smiles])
    except Exception as e:
        return {"available": False, "note": f"QSAR prediction failed: {e}"}
    if df.empty or not bool(df.iloc[0].get("Parsed_OK", False)):
        return {"available": False, "note": "Could not parse this SMILES for QSAR."}
    row = df.iloc[0]
    return {
        "available": True,
        "predicted_pIC50": row.get("Predicted_pIC50"),
        "in_domain": bool(row.get("In_AD")),
        "model": target.metrics.get("Best_Model") or target.metrics.get("best_model"),
        "test_r2": target.metrics.get("R2_Test"),
        "test_rmse": target.metrics.get("RMSE_Test"),
    }


def _admet(smiles):
    import admet as ADMET
    try:
        return {"available": True, "profile": ADMET.admet_profile(smiles)}
    except Exception as e:
        return {"available": False, "note": str(e)}


def _evidence_summary(target_id, docking, qsar, target_pred, off_target, admet, reported_activity):
    """Plain-language rollup, built only from fields that are actually
       present -- never invents a number that wasn't computed."""
    lines = []
    vina = docking.get("vina_score")
    if vina is not None:
        lines.append(f"Docking predicts binding to {target_id} with a Vina score of {vina:.2f} kcal/mol.")
    else:
        lines.append(f"Docking did not produce a usable pose for this compound against {target_id}.")

    if qsar.get("available") and qsar.get("predicted_pIC50") is not None:
        dom = "within" if qsar.get("in_domain") else "outside"
        lines.append(
            f"The QSAR model ({qsar.get('model') or 'best available'}, test R² = {qsar.get('test_r2')}) "
            f"predicts a pIC50 of {qsar['predicted_pIC50']:.2f}, {dom} the model's applicability domain."
        )

    if target_pred.get("available"):
        if target_pred.get("on_target_supported"):
            hit = target_pred["hit"]
            lines.append(
                f"Independent ligand-based evidence supports this target: {hit['n_similar_actives']} "
                f"structurally similar known active(s) exist for {target_id} (best Tanimoto {hit['best_similarity']})."
            )
        else:
            lines.append(f"No structurally similar known actives were found for {target_id} in the curated bioactivity data.")

    if off_target.get("available") and off_target.get("n_off_targets"):
        names = ", ".join(h["target_id"] for h in off_target["hits"][:5])
        lines.append(f"Similarity-based off-target signal was also found for {off_target['n_off_targets']} other target(s): {names}.")

    if admet.get("available"):
        prof = admet["profile"]
        flags = prof.get("drug_likeness_flags") or {}
        n_alerts = prof.get("n_alerts")
        pass_txt = "passes" if flags.get("lipinski_pass") else "does not pass"
        lines.append(f"ADMET profile: {pass_txt} Lipinski's rule of five" + (f", {n_alerts} structural alert(s) flagged." if n_alerts is not None else "."))

    if reported_activity.get("available") and reported_activity.get("n_results"):
        lines.append(f"{reported_activity['n_results']} related paper(s) were found in PubMed for \"{reported_activity['query']}\".")

    return " ".join(lines)


def _methods_draft(target_id, run_metadata, qsar):
    rm = run_metadata or {}
    mode = "site-specific" if rm.get("docking_mode") == "site_specific" else "blind (whole-protein)"
    parts = [
        f"Molecular docking against {target_id} was performed using {rm.get('engine') or 'AutoDock Vina'} "
        f"(exhaustiveness = {rm.get('exhaustiveness')}, {rm.get('n_poses')} poses per compound, {mode} search).",
    ]
    if rm.get("pdb_source"):
        parts.append(f"The receptor structure was sourced from {rm['pdb_source']}.")
    if qsar.get("available") and qsar.get("model"):
        parts.append(
            f"Bioactivity was additionally predicted with a {qsar['model']} QSAR model "
            f"(test R² = {qsar.get('test_r2')}, test RMSE = {qsar.get('test_rmse')})."
        )
    versions = rm.get("software_versions") or {}
    if versions:
        vtxt = ", ".join(f"{k} {v}" for k, v in versions.items() if v)
        if vtxt:
            parts.append(f"Software versions: {vtxt}.")
    return " ".join(parts)


def build_report(target_id, smiles, dock_row, run_metadata=None, plant_source=None, qsar_precomputed=None,
                  literature_query=None, include_literature=True):
    """dock_row: one DockResultRow dict (already computed, from a
       finished Docking or Screen job) -- this function never re-runs
       docking, it only augments an existing result. qsar_precomputed:
       a Screen job's own `row["qsar"]` dict, if available."""
    chem = _chemical_identity(smiles)
    reported = (
        _reported_activity(literature_query or plant_source or target_id)
        if include_literature else {"available": False, "note": "Literature search skipped for this report."}
    )
    target_pred, off_target = _target_and_off_target(smiles, target_id)
    qsar = _qsar(target_id, smiles, precomputed=qsar_precomputed)
    admet = _admet(smiles)
    docking = dock_row or {}

    summary = _evidence_summary(target_id, docking, qsar, target_pred, off_target, admet, reported)
    methods = _methods_draft(target_id, run_metadata, qsar)

    return {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "target_id": target_id,
        "pipeline_stages": PIPELINE_STAGES,
        "natural_source": {"plant_source": plant_source},
        "chemical_identity": chem,
        "reported_activity": reported,
        "target_prediction": target_pred,
        "qsar_prediction": qsar,
        "docking": docking,
        "interaction_analysis": {"interactions": docking.get("interactions") or [], "residue_overlap_pct": docking.get("residue_overlap_pct")},
        "admet": admet,
        "off_target_analysis": off_target,
        "evidence_summary": summary,
        "methods_draft": methods,
        "reproducibility": run_metadata or {},
        "disclaimer": "Computational predictions only — not validated experimentally. Treat as a hypothesis-generation aid, not a diagnostic or clinical claim.",
    }


def to_markdown(report):
    """Renders the structured report dict as a single Markdown document —
       the downloadable form of the same data the UI renders as sections."""
    r = report
    lines = [f"# Computational Evidence Report — {r['target_id']}", ""]
    lines.append(f"_Generated {r['generated_at']}_")
    lines.append("")
    lines.append(f"> {r['disclaimer']}")
    lines.append("")
    lines.append("## Pipeline")
    lines.append(" → ".join(r["pipeline_stages"]))
    lines.append("")

    lines.append("## 1. Natural source")
    lines.append(r["natural_source"].get("plant_source") or "_Not specified._")
    lines.append("")

    lines.append("## 2. Chemical identity")
    ci = r["chemical_identity"]
    if ci.get("available"):
        lines.append(f"- SMILES: `{ci['smiles']}`")
        lines.append(f"- Molecular formula: {ci['molecular_formula']}")
        lines.append(f"- Molecular weight: {ci['molecular_weight']}")
        lines.append(f"- LogP: {ci['logp']}, TPSA: {ci['tpsa']}, HBD/HBA: {ci['n_hbd']}/{ci['n_hba']}, rotatable bonds: {ci['n_rotatable_bonds']}")
    else:
        lines.append(f"_{ci.get('note', 'Unavailable.')}_")
    lines.append("")

    lines.append("## 3. Reported activity (literature)")
    ra = r["reported_activity"]
    if ra.get("available"):
        for p in ra["papers"]:
            lines.append(f"- [{p['title']}]({p.get('url') or '#'}) — {p.get('authors', '')} ({p.get('year', 'n.d.')}), {p.get('journal', '')}")
    else:
        lines.append(f"_{ra.get('note', 'Unavailable.')}_")
    lines.append("")

    lines.append("## 4. Target prediction")
    tp = r["target_prediction"]
    lines.append(
        f"On-target similarity evidence: {'supported' if tp.get('on_target_supported') else 'not found'}."
        if tp.get("available") else f"_{tp.get('note', 'Unavailable.')}_"
    )
    lines.append("")

    lines.append("## 5. QSAR prediction")
    qs = r["qsar_prediction"]
    if qs.get("available") and qs.get("predicted_pIC50") is not None:
        lines.append(f"- Predicted pIC50: {qs['predicted_pIC50']:.3f} ({'in' if qs.get('in_domain') else 'out of'} domain)")
        lines.append(f"- Model: {qs.get('model')} (test R² = {qs.get('test_r2')}, RMSE = {qs.get('test_rmse')})")
    else:
        lines.append(f"_{qs.get('note', 'Unavailable.')}_")
    lines.append("")

    lines.append("## 6. Docking & binding")
    dk = r["docking"]
    if dk.get("vina_score") is not None:
        lines.append(f"- Vina score: {dk['vina_score']:.2f} kcal/mol")
        lines.append(f"- Confidence: {dk.get('confidence')}")
    else:
        lines.append(f"_{dk.get('reason') or dk.get('error') or 'No docking result.'}_")
    lines.append("")

    lines.append("## 7. Interaction analysis")
    ia = r["interaction_analysis"]["interactions"]
    if ia:
        for i in ia[:15]:
            label = i.get("label") or i.get("name") or i.get("type") or "interaction"
            res = f"{i.get('resname', '')}{i.get('resid', '')}".strip() or i.get("residue", "")
            lines.append(f"- {label} — {res}")
    else:
        lines.append("_No interaction data available._")
    lines.append("")

    lines.append("## 8. ADMET profile")
    ad = r["admet"]
    if ad.get("available"):
        prof = ad["profile"]
        pc = prof.get("physicochemical") or {}
        flags = prof.get("drug_likeness_flags") or {}
        lines.append(f"- MW {pc.get('mw')}, LogP {pc.get('logp')}, QED {pc.get('qed')}")
        lines.append(f"- Lipinski: {'pass' if flags.get('lipinski_pass') else 'fail'} ({flags.get('lipinski_violations')} violations), {prof.get('n_alerts')} structural alert(s)")
    else:
        lines.append(f"_{ad.get('note', 'Unavailable.')}_")
    lines.append("")

    lines.append("## 9. Off-target analysis")
    ot = r["off_target_analysis"]
    if ot.get("available") and ot.get("hits"):
        for h in ot["hits"]:
            lines.append(f"- {h['target_id']} — best similarity {h['best_similarity']}, {h['n_similar_actives']} similar active(s)")
    elif ot.get("available"):
        lines.append("_No off-target similarity signal found._")
    else:
        lines.append(f"_{ot.get('note', 'Unavailable.')}_")
    lines.append("")

    lines.append("## Overall evidence summary")
    lines.append(r["evidence_summary"])
    lines.append("")

    lines.append("## Methods (draft)")
    lines.append(r["methods_draft"])
    lines.append("")

    lines.append("## Reproducibility")
    rp = r["reproducibility"]
    for k in ("target_id", "plant_source", "timestamp", "docking_mode", "engine", "exhaustiveness", "n_poses", "pdb_source"):
        if rp.get(k) is not None:
            lines.append(f"- {k}: {rp[k]}")
    lines.append("")

    return "\n".join(lines)
