"""
Classifies a dock_compound()/redock_reference() result's failure into one
of a small set of concrete categories, with a suggested action — the
existing status enum ("ligand_prep_failed"/"no_pose"/"ok"/
"engine_unavailable") is coarser than useful for a results table trying to
tell "this compound is chemically unusual" apart from "something about the
engine/receptor/box is wrong."

Built by pattern-matching the free-text error/reason strings the pipeline
already produces (see pipeline.py, ligand_prep.py, consensus.py) rather than
inventing detection for failure modes no code path here actually produces
— "other" is the honest fallback for anything unrecognised, not a made-up
category.
"""


def classify_failure(result: dict):
    """Returns {'category': str, 'suggested_action': str}, or None for a
       successful ('ok') result — nothing to diagnose there."""
    status = result.get("status")
    if status == "ok":
        return None

    error = (result.get("error") or "").lower()
    engine_errors = result.get("engine_errors") or {}
    reason = (result.get("reason") or "").lower()

    if status == "ligand_prep_failed":
        if "invalid smiles" in error:
            return {"category": "invalid_molecule",
                    "suggested_action": "RDKit could not parse this SMILES at all — check it's a valid structure."}
        if "3d embedding" in error:
            return {"category": "conformer_generation_failed",
                    "suggested_action": "Could not generate a 3D conformer for this molecule (often a strained or "
                                        "unusual ring system) — try a different tautomer or stereochemistry."}
        if "meeko" in error or "pdbqt" in error:
            return {"category": "ligand_conversion_failed",
                    "suggested_action": "Could not convert this structure to PDBQT (atom typing/charge assignment) "
                                        "— check for unusual atom types or a valence error."}
        return {"category": "ligand_prep_failed", "suggested_action": f"Ligand preparation failed: {result.get('error') or 'see error message'}."}

    if status == "engine_unavailable":
        return {"category": "engine_unavailable",
                "suggested_action": "AutoDock Vina isn't installed or wasn't found on PATH — check the Docking tab's status."}

    if status == "no_pose":
        if engine_errors:
            detail = "; ".join(f"{k}: {v}" for k, v in engine_errors.items())
            return {"category": "engine_error", "suggested_action": f"The docking engine reported an error: {detail}"}
        if "no physically valid poses" in reason:
            return {"category": "pose_validity_failed",
                    "suggested_action": "Vina produced poses, but none passed the physical-validity check (PoseBusters) "
                                        "— try a higher exhaustiveness, or confirm the binding box isn't clipping the ligand."}
        return {"category": "no_pose", "suggested_action": "No docking pose was produced — see the reason above for detail."}

    return {"category": "other", "suggested_action": "See the error message above for detail."}
