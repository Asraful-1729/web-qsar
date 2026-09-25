"""
Per-compound docking orchestration + per-target validation.
[orchestration WRITTEN; gated helpers TESTED; subprocess runs UNVALIDATED]

  dock_compound:
    ligand prep -> Vina dock -> PoseBusters gate -> best valid pose
                -> GNINA CNN rescore (second opinion)
                -> confidence (validity + self-consistency + CNNscore)
                -> optional 2D interaction diagram + reference-residue overlap
  redock_reference: re-dock the target's known ligand, RMSD to the crystal pose.
"""
import os
from rdkit import Chem

from .ligand_prep import prepare_ligand
from .validity import ValidityGate
from .consensus import select_pose, assign_confidence
from .engines import VinaEngine, GninaRescorer
from .rmsd import safe_rmsd
from .failure_diagnostics import classify_failure


def pose_pdb_with_hydrogens(pose_mol):
    """PDB block of the docked pose for 3D (ball-and-stick) display, with every
       hydrogen present — including polar ones like NH2/OH.

       VinaEngine.dock() hands back heavy-atom-only poses (Chem.RemoveHs, since
       that's what the interaction detector/GNINA rescoring already validated
       against), so there is no docked H position to reuse here. Chem.AddHs(...,
       addCoords=True) fills every hydrogen back in with standard bond geometry
       from the docked heavy-atom positions — the same approach visualization
       tools (PyMOL/Chimera 'add hydrogens') use on a heavy-atom structure.
       This is a geometry completion for display, not Vina's own optimised H
       placement — polar-H orientation here is a reasonable estimate, not the
       exact rotamer Vina searched."""
    vis_mol = Chem.AddHs(pose_mol, addCoords=True)
    return Chem.MolToPDBBlock(vis_mol)


def dock_compound(profile, smiles, engine=None, rescorer=None, n_poses=9,
                  rmsd_threshold=2.0, make_diagram=False, reference_interactions=None):
    lig = prepare_ligand(smiles)
    if not lig.ok:
        result = {"smiles": smiles, "status": "ligand_prep_failed", "error": lig.error}
        diag = classify_failure(result)
        if diag:
            result.update(diag)
        return result
    engine = engine or VinaEngine()
    errors = {}
    try:
        vina_poses = engine.dock(profile, lig, n_poses=n_poses) if engine.available() else []
    except Exception as e:
        vina_poses = []; errors["vina"] = str(e)

    gate = ValidityGate(receptor_pdb=profile.get("receptor_pdb"))
    sel = select_pose(vina_poses, gate, rmsd_threshold)

    result = {"smiles": smiles, "target": profile.get("target_id"),
              "n_valid": sel["n_valid"], "pose_self_consistency": sel.get("pose_self_consistency", 0),
              "vina_score": sel.get("vina_score"), "consensus_pose": sel.get("consensus_pose"),
              "engine_errors": errors, "reason": sel["reason"]}

    best = sel.get("best")
    gnina = None
    if best is not None:
        # GNINA CNN rescoring (second opinion)
        rescorer = rescorer or GninaRescorer()
        if rescorer.available() and profile.get("receptor_pdb"):
            gnina = rescorer.rescore(profile["receptor_pdb"], best.mol)
            result["gnina"] = gnina
        # interaction profiling + optional 2D diagram: PLIP (gold-standard
        # typing) when installed, else the built-in distance-based detector;
        # rendered LigPlot+-style (H-bond/pi-stacking/hydrophobic/salt-bridge/
        # halogen, polar/non-polar residue halos, legend — see CLAUDE.md §8).
        if profile.get("receptor_pdb") and os.path.exists(profile["receptor_pdb"]):
            try:
                from . import interaction_diagram as ID
                from . import interactions as I
                inter, source = ID.detect_interactions(profile["receptor_pdb"], best.mol)
                result["interactions"] = inter
                result["interaction_source"] = source
                if reference_interactions is not None:
                    result["residue_overlap_pct"] = I.residue_overlap(inter, reference_interactions)
                if make_diagram:
                    ref_res = {h["residue"] for h in (reference_interactions or [])}
                    result["interaction_png"] = ID.diagram_png(
                        best.mol, inter, title=f"{smiles[:30]}", source=source, ref_residues=ref_res)
            except Exception as e:
                result["interaction_error"] = str(e)
        if make_diagram:
            try:
                result["pose_pdb"] = pose_pdb_with_hydrogens(best.mol)
            except Exception as e:
                result["pose_pdb_error"] = str(e)

    result["confidence"] = assign_confidence(sel, gnina)
    result["status"] = "ok" if sel.get("consensus_pose") else "no_pose"
    diag = classify_failure(result)
    if diag:
        result.update(diag)
    return result


def redock_reference(profile, reference_smiles, crystal_sdf=None, engine=None, rmsd_threshold=2.0, n_poses=10,
                     seed=0xf00d):
    """VALIDATION: re-dock the reference ligand; if a crystal SDF is given, compute
       RMSD of the top pose to the crystal pose (atom-map-safe).

       seed is fixed (not None) by default — unlike ordinary compound
       docking, a "standardized" redocking-validation RMSD needs to
       reproduce the same result for the same PDB/ligand/parameters on
       every run, not a fresh random Vina search each time."""
    engine = engine or VinaEngine()
    if not engine.available():
        return {"status": "engine_unavailable", "engine": engine.name}
    lig = prepare_ligand(reference_smiles)
    if not lig.ok:
        return {"status": "ligand_prep_failed", "error": lig.error}
    poses = engine.dock(profile, lig, n_poses=n_poses, seed=seed)
    if not poses:
        return {"status": "no_pose"}
    top = min(poses, key=lambda p: p.score)
    out = {"status": "ok", "engine": engine.name, "top_score": round(top.score, 2), "n_poses": len(poses)}
    # The redocked pose itself — same H-completion treatment dock_compound's
    # pose_pdb uses, so an overlay view can show it alongside the crystal
    # pose instead of the caller only ever seeing the RMSD *number*.
    try:
        out["redocked_pose_pdb"] = pose_pdb_with_hydrogens(top.mol)
    except Exception:
        pass  # non-fatal — the RMSD/validated result below is still meaningful without it
    if crystal_sdf and os.path.exists(crystal_sdf):
        crystal = next((m for m in Chem.SDMolSupplier(crystal_sdf, removeHs=True) if m), None)
        if crystal is not None:
            try:
                out["reference_rmsd"] = round(safe_rmsd(top.mol, crystal), 3)
                out["validated"] = out["reference_rmsd"] < rmsd_threshold
            except ValueError as e:
                out["rmsd_error"] = str(e)
            # The experimental (crystal) pose itself, verbatim — lets a
            # caller render an "experimental vs redocked" 3D overlay
            # alongside the RMSD number, not just report the number alone.
            with open(crystal_sdf) as f:
                out["crystal_pose_sdf"] = f.read()
    return out


def redock_reference_for_profile(profile, exhaustiveness=8, n_poses=10):
    """Runs automatically alongside a docking/screen submission, using the
       SAME PDB structure and docking parameters (exhaustiveness/n_poses)
       the caller already resolved for their compounds — there is no
       separate "validate this structure" action any more; whatever
       structure/settings a user picks for docking is what gets validated.

       Forces --cpu 1 (see redock_reference/VinaEngine.dock's own notes):
       a fixed --seed alone was not enough for run-to-run reproducibility,
       since Vina's multi-threaded search still drifted by ~0.1 A between
       identical runs — single-threaded removes that.

       Returns None when this profile has no reference ligand at all
       (site_source == "none_validated" — e.g. Blind-only targets), or a
       redock_reference()-shaped dict otherwise. status
       "no_crystal_reference" means a reference ligand IS known but no
       bond-order-correct crystal pose exists to compare against (most
       commonly: a manually-picked structure whose crystal-SDF build
       failed) — informational, never blocks the compound docking that
       triggered this call."""
    if not profile:
        return None
    resname = profile.get("reference_ligand_resname")
    if not resname:
        return None
    crystal_sdf = profile.get("crystal_sdf")
    reference_smiles = profile.get("reference_smiles")
    if not crystal_sdf or not os.path.exists(crystal_sdf):
        # A manually-picked structure (docking_receptor_custom /
        # docking_alternate_ligand_build) always tries to build its OWN
        # crystal_sdf and, on failure, records why in crystal_sdf_error
        # instead — see those endpoints. Falling through to the registry
        # DEFAULT target's crystal_ligand.sdf below in that case was a real
        # bug: that file is a DIFFERENT PDB entry/ligand (whatever the
        # target's automatic default happens to be), in a completely
        # different receptor's coordinate frame, so "redocking" against it
        # produced either a meaningless/erroring comparison or, worse, an
        # apparently-successful RMSD for the WRONG molecule entirely — this
        # is what made "redock a different/custom PDB" look broken (while
        # the automatic-default path, which legitimately has no crystal_sdf
        # key at all, correctly used the fallback). A profile that never
        # attempted a custom build (the plain registry default) has neither
        # key set, so it still falls through to the fallback below exactly
        # as before.
        if "crystal_sdf_error" in profile:
            return {"status": "no_crystal_reference", "reference_ligand_resname": resname,
                    "crystal_sdf_error": profile["crystal_sdf_error"]}
        from . import profile as DOCK_PROFILE
        target_id = profile.get("target_id")
        # The registry's own default structure keeps its crystal pose at
        # this fixed, conventional path (see scripts/validate_target.py)
        # rather than in the profile dict itself.
        cand = os.path.join(DOCK_PROFILE.DOCKING_TARGETS_DIR, target_id, "crystal_ligand.sdf") if target_id else None
        crystal_sdf = cand if cand and os.path.exists(cand) else None
    if not crystal_sdf:
        return {"status": "no_crystal_reference", "reference_ligand_resname": resname}
    if not reference_smiles:
        crystal_mol = next((m for m in Chem.SDMolSupplier(crystal_sdf, removeHs=True) if m), None)
        if crystal_mol is None:
            return {"status": "no_crystal_reference", "reference_ligand_resname": resname}
        reference_smiles = Chem.MolToSmiles(crystal_mol)
    engine = VinaEngine(exhaustiveness=exhaustiveness, cpu=1)
    result = redock_reference(profile, reference_smiles, crystal_sdf=crystal_sdf, engine=engine, n_poses=n_poses)
    result["reference_ligand_resname"] = resname
    return result