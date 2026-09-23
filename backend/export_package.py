"""
B13's "one-click experiment package" — bundles EVERYTHING about one
completed docking/screen run into a single ZIP: metadata (A6's
reproducibility snapshot), the receptor, per-compound complex poses,
interaction diagrams + a consolidated interactions table, GNINA CNN
rescoring, Fresh Decoy Validation (when run for a compound), and the
target-level redocking validation. Shared by both /api/docking/job/{jid}/
export_package and /api/screen/job/{jid}/export_package since the
per-compound result shape (DockResultRow) is the same either way.
"""
import base64
import csv
import io
import json
import os
import re
import zipfile


def _safe_name(smiles, idx):
    s = re.sub(r"[^A-Za-z0-9]+", "_", smiles)[:40] or "compound"
    return f"{idx:03d}_{s}"


def _combine_pdb_text(receptor_pdb, pose_pdb):
    """Same strip-terminators-then-TER-then-END logic as the frontend's
       mol3d.ts combinePdbText — duplicated here (not imported, there's no
       shared JS/Python boundary) because a server-side ZIP export can't
       depend on client-side JS having run."""
    def strip(s):
        lines = [ln for ln in s.split("\n") if not re.match(r"^(END|ENDMDL)\s*$", ln.strip())]
        return "\n".join(lines).rstrip("\n")
    return f"{strip(receptor_pdb)}\nTER\n{strip(pose_pdb)}\nEND\n"


def build_zip(run_metadata: dict, receptor_pdb_path, results: list,
              fresh_decoy_results: dict | None = None, redocking_validation: dict | None = None) -> bytes:
    """fresh_decoy_results: {smiles: fresh_decoy_validation()'s result dict}
       — only compounds a user actually clicked "Run Fresh Decoy Validation"
       for will have an entry; not run for the rest, same as in the UI.
       redocking_validation: the target-level (not per-compound) redocking_
       reference_for_profile() result, if this run's structure had one."""
    buf = io.BytesIO()
    receptor_pdb_text = None
    if receptor_pdb_path and os.path.exists(receptor_pdb_path):
        with open(receptor_pdb_path) as f:
            receptor_pdb_text = f.read()
    fresh_decoy_results = fresh_decoy_results or {}

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("metadata.json", json.dumps(run_metadata, indent=2, default=str))

        if receptor_pdb_text:
            z.writestr("receptor.pdb", receptor_pdb_text)

        csv_buf = io.StringIO()
        w = csv.writer(csv_buf)
        w.writerow(["smiles", "status", "vina_score", "confidence", "n_valid", "category", "reason", "suggested_action",
                   "gnina_cnn_score", "gnina_cnn_affinity", "gnina_affinity",
                   "fresh_decoy_percentile", "fresh_decoy_discrimination", "fresh_decoy_compound_score"])
        for r in results:
            g = r.get("gnina") or {}
            fd = fresh_decoy_results.get(r.get("smiles")) or {}
            w.writerow([r.get("smiles"), r.get("status"), r.get("vina_score"), r.get("confidence"),
                       r.get("n_valid"), r.get("category"), r.get("reason"), r.get("suggested_action"),
                       g.get("cnn_score"), g.get("cnn_affinity"), g.get("gnina_affinity"),
                       fd.get("percentile"), fd.get("discrimination"), fd.get("compound_score")])
        z.writestr("results.csv", csv_buf.getvalue())

        # Consolidated interaction log — one row per interaction EVENT
        # across every compound (same shape as the UI's own "Interaction
        # table"), not just the per-compound PNG diagram.
        inter_buf = io.StringIO()
        iw = csv.writer(inter_buf)
        iw.writerow(["compound", "interacting_residue", "bond_type", "length_angstrom"])
        for r in results:
            for h in (r.get("interactions") or []):
                residue = h.get("residue") or (f"{h['chain']}:{h['resid']}" if h.get("chain") and h.get("resid") is not None else "")
                iw.writerow([r.get("smiles"), residue, h.get("label") or h.get("category") or h.get("type"), h.get("distance")])
        z.writestr("interactions.csv", inter_buf.getvalue())

        for i, r in enumerate(results):
            name = _safe_name(r.get("smiles") or "", i)
            pose_pdb = r.get("pose_pdb")
            if pose_pdb:
                text = _combine_pdb_text(receptor_pdb_text, pose_pdb) if receptor_pdb_text else pose_pdb
                z.writestr(f"poses/{name}.pdb", text)
            png_b64 = r.get("interaction_png")
            if png_b64:
                try:
                    z.writestr(f"interactions/{name}.png", base64.b64decode(png_b64))
                except Exception:
                    pass

            fd = fresh_decoy_results.get(r.get("smiles"))
            if fd:
                z.writestr(f"fresh_decoy/{name}.json", json.dumps(fd, indent=2, default=str))
                decoys = fd.get("decoys") or []
                if decoys:
                    fd_buf = io.StringIO()
                    fw = csv.writer(fd_buf)
                    fw.writerow(["name", "smiles", "source_target", "score"])
                    for d in decoys:
                        fw.writerow([d.get("name"), d.get("smiles"), d.get("source_target"), d.get("score")])
                    z.writestr(f"fresh_decoy/{name}_decoys.csv", fd_buf.getvalue())

        if redocking_validation:
            rv = dict(redocking_validation)
            redocked_pose_pdb = rv.pop("redocked_pose_pdb", None)
            crystal_pose_sdf = rv.pop("crystal_pose_sdf", None)
            z.writestr("redocking_validation/summary.json", json.dumps(rv, indent=2, default=str))
            if redocked_pose_pdb:
                z.writestr("redocking_validation/redocked_pose.pdb", redocked_pose_pdb)
            if crystal_pose_sdf:
                z.writestr("redocking_validation/experimental_pose.sdf", crystal_pose_sdf)
    return buf.getvalue()
