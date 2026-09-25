import type { RedockingValidationResult } from "../lib/types";
import { PoseOverlayViewer } from "./PoseOverlayViewer";

const STATUS_MESSAGES: Record<string, string> = {
  no_crystal_reference: "No crystal reference pose available for this structure — redocking validation skipped.",
  engine_unavailable: "Docking engine unavailable.",
  ligand_prep_failed: "Could not prepare the reference ligand for docking.",
  no_pose: "Vina returned no valid pose for the reference ligand.",
};

/** Redocking-pose validation result, shown alongside a docking/screen
    job's compound results — it ran automatically as part of THIS
    submission (same PDB/exhaustiveness/num poses used for the compounds
    above), not as a separate action. Per-structure, not per-compound: one
    number for the whole run, regardless of how many compounds were
    submitted. See docking/pipeline.py's redock_reference_for_profile. */
export function RedockingValidationNote({ result }: { result: RedockingValidationResult | null | undefined }) {
  if (!result) return null;
  const hasRmsd = result.status === "ok" && result.reference_rmsd != null;
  return (
    <div className="mb-3 rounded-lg border border-line bg-surface px-3 py-2 text-[12.5px]">
      <div className="mb-0.5 font-semibold text-ink">Redocking validation</div>
      {hasRmsd ? (
        <>
          <div className="flex flex-wrap items-center gap-x-4 gap-y-0.5">
            <span className="text-inkmut">
              Re-docked reference ligand{result.reference_ligand_resname ? <> <b className="text-ink">{result.reference_ligand_resname}</b></> : ""} vs. its crystal pose:
            </span>
            <span>
              RMSD <b>{result.reference_rmsd!.toFixed(3)} Å</b>
            </span>
            <span className={result.validated ? "text-brand-800" : "text-clay"}>
              {result.validated ? "validated (< 2 Å)" : "not validated (≥ 2 Å)"}
            </span>
          </div>
          <PoseOverlayViewer crystalPoseSdf={result.crystal_pose_sdf} redockedPosePdb={result.redocked_pose_pdb} />
        </>
      ) : (
        <div className="text-inkmut">
          {result.status === "no_crystal_reference" && result.crystal_sdf_error
            ? `Could not build a crystal reference pose for this structure's ${result.reference_ligand_resname || "reference"} ligand — ${result.crystal_sdf_error}`
            : STATUS_MESSAGES[result.status] || result.error || "Redocking validation could not run."}
        </div>
      )}
      <div className="mt-0.5 text-inkmut">
        Sanity check on this structure's own docking setup (same PDB/exhaustiveness/num poses as the run above) — not a score for any
        compound you submitted.
      </div>
    </div>
  );
}
