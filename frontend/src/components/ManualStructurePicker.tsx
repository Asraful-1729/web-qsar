import { useEffect, useState } from "react";
import type { AdvancedDockingState } from "../lib/useAdvancedDocking";
import type { ReceptorProfile } from "../lib/types";
import { ReceptorBeforeAfter } from "./ReceptorBeforeAfter";
import { Modal } from "./Modal";
import { SplitIcon, DownloadIcon } from "./Icons";

/** B13 — client-side download of the already-fetched prep_report data
    (same pattern as the enrichment plots' "Download PNG": no new
    endpoint needed, the facts are already in the profile response). */
function downloadPrepReport(profile: ReceptorProfile) {
  const lines = [
    `Receptor preparation report`,
    `Target: ${profile.target_id ?? "?"}`,
    `Structure: ${profile.pdb_source ?? "?"}`,
    "",
    ...(profile.prep_report || []).map((s, i) => `${i + 1}. ${s.label}\n   ${s.detail}`),
  ].filter((l): l is string => l != null);
  const blob = new Blob([lines.join("\n")], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `prep_report_${(profile.target_id || "receptor").replace(/[^A-Za-z0-9]+/g, "_")}.txt`;
  a.click();
  URL.revokeObjectURL(url);
}

/** Real per-phase wall time for this structure-prep run — makes a slow
    run (network fetch vs. the offline build itself) diagnosable at a
    glance instead of one opaque timed-but-unbroken-down wait. */
function formatTiming(t: NonNullable<ReceptorProfile["timing"]>): string {
  const parts: string[] = [];
  if (t.fetch_pdb_seconds != null) parts.push(`fetch structure ${t.fetch_pdb_seconds}s`);
  if (t.build_seconds != null) parts.push(`strip/repair/PDBQT ${t.build_seconds}s`);
  const total = Object.values(t).reduce((a, b) => a + (b || 0), 0);
  return `${parts.join(" · ")} (total ${total.toFixed(1)}s)`;
}

const STATUS_CLS: Record<string, string> = {
  muted: "text-inkmut",
  ok: "border-brand-300/50 bg-brand-500/[0.08] text-brand-800 rounded-lg border px-2.5 py-2",
  warn: "border-amber/30 bg-amber/10 text-amber rounded-lg border px-2.5 py-2",
  err: "text-clay",
};

/** Manual structure override — lives right below the Target field
    (TargetBrowser.tsx), not tucked inside the collapsed Advanced Settings
    panel, since picking a different structure is a common enough action
    to want visible without opening that panel first. */
export function ManualStructurePicker({ adv }: { adv: AdvancedDockingState }) {
  const [beforeAfterOpen, setBeforeAfterOpen] = useState(false);

  // Close on a target switch — it would otherwise keep showing the
  // PREVIOUS target's before/after comparison while the new one loads.
  useEffect(() => {
    setBeforeAfterOpen(false);
  }, [adv.targetId]);

  return (
    <div className="mt-3">
      <label className="field-label">Manual structure (overrides the automatic recommendation)</label>
      <div className="max-h-[200px] overflow-y-auto rounded-lg border border-line">
        {adv.candidatesLoading && <div className="p-2 text-[12.5px] text-inkmut">Loading structural evidence…</div>}
        {!adv.candidatesLoading && adv.candidates?.length === 0 && (
          <div className="p-2 text-[12.5px] text-inkmut">No qualifying structures on record.</div>
        )}
        {!adv.candidatesLoading &&
          adv.candidates?.map((c) => {
            const q = [
              c.resolution != null ? `${c.resolution} Å` : null,
              c.ligand_RSCC != null ? `RSCC ${c.ligand_RSCC}` : null,
              c.ligand_RSR != null ? `RSR ${c.ligand_RSR}` : null,
            ]
              .filter(Boolean)
              .join(" · ");
            // No explicit manual pick yet -> visually treat the
            // registry's current automatic default as the selected
            // one (it already IS what a submit would use), so the
            // user sees what's actually going to run without having
            // to click anything first.
            const on = adv.pickedPdb ? adv.pickedPdb === c.pdb_id : c.is_current_default;
            return (
              <div
                key={c.pdb_id}
                // A blind_only candidate has no real ligand to choose
                // from at all — build it directly. Every other candidate
                // just SELECTS the PDB (lists its co-crystallized ligands
                // below); nothing gets built until the user explicitly
                // picks one of them.
                onClick={() => (c.blind_only ? adv.pickLigand(c.pdb_id, null, null) : adv.selectPdb(c.pdb_id))}
                className={`flex cursor-pointer items-center gap-2 border-b border-line/70 px-2.5 py-1.5 text-[12.5px] last:border-0 hover:bg-surface2/60 ${on ? "bg-brand-500/[0.08]" : ""}`}
              >
                <span className="min-w-[48px] font-bold text-ink">{c.pdb_id}</span>
                <span className="flex-1 text-inkmut">
                  rank #{c.csv_rank ?? "?"} · {c.blind_only ? "no validated ligand — Blind mode only" : `suggested: ${c.resname}`} · {q}
                </span>
                {c.is_current_default && (
                  <span className="flex items-center gap-1">
                    <span className="badge bg-brand-500/15 text-brand-800">automatic default</span>
                  </span>
                )}
                {c.blind_only && (
                  <span className="flex items-center gap-1">
                    <span className="badge bg-amber/15 text-amber">blind only</span>
                  </span>
                )}
              </div>
            );
          })}
      </div>
      {adv.pickedPdb && (adv.ligandOptionsLoading || adv.ligandOptions || adv.ligandOptionsError) && (
        <div className="mt-2">
          <label className="field-label">Co-crystallized ligands — {adv.pickedPdb}</label>
          <p className="field-hint mb-1">
            Every real ligand found in this structure. Pick one to build the receptor and binding box around — which one is biologically
            relevant depends on this structure, so none is pre-selected for you.
          </p>
          {adv.ligandOptionsLoading && <div className="p-2 text-[12.5px] text-inkmut">Fetching {adv.pickedPdb} and listing its ligands…</div>}
          {adv.ligandOptionsError && <div className="p-2 text-[12.5px] text-clay">{adv.ligandOptionsError}</div>}
          {!adv.ligandOptionsLoading && adv.ligandOptions && adv.ligandOptions.length === 0 && (
            <div className="p-2 text-[12.5px] text-inkmut">No real co-crystallized ligand found in {adv.pickedPdb} — use Blind mode.</div>
          )}
          {!adv.ligandOptionsLoading && !!adv.ligandOptions?.length && (
            <div className="max-h-[160px] overflow-y-auto rounded-lg border border-line">
              {adv.ligandOptions.map((lig) => {
                const key = `${lig.chain}:${lig.resnum}:${lig.resname}`;
                const on = adv.customProfile?.reference_ligand_resname === lig.resname && adv.customProfile?.chain === lig.chain;
                return (
                  <div
                    key={key}
                    onClick={() => adv.pickedPdb && adv.pickLigand(adv.pickedPdb, lig.resname, lig.chain)}
                    className={`flex cursor-pointer items-center gap-2 border-b border-line/70 px-2.5 py-1.5 text-[12.5px] last:border-0 hover:bg-surface2/60 ${on ? "bg-brand-500/[0.08]" : ""}`}
                  >
                    <span className="min-w-[64px] font-bold text-ink">{lig.resname}</span>
                    <span className="flex-1 text-inkmut">
                      chain {lig.chain} · residue {lig.resnum} · {lig.n_atoms} heavy atoms
                    </span>
                    {on && <span className="badge bg-brand-500/15 text-brand-800">selected</span>}
                  </div>
                );
              })}
            </div>
          )}
        </div>
      )}
      {adv.structureStatus && <div className={`field-hint ${STATUS_CLS[adv.structureStatus.kind]}`}>{adv.structureStatus.text}</div>}
      {!!adv.customProfile?.ligand_sanity_warnings?.length && (
        <div className="mt-1.5 rounded-lg border border-amber/30 bg-amber/10 px-2.5 py-2 text-[12px] text-amber">
          <b>Ligand sanity check flagged this pick:</b>
          <ul className="mt-1 list-disc pl-4">
            {adv.customProfile.ligand_sanity_warnings.map((w, i) => (
              <li key={i}>{w}</li>
            ))}
          </ul>
        </div>
      )}
      {adv.customProfile?.raw_pdb_path && adv.customProfile?.receptor_pdb && (
        <button type="button" className="btn-ghost mt-1.5" onClick={() => setBeforeAfterOpen(true)}>
          <SplitIcon className="h-3.5 w-3.5" />
          Compare before/after receptor prep
        </button>
      )}
      {beforeAfterOpen && adv.customProfile?.raw_pdb_path && adv.customProfile?.receptor_pdb && (
        <Modal
          title="Receptor preparation — before/after"
          sub={`${adv.customProfile.target_id ?? ""} — ${adv.customProfile.pdb_source ?? "?"}`}
          onClose={() => setBeforeAfterOpen(false)}
        >
          <div className="h-full p-3.5">
            <ReceptorBeforeAfter
              rawPdbPath={adv.customProfile.raw_pdb_path}
              cleanPdbPath={adv.customProfile.receptor_pdb}
              heightClassName="h-full"
            />
          </div>
        </Modal>
      )}
      {!!adv.customProfile?.prep_report?.length && (
        <button type="button" className="btn-ghost mt-1.5" onClick={() => downloadPrepReport(adv.customProfile!)}>
          <DownloadIcon className="h-3.5 w-3.5" />
          Download prep report
        </button>
      )}
      {!!adv.customProfile?.timing && <div className="field-hint mt-1">Timing: {formatTiming(adv.customProfile.timing)}</div>}
    </div>
  );
}
