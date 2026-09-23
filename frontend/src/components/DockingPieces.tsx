import { useState } from "react";
import * as api from "../lib/api";
import { apiUrl } from "../lib/api";
import { combinePdbText, fetchTextCached } from "../lib/mol3d";
import type { AdvancedDockingBody, AlternateLigand, DockResultRow } from "../lib/types";
import { PoseViewer } from "./PoseViewer";
import { SwapIcon, DownloadIcon, ZoomIcon, CubeIcon } from "./Icons";
import { Modal } from "./Modal";

// Matches backend/docking/failure_diagnostics.py's category slugs exactly.
const FAILURE_CATEGORY_LABELS: Record<string, string> = {
  invalid_molecule: "Invalid molecule",
  conformer_generation_failed: "3D conformer generation failed",
  ligand_conversion_failed: "Ligand-to-PDBQT conversion failed",
  ligand_prep_failed: "Ligand preparation failed",
  engine_unavailable: "Docking engine unavailable",
  engine_error: "Docking engine error",
  pose_validity_failed: "No physically valid pose",
  no_pose: "No pose produced",
  other: "Failed",
};

/** Downloads the receptor+pose "complex" PDB for one docked compound —
    the same two structures PoseViewer already renders together in 3D,
    just written out as a real file. receptorPdbPath is the file that was
    ACTUALLY docked against for this job (from the job's own response),
    not necessarily whatever the target's current default happens to be
    now. */
export function DownloadComplexButton({
  smiles,
  posePdb,
  receptorPdbPath,
}: {
  smiles: string;
  posePdb?: string | null;
  receptorPdbPath?: string | null;
}) {
  const [busy, setBusy] = useState(false);
  if (!posePdb) return null;

  const run = async () => {
    setBusy(true);
    try {
      const receptorPdb = receptorPdbPath
        ? await fetchTextCached(apiUrl(`/api/docking/receptor_file?path=${encodeURIComponent(receptorPdbPath)}`))
        : null;
      const text = receptorPdb ? combinePdbText(receptorPdb, posePdb) : posePdb;
      const blob = new Blob([text], { type: "chemical/x-pdb" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      const safeName = smiles.replace(/[^A-Za-z0-9]+/g, "_").slice(0, 40) || "compound";
      a.href = url;
      a.download = `${safeName}${receptorPdb ? "_complex" : "_pose"}.pdb`;
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
    } finally {
      setBusy(false);
    }
  };

  return (
    <button type="button" className="btn-link" disabled={busy} onClick={run}>
      <DownloadIcon className="h-3.5 w-3.5" />
      {busy ? "Preparing…" : receptorPdbPath ? "Download complex (PDB)" : "Download pose (PDB)"}
    </button>
  );
}

export function InteractionTable({ interactions }: { interactions?: DockResultRow["interactions"] }) {
  if (!interactions || !interactions.length) return null;
  const rows = [...interactions].sort((a, b) => (a.distance ?? 0) - (b.distance ?? 0));
  return (
    <div className="card p-3.5">
      <div className="mb-2.5 text-[10.5px] font-bold uppercase tracking-wide text-brand-700">
        Protein-ligand nonbonding interactions ({rows.length})
      </div>
      <div className="max-h-[300px] overflow-y-auto rounded-lg border border-line">
        <table className="w-full border-collapse text-[12.5px]">
          <thead>
            <tr>
              {["Name", "Category", "Type", "Distance (Å)"].map((h) => (
                <th key={h} className="sticky top-0 z-10 border-b border-line bg-surface2 px-2.5 py-1.5 text-left text-[10px] font-semibold uppercase tracking-wide text-inkmut">
                  {h}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((r, i) => (
              <tr key={i} className="hover:bg-canvas">
                <td className="border-b border-surface2 px-2.5 py-1.5">{r.name || r.residue}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5">{r.category || ""}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5">{r.label || r.type || ""}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5">{r.distance ?? "—"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export function FreshDecoyButton({
  smiles,
  targetId,
  advanced,
  parentKind,
  parentJobId,
}: {
  smiles: string;
  targetId: string | null;
  advanced: AdvancedDockingBody | null;
  /** The docking/screen job this compound's row belongs to — when given,
      the finished result is ALSO stored on that job server-side so
      "Download full experiment package" includes it, not just this
      inline summary text. */
  parentKind?: "docking" | "screen";
  parentJobId?: string | null;
}) {
  const [busy, setBusy] = useState(false);
  const [label, setLabel] = useState("Run Fresh Decoy Validation");
  const [status, setStatus] = useState<{ kind: "muted" | "ok" | "err"; text: string } | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);

  const run = async () => {
    if (!targetId) {
      setStatus({ kind: "err", text: "No target context for this result." });
      return;
    }
    setBusy(true);
    setLabel("Generating decoys & docking (~2-5 min)…");
    setStatus({ kind: "muted", text: "Generating ~50 decoys matched to this compound…" });
    try {
      const sub = await api.submitFreshDecoy(
        targetId,
        smiles,
        advanced,
        parentKind && parentJobId ? { kind: parentKind, jobId: parentJobId } : null
      );
      setJobId(sub.job_id);
      let j;
      while (true) {
        await api.sleep(2000);
        j = await api.pollRetry(() => api.freshDecoyJob(sub.job_id));
        if (j.status === "done" || j.status === "error" || j.status === "cancelled") break;
        const done = j.done || 0;
        const total = j.total || "?";
        setStatus({ kind: "muted", text: `Docking ${done}/${total}…` });
        setLabel(`Docking ${done}/${total}…`);
      }
      if (j.status === "cancelled") {
        setStatus({ kind: "muted", text: "Stopped by user." });
        return;
      }
      if (j.status === "error") {
        setStatus({ kind: "err", text: j.error || "failed" });
        return;
      }
      const res = j.result!;
      if (res.error) {
        setStatus({ kind: "muted", text: res.error });
        return;
      }
      const ds = res.decoy_stats;
      const rs = res.run_settings;
      setStatus({
        kind: "ok",
        text:
          `Fresh percentile: ${res.percentile}% · Decoy discrimination: ${res.discrimination} — compound score ${res.compound_score} kcal/mol vs ${res.n_decoys_docked} freshly-docked, property-matched & topologically-dissimilar decoys${res.n_decoys_failed ? ` (${res.n_decoys_failed} failed to dock)` : ""}.` +
          (ds ? ` Decoy scores: mean ${ds.mean}, median ${ds.median}, SD ${ds.sd} (range ${ds.min} to ${ds.max}).` : "") +
          (rs ? ` Run: ${rs.docking_mode === "blind" ? "blind" : "site-specific"}${rs.pdb_source ? `, ${rs.pdb_source}` : ""}${rs.exhaustiveness != null ? `, exhaustiveness ${rs.exhaustiveness}` : ""}.` : ""),
      });
    } catch (e: any) {
      setStatus({ kind: "err", text: e.message || "Error" });
    } finally {
      setBusy(false);
      setJobId(null);
      setLabel("Run Fresh Decoy Validation");
    }
  };

  const stop = async () => {
    if (!jobId) return;
    try {
      await api.cancelFreshDecoy(jobId);
    } catch {
      /* the poll loop will still surface a final status either way */
    }
  };

  return (
    <div onClick={(e) => e.stopPropagation()}>
      <button type="button" className="btn-link" disabled={busy} onClick={run}>
        {label}
      </button>
      {busy && jobId && (
        <button type="button" className="btn-link ml-1.5" onClick={stop}>
          Stop
        </button>
      )}
      {status && (
        <div
          className={`mt-1.5 max-w-[220px] text-[12px] ${
            status.kind === "ok" ? "font-medium text-brand-700" : status.kind === "err" ? "text-clay" : "text-inkmut"
          }`}
        >
          {status.text}
        </div>
      )}
    </div>
  );
}

/** One card-style group in the consolidated per-compound detail view.
    Callers decide whether a card applies at all (conditionally rendering
    the whole <DetailCard> rather than this component guessing from its
    children) — a generic "are my children empty" heuristic can't tell a
    real wrapper <div> with nothing conditionally rendered inside it from
    one with real content, so that decision belongs with the caller, who
    actually knows. */
function DetailCard({
  title,
  icon,
  actions,
  children,
  className = "",
}: {
  title: string;
  icon?: React.ReactNode;
  actions?: React.ReactNode;
  children: React.ReactNode;
  className?: string;
}) {
  return (
    <div className={`card flex flex-col p-3.5 ${className}`}>
      <div className="mb-2.5 flex shrink-0 items-center justify-between gap-2">
        <h5 className="flex items-center gap-1.5 text-[10.5px] font-bold uppercase tracking-wide text-brand-700">
          {icon}
          {title}
        </h5>
        {actions}
      </div>
      {children}
    </div>
  );
}

/** A single labeled number/fact in a compact grid — the "Summary" card's
    building block, replacing a wall of inline "label: value ·" text with
    scannable tiles. */
function StatTile({ label, value, tone = "ink" }: { label: string; value: React.ReactNode; tone?: "ink" | "amber" }) {
  return (
    <div className="rounded-lg border border-line bg-surface2/50 px-2.5 py-2">
      <div className="text-[10px] font-semibold uppercase tracking-wide text-inkmut">{label}</div>
      <div className={`mt-0.5 text-[13.5px] font-semibold leading-snug ${tone === "amber" ? "text-amber" : "text-ink"}`}>{value}</div>
    </div>
  );
}

export function DockDetailPanel({
  r,
  receptorPdbPath,
  jobId,
  reportKind = "docking",
}: {
  r: DockResultRow;
  receptorPdbPath?: string | null;
  /** Enables the A5 "Generate research report" button when the caller
      has a completed job id to report against — omitted (e.g. from
      TargetInfoTab's ad-hoc panels, which have no job) hides it. */
  jobId?: string | null;
  reportKind?: "docking" | "screen";
}) {
  const canView = !!r.interaction_png;
  const g = r.gnina;
  const hasDocking = r.vina_score != null || !!r.confidence || (g && (g.cnn_score != null || g.cnn_affinity != null || g.gnina_affinity != null));
  const hasCompInfo = r.n_valid != null || r.pose_self_consistency != null || (!!r.status && r.status !== "ok");
  const hasSummary = hasDocking || hasCompInfo;
  const [lightboxOpen, setLightboxOpen] = useState(false);
  const downloadActions = canView && jobId && (
    <div className="flex gap-2">
      <a
        className="btn-link px-2 py-1 text-[11px]"
        href={api.apiUrl(`/api/${reportKind}/job/${jobId}/interaction_diagram?smiles=${encodeURIComponent(r.smiles)}&fmt=svg`)}
        download
      >
        <DownloadIcon className="h-3 w-3" />
        SVG
      </a>
      <a
        className="btn-link px-2 py-1 text-[11px]"
        href={api.apiUrl(`/api/${reportKind}/job/${jobId}/interaction_diagram?smiles=${encodeURIComponent(r.smiles)}&fmt=tiff`)}
        download
      >
        <DownloadIcon className="h-3 w-3" />
        TIFF
      </a>
    </div>
  );
  return (
    <div className="bg-canvas/60 px-5 py-4">
      {r.suggested_action && (
        <div className="mb-4 rounded-lg border border-amber/30 bg-amber/10 px-3.5 py-2.5 text-[12.5px]">
          <div className="text-ink">
            <b>{FAILURE_CATEGORY_LABELS[r.category || ""] || r.category || "Failed"}</b>
            {r.reason || r.error ? ` — ${r.reason || r.error}` : ""}
          </div>
          <div className="mt-1 text-amber">{r.suggested_action}</div>
        </div>
      )}
      <div className="grid grid-cols-1 gap-4 lg:grid-cols-[1.15fr_1fr]">
        <div className="flex flex-col">
          <DetailCard title="Binding site & interactions" actions={downloadActions} className="flex-1">
            {(r.residue_overlap_pct != null || (canView && r.interaction_source)) && (
              <div className="mb-2 flex shrink-0 flex-wrap gap-x-4 gap-y-0.5 text-[12px] text-inkmut">
                {r.residue_overlap_pct != null && <span>Shares {r.residue_overlap_pct}% of the reference drug's contact residues.</span>}
                {canView && r.interaction_source && <span>Interaction detection: {r.interaction_source}.</span>}
              </div>
            )}
            {canView ? (
              <button
                type="button"
                onClick={() => setLightboxOpen(true)}
                className="group relative flex min-h-[260px] flex-1 cursor-zoom-in items-center justify-center overflow-hidden rounded-lg border border-line bg-white"
              >
                <img src={`data:image/png;base64,${r.interaction_png}`} className="h-full w-full object-contain" />
                <span className="pointer-events-none absolute inset-0 flex items-center justify-center bg-ink/0 opacity-0 transition group-hover:bg-ink/5 group-hover:opacity-100">
                  <span className="flex items-center gap-1.5 rounded-lg bg-surface/95 px-3 py-1.5 text-[12px] font-semibold text-ink shadow-card">
                    <ZoomIcon className="h-3.5 w-3.5" />
                    Enlarge
                  </span>
                </span>
              </button>
            ) : (
              <div className="py-2 text-[13px] text-inkmut">No interaction diagram for this pose.</div>
            )}
          </DetailCard>
        </div>

        <div className="space-y-4">
          {hasSummary && (
            <DetailCard title="Summary">
              <div className="grid grid-cols-2 gap-2">
                {r.vina_score != null && <StatTile label="Vina score" value={`${r.vina_score} kcal/mol`} />}
                {r.confidence && <StatTile label="Confidence" value={r.confidence} />}
                {g?.cnn_score != null && <StatTile label="GNINA CNN score" value={g.cnn_score} />}
                {g?.cnn_affinity != null && <StatTile label="GNINA CNN affinity" value={g.cnn_affinity} />}
                {g?.gnina_affinity != null && <StatTile label="GNINA affinity" value={`${g.gnina_affinity} kcal/mol`} />}
                {r.n_valid != null && <StatTile label="PoseBusters-valid poses" value={r.n_valid} />}
                {r.pose_self_consistency != null && <StatTile label="Pose self-consistency" value={r.pose_self_consistency} />}
                {r.status && r.status !== "ok" && (
                  <StatTile label="Status" value={r.reason ? `${r.status} — ${r.reason}` : r.status} tone="amber" />
                )}
              </div>
            </DetailCard>
          )}

          {r.pose_pdb && (
            <DetailCard title="3D pose" icon={<CubeIcon className="h-3.5 w-3.5" />}>
              <div className="flex flex-wrap items-center gap-2">
                <PoseViewer posePdb={r.pose_pdb} receptorPdbPath={receptorPdbPath} interactions={r.interactions} />
                <DownloadComplexButton smiles={r.smiles} posePdb={r.pose_pdb} receptorPdbPath={receptorPdbPath} />
              </div>
            </DetailCard>
          )}

          <InteractionTable interactions={r.interactions} />
        </div>
      </div>
      {lightboxOpen && canView && (
        <Modal title="Ligand interaction diagram" sub={r.smiles} onClose={() => setLightboxOpen(false)}>
          <div className="flex flex-col items-center gap-4 p-5">
            <img src={`data:image/png;base64,${r.interaction_png}`} className="max-w-full rounded-lg border border-line bg-white" />
            {downloadActions}
          </div>
        </Modal>
      )}
    </div>
  );
}

/** A5 — assembles and shows the full evidence-chain report for one
    compound (natural source, chemical identity, literature, target
    prediction, QSAR, docking, interactions, ADMET, off-target,
    summary) on demand, since it's a slower call (a live PubMed request
    is part of it) that most users won't want for every row. */
export function ResearchReportButton({ jobId, smiles, kind }: { jobId: string; smiles: string; kind: "docking" | "screen" }) {
  const [state, setState] = useState<"idle" | "loading" | "error" | "done">("idle");
  const [error, setError] = useState("");
  const [data, setData] = useState<{ report: any; markdown: string } | null>(null);

  const generate = async () => {
    setState("loading");
    setError("");
    try {
      const r = await api.researchReport(kind, jobId, smiles, true);
      setData(r);
      setState("done");
    } catch (e: any) {
      setError(e.message || "Error");
      setState("error");
    }
  };

  const downloadMarkdown = () => {
    if (!data) return;
    const blob = new Blob([data.markdown], { type: "text/markdown" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `research_report_${(smiles || "compound").slice(0, 24).replace(/[^A-Za-z0-9]+/g, "_")}.md`;
    a.click();
    URL.revokeObjectURL(url);
  };

  if (state === "idle" || state === "error") {
    return (
      <div>
        <button type="button" className="btn-link" onClick={generate}>
          Generate research report
        </button>
        {state === "error" && <div className="field-hint text-clay">{error}</div>}
      </div>
    );
  }
  if (state === "loading") {
    return <div className="text-[12.5px] text-inkmut">Assembling evidence chain (includes a live PubMed lookup)…</div>;
  }
  const r = data!.report;
  return (
    <div>
      <div className="mb-2.5 flex items-center justify-between">
        <div className="text-[11px] text-inkmut">{r.pipeline_stages.join(" → ")}</div>
        <button type="button" className="btn-link shrink-0" onClick={downloadMarkdown}>
          Download report (.md)
        </button>
      </div>
      <div className="rounded-lg border border-line bg-surface1 p-3 text-[12.5px] leading-relaxed text-ink">{r.evidence_summary}</div>
      <ReportField label="Natural source" value={r.natural_source?.plant_source} />
      <ReportField
        label="Chemical identity"
        value={r.chemical_identity?.available && `${r.chemical_identity.molecular_formula}, MW ${r.chemical_identity.molecular_weight}, LogP ${r.chemical_identity.logp}`}
      />
      <ReportField
        label="Reported activity"
        value={
          r.reported_activity?.available
            ? `${r.reported_activity.n_results} paper(s) for "${r.reported_activity.query}"`
            : r.reported_activity?.note
        }
      />
      <ReportField
        label="Target prediction"
        value={r.target_prediction?.available && (r.target_prediction.on_target_supported ? "Supported by similar known actives" : "No similar known actives found")}
      />
      <ReportField
        label="Off-target analysis"
        value={r.off_target_analysis?.available && `${r.off_target_analysis.n_off_targets ?? 0} other target(s) with similarity signal`}
      />
      <details className="mt-2">
        <summary className="cursor-pointer text-[11.5px] font-semibold text-brand-700">Methods (draft)</summary>
        <p className="mt-1 text-[12px] leading-relaxed text-inkmut">{r.methods_draft}</p>
      </details>
    </div>
  );
}

function ReportField({ label, value }: { label: string; value?: string | null | false }) {
  if (!value) return null;
  return (
    <div className="mt-1.5 text-[12px] text-inkmut">
      <b className="text-ink">{label}:</b> {value}
    </div>
  );
}

/** "Dock again with a different ligand": the same raw PDB structure a job
    docked against often has MORE than one real co-crystallized ligand
    (a second binding site, or one copy per chain in a crystallographic
    dimer) — the pipeline always silently centers the box on just the
    single largest one. This lets the user pick a different real ligand
    from the same structure and redock, with every other setting
    (exhaustiveness, poses, GNINA, compound list) held identical to the
    original run — the actual "build a new receptor + resubmit" work is
    owned by the caller (onRedock), since that has to replace the whole
    results table the same way "Reproduce this analysis" does; this
    component only owns fetching the candidate list and letting the user
    pick one. */
export function AlternateLigandButton({
  jobId,
  kind,
  onRedock,
  busy,
}: {
  jobId: string;
  kind: "docking" | "screen";
  onRedock: (lig: AlternateLigand) => void;
  busy?: boolean;
}) {
  const [state, setState] = useState<"idle" | "loading" | "picking" | "none" | "error">("idle");
  const [error, setError] = useState("");
  const [data, setData] = useState<{ current: AlternateLigand | null; ligands: AlternateLigand[] } | null>(null);
  const [picked, setPicked] = useState("");

  const open = async () => {
    setState("loading");
    setError("");
    try {
      const d = kind === "docking" ? await api.dockingAlternateLigands(jobId) : await api.screenAlternateLigands(jobId);
      if (!d.available || d.ligands.length < 2) {
        setState("none");
        return;
      }
      setData({ current: d.current, ligands: d.ligands });
      setPicked("");
      setState("picking");
    } catch (e: any) {
      setError(e.message || "Error");
      setState("error");
    }
  };

  const confirm = () => {
    if (!data || !picked) return;
    const lig = data.ligands.find((l) => `${l.chain}:${l.resnum}` === picked);
    if (lig) onRedock(lig);
    setState("idle");
  };

  if (state === "idle" || state === "error") {
    return (
      <span>
        <button type="button" className="btn-link" onClick={open} disabled={busy}>
          <SwapIcon className="h-3.5 w-3.5" />
          Dock again with a different ligand
        </button>
        {state === "error" && <div className="field-hint text-clay">{error}</div>}
      </span>
    );
  }
  if (state === "loading") return <span className="text-[12.5px] text-inkmut">Checking for other ligands in this structure…</span>;
  if (state === "none") return <span className="text-[12.5px] text-inkmut">Only one real ligand found in this structure.</span>;

  return (
    <span className="inline-flex flex-wrap items-center gap-2">
      <select className="field-input inline-block w-auto text-[12.5px]" value={picked} onChange={(e) => setPicked(e.target.value)}>
        <option value="" disabled>
          Pick a different ligand…
        </option>
        {data!.ligands.map((l) => {
          const key = `${l.chain}:${l.resnum}`;
          const isCurrent = !!data!.current && l.chain === data!.current.chain && l.resnum === data!.current.resnum;
          return (
            <option key={key} value={key} disabled={isCurrent}>
              {l.resname} · chain {l.chain} · residue {l.resnum}
              {isCurrent ? " (current)" : ""}
            </option>
          );
        })}
      </select>
      <button type="button" className="btn-link" onClick={confirm} disabled={!picked || busy}>
        {busy ? "Docking…" : "Dock"}
      </button>
      <button type="button" className="btn-link" onClick={() => setState("idle")} disabled={busy}>
        Cancel
      </button>
    </span>
  );
}
