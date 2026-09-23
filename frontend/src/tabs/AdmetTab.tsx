import { useState } from "react";
import * as api from "../lib/api";
import { useMoleculeInput } from "../lib/useMoleculeInput";
import { MoleculeInputPanel } from "../components/MoleculeInputPanel";
import { SectionIntro, SidebarLayout, useSidebarCollapsed } from "../components/Shell";
import { Disclaimer, EmptyState, ErrorBox, Notice, ProgressBar, Spinner } from "../components/Feedback";
import type { AdmetProfile, AdmetResponse } from "../lib/types";

function csvCell(v: unknown): string {
  if (v === null || v === undefined) return "";
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

/** Every ADMET-AI column across every profile (not just the curated
    subset the results table shows), one row per submitted compound —
    see admet.py's grouped_learned() "raw" field. Unparsed compounds or
    ones the worker had no prediction for still get a row (smiles only,
    columns blank) rather than being silently dropped. */
function downloadAdmetCsv(d: AdmetResponse) {
  const rawKeys: string[] = [];
  const seen = new Set<string>();
  for (const p of d.profiles) {
    for (const k of Object.keys(p.learned?.raw || {})) {
      if (!seen.has(k)) {
        seen.add(k);
        rawKeys.push(k);
      }
    }
  }
  const header = ["smiles", ...rawKeys];
  const lines = [header.map(csvCell).join(",")];
  for (const p of d.profiles) {
    const smiles = p.standardised_smiles || p.input_smiles;
    const raw = p.learned?.raw || {};
    lines.push([smiles, ...rawKeys.map((k) => raw[k])].map(csvCell).join(","));
  }
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "admet_full.csv";
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

const GROUP_ORDER = ["Absorption", "Distribution", "Metabolism", "Excretion", "Toxicity", "Other"];
const TONE_BORDER: Record<string, string> = {
  good: "border-l-brand-500",
  warn: "border-l-amber",
  bad: "border-l-clay",
  neutral: "border-l-slateout",
};

export function AdmetTab() {
  const mol = useMoleculeInput();
  const [state, setState] = useState<"idle" | "loading" | "progress" | "error" | "done">("idle");
  const [error, setError] = useState("");
  const [pct, setPct] = useState(0);
  const [total, setTotal] = useState(0);
  const [data, setData] = useState<AdmetResponse | null>(null);
  const [sidebarCollapsed, toggleSidebar] = useSidebarCollapsed();

  const run = async () => {
    setState("loading");
    setError("");
    try {
      const smiles = await mol.resolve();
      if (!smiles.length) throw new Error("Enter at least one SMILES.");
      const d = await api.admet(smiles);
      if (d.mode === "job" && d.job_id) {
        const jobId = d.job_id;
        setTotal(d.total || 0);
        setState("progress");
        setPct(2);
        while (true) {
          await api.sleep(1500);
          const s = await api.pollRetry(() => api.admetJob(jobId));
          if (s.status === "done") {
            setData(s);
            setState("done");
            return;
          }
          if (s.status === "error") throw new Error("ADMET job failed.");
          setPct(Math.round((100 * (s.done || 0)) / (s.total || d.total || 1)));
        }
      } else {
        setData(d);
        setState("done");
      }
    } catch (e: any) {
      setError(e.message || "Error");
      setState("error");
    }
  };

  return (
    <SidebarLayout
      collapsed={sidebarCollapsed}
      onToggle={toggleSidebar}
      sidebar={
        <aside className="card sticky top-[78px] max-h-[calc(100vh-96px)] overflow-y-auto p-[18px]">
          <SectionIntro title="ADMET profile" sub="Drug-likeness & structural alerts. Target-independent." />
          <MoleculeInputPanel state={mol} />
          <button className="btn-primary mt-[18px]" onClick={run} disabled={state === "loading" || state === "progress"}>
            Profile compounds
          </button>
        </aside>
      }
    >
      <main className="card min-h-[60vh] overflow-hidden">
        {state === "idle" && <EmptyState title="Profile your compounds" hint="Enter molecules to compute their ADMET profile." />}
        {state === "loading" && <Spinner />}
        {state === "progress" && <ProgressBar pct={pct} label={`Profiling ${total} compounds in the ADMET-AI worker…`} />}
        {state === "error" && <ErrorBox message={error} />}
        {state === "done" && data && <AdmetResults d={data} />}
      </main>
    </SidebarLayout>
  );
}

function AdmetResults({ d }: { d: AdmetResponse }) {
  const learnedOn = d.learned?.available;
  const [open, setOpen] = useState<Set<number>>(new Set());
  const toggle = (i: number) =>
    setOpen((s) => {
      const n = new Set(s);
      n.has(i) ? n.delete(i) : n.add(i);
      return n;
    });

  return (
    <div>
      {!learnedOn && <Notice>{d.learned?.note || "Learned ADMET endpoints unavailable — showing drug-likeness layer only."}</Notice>}
      <Disclaimer>
        Drug-likeness flags are informational — natural products often violate them while remaining bioactive. They are never used to filter compounds.
        {learnedOn ? " Learned endpoints from ADMET-AI." : ""}
      </Disclaimer>
      <div className="flex justify-end px-5 pt-2">
        <button type="button" className="btn-link" onClick={() => downloadAdmetCsv(d)}>
          Download all data (CSV)
        </button>
      </div>
      <div className="max-h-[calc(100vh-260px)] overflow-y-auto overflow-x-hidden">
        <table className="w-full table-fixed border-collapse text-[13px]">
          <thead>
            <tr>
              {[
                { h: "", w: "w-6" },
                { h: "Compound", w: "w-64" },
                { h: "MW (Da)", w: "w-16", tip: "Molecular weight of the standardized structure (RDKit Descriptors.MolWt)." },
                { h: "LogP", w: "w-16", tip: "Calculated octanol-water partition coefficient (RDKit Crippen.MolLogP) — higher = more lipophilic. Dimensionless." },
                { h: "QED (0–1)", w: "w-16", tip: "Quantitative Estimate of Drug-likeness (RDKit QED.qed), combining several properties into one score — closer to 1 is more drug-like." },
                { h: "Lipinski", w: "w-20", tip: "Count of Lipinski Rule-of-Five violations (MW > 500 Da, LogP > 5, H-bond donors > 5, H-bond acceptors > 10)." },
                { h: "Alerts", w: "w-16", tip: "Count of matched structural-alert substructure patterns (PAINS / Brenk / NIH filters) — informational, never used to exclude a compound." },
                ...(learnedOn
                  ? [{ h: "Tox flags", w: "w-20", tip: "Count of ADMET-AI learned endpoints flagged as a concern for this compound (see the expanded row for which ones)." }]
                  : []),
              ].map((c, i) => (
                <th
                  key={i}
                  title={c.tip}
                  className={`sticky top-0 z-10 border-b border-line bg-surface2 px-2.5 py-2.5 text-left text-[10.5px] font-semibold uppercase tracking-wide text-inkmut ${c.tip ? "cursor-help decoration-dotted underline decoration-1 underline-offset-2" : ""} ${c.w || ""}`}
                >
                  {c.h}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {d.profiles.map((p, i) => (
              <AdmetRow key={i} p={p} i={i} learnedOn={!!learnedOn} isOpen={open.has(i)} toggle={toggle} />
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function AdmetRow({ p, i, learnedOn, isOpen, toggle }: { p: AdmetProfile; i: number; learnedOn: boolean; isOpen: boolean; toggle: (i: number) => void }) {
  if (!p.parsed_ok) {
    return (
      <tr>
        <td className="border-b border-surface2 px-3 py-2.5" />
        <td className="smi-mono border-b border-surface2 px-3 py-2.5" title={p.input_smiles}>{p.input_smiles}</td>
        <td colSpan={learnedOn ? 6 : 5} className="border-b border-surface2 px-3 py-2.5 text-inkmut">
          could not parse
        </td>
      </tr>
    );
  }
  const pc = p.physicochemical!;
  const fl = p.drug_likeness_flags!;
  const nflag = learnedOn && p.learned?.flags ? p.learned.flags.length : 0;
  const canExpand = learnedOn && !!p.learned?.groups;
  return (
    <>
      <tr className={canExpand ? "cursor-pointer hover:bg-canvas" : ""} onClick={() => canExpand && toggle(i)}>
        <td className="border-b border-surface2 px-3 py-2.5 text-brand-600">{canExpand ? (isOpen ? "▾" : "▸") : ""}</td>
        <td className="smi-mono border-b border-surface2 px-3 py-2.5" title={p.standardised_smiles}>{p.standardised_smiles}</td>
        <td className="border-b border-surface2 px-3 py-2.5">{pc.mw}</td>
        <td className="border-b border-surface2 px-3 py-2.5">{pc.logp}</td>
        <td className="border-b border-surface2 px-3 py-2.5">{pc.qed}</td>
        <td className="border-b border-surface2 px-3 py-2.5">
          <span className={`badge ${fl.lipinski_pass ? "bg-brand-500/15 text-brand-800" : "bg-amber/15 text-amber"}`}>{fl.lipinski_violations} viol</span>
        </td>
        <td className="border-b border-surface2 px-3 py-2.5">
          <span className={`badge ${p.n_alerts ? "bg-amber/15 text-amber" : "bg-brand-500/15 text-brand-800"}`}>{p.n_alerts}</span>
        </td>
        {learnedOn && (
          <td className="border-b border-surface2 px-3 py-2.5">
            <span className={`badge ${nflag ? "bg-amber/15 text-amber" : "bg-brand-500/15 text-brand-800"}`}>{nflag}</span>
          </td>
        )}
      </tr>
      {canExpand && isOpen && (
        <tr>
          <td className="border-b border-surface2 px-3 py-2.5" />
          <td colSpan={learnedOn ? 6 : 5} className="border-b border-surface2 bg-surface2/40 p-0">
            <div className="px-5 py-3.5">{epDetail(p.learned!.groups!)}</div>
          </td>
        </tr>
      )}
    </>
  );
}

function epDetail(groups: NonNullable<AdmetProfile["learned"]>["groups"]) {
  if (!groups) return null;
  return (
    <>
      {GROUP_ORDER.filter((g) => groups[g]).map((g) => (
        <div key={g} className="mb-3">
          <h5 className="mb-1.5 text-[11px] font-bold uppercase tracking-wide text-brand-700">{g}</h5>
          {groups[g].map((e, j) => {
            const tip =
              e.task === "class"
                ? `${e.label} — learned prediction from ADMET-AI (Chemprop), shown as the predicted probability of a positive/risk classification.`
                : `${e.label} — learned prediction from ADMET-AI (Chemprop), reported in ${e.unit}.`;
            return (
              <span
                key={j}
                title={e.percentile != null ? `${tip} Percentile is this compound's rank among DrugBank-approved drugs on this endpoint.` : tip}
                className={`mb-1.5 mr-1.5 inline-flex cursor-help items-center gap-1.5 rounded-lg border-l-[3px] border border-line bg-surface px-2.5 py-1.5 text-[12.5px] decoration-dotted underline-offset-4 hover:underline ${TONE_BORDER[e.tone] || "border-l-slateout"}`}
              >
                <span>{e.label}</span>
                <span className="font-bold">{e.display}</span>
                {e.percentile != null && <span className="text-[11px] text-inkmut">{e.percentile}%ile</span>}
              </span>
            );
          })}
        </div>
      ))}
    </>
  );
}
