import { useEffect, useState } from "react";
import * as api from "../lib/api";
import { SectionIntro, SidebarLayout, useSidebarCollapsed } from "../components/Shell";
import { EmptyState, ErrorBox, Notice, Spinner } from "../components/Feedback";
import type { SimilarityResult } from "../lib/types";

export function SimilarityTab() {
  const [status, setStatus] = useState<"loading" | "unavailable" | "downloading" | "ready">("loading");
  const [dlJobId, setDlJobId] = useState<string | null>(null);
  const [dlProgress, setDlProgress] = useState<{ done: number; total: number } | null>(null);
  const [dlError, setDlError] = useState<string | null>(null);
  const [sidebarCollapsed, toggleSidebar] = useSidebarCollapsed();

  const checkStatus = async () => {
    try {
      const s = await api.similarityStatus();
      setStatus(s.available ? "ready" : "unavailable");
    } catch {
      setStatus("unavailable");
    }
  };

  useEffect(() => {
    checkStatus();
  }, []);

  const startDownload = async () => {
    setDlError(null);
    try {
      const r = await api.similarityDownloadStart();
      if (r.already_installed) {
        setStatus("ready");
        return;
      }
      if (!r.job_id) return;
      setDlJobId(r.job_id);
      setStatus("downloading");
      while (true) {
        await api.sleep(1000);
        const p = await api.similarityDownloadProgress(r.job_id);
        setDlProgress({ done: p.done, total: p.total });
        if (p.state === "done") {
          setDlJobId(null);
          setStatus("ready");
          return;
        }
        if (p.state === "error") {
          setDlError(p.error || "Download failed.");
          setDlJobId(null);
          setStatus("unavailable");
          return;
        }
        if (p.state === "cancelled") {
          setDlJobId(null);
          setStatus("unavailable");
          return;
        }
      }
    } catch (e: any) {
      setDlError(e.message || "Error");
      setStatus("unavailable");
    }
  };

  const stopDownload = async () => {
    if (!dlJobId) return;
    try {
      await api.similarityDownloadCancel(dlJobId);
    } catch {
      /* the poll loop will still surface a final status either way */
    }
  };

  return (
    <SidebarLayout
      collapsed={sidebarCollapsed}
      onToggle={toggleSidebar}
      sidebar={
        <aside className="card sticky top-[78px] max-h-[calc(100vh-96px)] overflow-y-auto p-[18px]">
        <SectionIntro
          title="Similarity search"
          sub="Upload a known drug or lead compound and find structurally similar natural products (COCONUT database, ~739K compounds)."
        />
        {status === "loading" && <Spinner />}
        {status === "unavailable" && (
          <div className="mt-3">
            <Notice>
              The similarity index (~87 MB, downloaded once) isn't installed yet.
              {dlError && <div className="mt-1.5 text-clay">{dlError}</div>}
            </Notice>
            <button type="button" className="btn-primary mt-2.5" onClick={startDownload}>
              Download similarity index
            </button>
          </div>
        )}
        {status === "downloading" && dlProgress && (
          <div className="mt-3">
            <div className="h-1.5 overflow-hidden rounded-full bg-surface2">
              <div
                className="h-full rounded-full bg-brand-500 transition-all duration-300"
                style={{ width: `${dlProgress.total ? Math.round((100 * dlProgress.done) / dlProgress.total) : 0}%` }}
              />
            </div>
            <div className="mt-1 flex items-center gap-2 text-[11.5px] text-inkmut">
              <span>
                Downloading… {dlProgress.total ? `${(dlProgress.done / 1e6).toFixed(0)} / ${(dlProgress.total / 1e6).toFixed(0)} MB` : ""}
              </span>
              <button type="button" className="font-semibold text-brand-700 underline" onClick={stopDownload}>
                Stop
              </button>
            </div>
          </div>
        )}
        {status === "ready" && <SearchForm />}
      </aside>
      }
    >
      <main className="card min-h-[60vh] overflow-hidden">
        {status !== "ready" && <EmptyState title="Similarity search" hint="Download the index in the sidebar to get started." />}
      </main>
    </SidebarLayout>
  );
}

function SearchForm() {
  const [smiles, setSmiles] = useState("");
  const [threshold, setThreshold] = useState("0.4");
  const [state, setState] = useState<"idle" | "loading" | "error" | "done">("idle");
  const [error, setError] = useState("");
  const [result, setResult] = useState<SimilarityResult | null>(null);

  const run = async () => {
    if (!smiles.trim()) {
      setError("Enter a SMILES string.");
      setState("error");
      return;
    }
    setState("loading");
    setError("");
    try {
      const r = await api.similaritySearch(smiles.trim(), parseFloat(threshold) || 0.4);
      setResult(r);
      setState("done");
    } catch (e: any) {
      setError(e.message || "Error");
      setState("error");
    }
  };

  return (
    <>
      <label className="field-label">Query SMILES</label>
      <input className="field-input font-mono text-[12.5px]" placeholder="e.g. known drug/lead compound" value={smiles} onChange={(e) => setSmiles(e.target.value)} />
      <label className="field-label" style={{ marginTop: 12 }}>
        Minimum Tanimoto similarity
      </label>
      <input type="number" min={0} max={1} step={0.05} className="field-input" value={threshold} onChange={(e) => setThreshold(e.target.value)} />
      <button className="btn-primary mt-[18px]" onClick={run} disabled={state === "loading"}>
        Search
      </button>
      {state === "loading" && <div className="field-hint">Searching ~739K natural products…</div>}
      {state === "error" && <div className="field-hint text-clay">{error}</div>}
      {state === "done" && result && (
        <div className="field-hint">
          {result.n_results} match(es) above {threshold} similarity, out of {result.n_indexed.toLocaleString()} indexed compounds.
        </div>
      )}
      {state === "done" && result && <ResultsPreview result={result} />}
    </>
  );
}

/** The sidebar only shows a compact preview + "see full results" hint —
    the real table lives in <main> via a shared context would be cleaner,
    but for a single-query tool a simple top-N preview here is enough;
    full results table renders inline in the sidebar's own scroll area
    since SimilarityTab's <main> stays a fixed empty-state otherwise. */
function ResultsPreview({ result }: { result: SimilarityResult }) {
  return (
    <div className="mt-3 max-h-[400px] overflow-y-auto rounded-lg border border-line">
      {!result.results.length && <div className="p-2.5 text-[12.5px] text-inkmut">No matches above this threshold.</div>}
      {result.results.map((r) => (
        <div key={r.id} className="border-b border-line/70 px-2.5 py-2 text-[12px] last:border-0">
          <div className="flex items-center justify-between">
            <span className="font-semibold text-ink">{r.name || r.id}</span>
            <span className="badge bg-brand-500/15 text-brand-800">{r.tanimoto}</span>
          </div>
          <div className="smi-mono mt-0.5 truncate text-inkmut" title={r.smiles}>
            {r.smiles}
          </div>
          <div className="mt-0.5 text-inkmut">
            {r.scaffold_match ? "✓ same Murcko scaffold" : ""}
            {r.chemical_super_class ? `${r.scaffold_match ? " · " : ""}${r.chemical_super_class}` : ""}
            {r.mcs_n_atoms != null ? ` · MCS ${r.mcs_n_atoms} atoms` : ""}
          </div>
        </div>
      ))}
    </div>
  );
}
