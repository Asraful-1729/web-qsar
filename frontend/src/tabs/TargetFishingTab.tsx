import { useEffect, useRef, useState } from "react";
import * as api from "../lib/api";
import { SectionIntro, SidebarLayout, useSidebarCollapsed } from "../components/Shell";
import { EmptyState, Notice, Spinner } from "../components/Feedback";
import type { CuratedCompoundSuggestion, TargetPredictionV2Result } from "../lib/types";

/** Density-adaptive retrieval (target_prediction_v2.py) — pools evidence
    across the 10 nearest similar compounds (potency- and orthologue-
    weighted) when enough exist, else falls back to a single best match.
    An earlier v1 method (best-similarity + a heuristic evidence score)
    existed behind a method toggle here and was removed once v2 was
    validated to significantly beat it on genuinely held-out data — see
    target_prediction_v2/METHODS_AND_VALIDATION.md and
    documentation/TARGET_PREDICTION.md. */
export function TargetFishingTab() {
  const [available, setAvailable] = useState<"loading" | "no" | "yes">("loading");
  const [smiles, setSmiles] = useState("");
  const [state, setState] = useState<"idle" | "loading" | "error" | "done">("idle");
  const [error, setError] = useState("");
  const [result, setResult] = useState<TargetPredictionV2Result | null>(null);

  const [suggestions, setSuggestions] = useState<CuratedCompoundSuggestion[]>([]);
  const [suggestOpen, setSuggestOpen] = useState(false);
  const suggestTimer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const boxRef = useRef<HTMLDivElement>(null);
  const [sidebarCollapsed, toggleSidebar] = useSidebarCollapsed();

  useEffect(() => {
    api
      .targetPredictionV2Status()
      .then((s) => setAvailable(s.available ? "yes" : "no"))
      .catch(() => setAvailable("no"));
  }, []);

  // As-you-type suggestions from our own curated compound pool — matches
  // a SMILES fragment (e.g. a ring system just pasted in) or a target's
  // ChEMBL id (e.g. "CHEMBL1862" to browse a known active for that target
  // specifically). Debounced so it fires once typing pauses, not per
  // keystroke.
  useEffect(() => {
    if (suggestTimer.current) clearTimeout(suggestTimer.current);
    const q = smiles.trim();
    if (q.length < 2) {
      setSuggestions([]);
      return;
    }
    suggestTimer.current = setTimeout(() => {
      api
        .suggestCuratedCompounds(q)
        .then((r) => setSuggestions(r.results))
        .catch(() => setSuggestions([]));
    }, 300);
    return () => {
      if (suggestTimer.current) clearTimeout(suggestTimer.current);
    };
  }, [smiles]);

  useEffect(() => {
    const onDocClick = (e: MouseEvent) => {
      if (boxRef.current && !boxRef.current.contains(e.target as Node)) setSuggestOpen(false);
    };
    document.addEventListener("mousedown", onDocClick);
    return () => document.removeEventListener("mousedown", onDocClick);
  }, []);

  const selectSuggestion = (s: CuratedCompoundSuggestion) => {
    setSmiles(s.smiles);
    setSuggestOpen(false);
  };

  const run = async () => {
    if (!smiles.trim()) {
      setError("Enter a SMILES string.");
      setState("error");
      return;
    }
    setSuggestOpen(false);
    setState("loading");
    setError("");
    try {
      const r = await api.targetPredictionV2Predict(smiles.trim());
      setResult(r);
      setState("done");
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
        <SectionIntro
          title="Target prediction"
          sub="Given one compound, which protein targets is it likely to hit? Searches for structurally similar known bioactive compounds across every curated target — the reverse of disease-first browsing."
        />
        {available === "loading" && <Spinner />}
        {available === "no" && <Notice>Target Prediction data isn't available in this build.</Notice>}
        {available === "yes" && (
          <>
            <label className="field-label" style={{ marginTop: 12 }}>
              Query SMILES
            </label>
            <div ref={boxRef}>
              <input
                className="field-input font-mono text-[12.5px]"
                placeholder="e.g. a natural product, lead compound, or a target's ChEMBL id"
                value={smiles}
                onFocus={() => setSuggestOpen(true)}
                onChange={(e) => {
                  setSmiles(e.target.value);
                  setSuggestOpen(true);
                }}
              />
              {suggestOpen && suggestions.length > 0 && (
                <div className="mt-1.5 max-h-[280px] overflow-y-auto rounded-lg border border-line bg-surface">
                  {suggestions.map((s, i) => (
                    <div
                      key={i}
                      className="cursor-pointer border-b border-line/70 px-2.5 py-1.5 text-[12px] last:border-0 hover:bg-surface2/60"
                      onClick={() => selectSuggestion(s)}
                    >
                      <div className="smi-mono truncate text-ink" title={s.smiles}>
                        {s.smiles}
                      </div>
                      <div className="mt-0.5 text-[11px] text-inkmut">
                        known active for <b>{s.target_pref_name || s.target_chembl}</b>
                        {s.pchembl_value != null ? ` · pChEMBL ${s.pchembl_value}` : ""}
                      </div>
                    </div>
                  ))}
                </div>
              )}
            </div>
            <div className="field-hint">Suggestions come from the bioactivity index as you type (2+ characters).</div>
            <div className="field-hint" style={{ marginTop: 12 }}>
              Automatically decides whether to pool multiple similar compounds' evidence or fall back to a single
              best match, based on how many similar compounds exist for your query — no threshold to set.
            </div>
            <button className="btn-primary mt-[18px]" onClick={run} disabled={state === "loading"}>
              Search
            </button>
            {state === "loading" && <div className="field-hint">Searching curated bioactivity data…</div>}
            {state === "error" && <div className="field-hint text-clay">{error}</div>}
            {state === "done" && result && (
              <div className="field-hint">
                {result.results.length} target(s) ranked, out of {result.n_targets_indexed.toLocaleString()} indexed. Regime:{" "}
                {result.regime === "pooled" ? "pooled evidence" : "best-match fallback"} ({result.density} similar compound
                {result.density === 1 ? "" : "s"} found, threshold {result.density_threshold}).
              </div>
            )}
          </>
        )}
      </aside>
      }
    >
      <main className="card min-h-[60vh] overflow-hidden p-[18px]">
        {state !== "done" || !result ? (
          <EmptyState title="Target prediction" hint="Enter a compound in the sidebar to see its likely targets, with supporting evidence." />
        ) : (
          <ResultsTable result={result} />
        )}
      </main>
    </SidebarLayout>
  );
}

/** Every row carries an explicit confidence_label (never presented as a
    probability — Phase 4 measured calibration and it failed its own
    accuracy gate, see target_prediction_v2/METHODS_AND_VALIDATION.md §4)
    and full supporting evidence, since aggregate benchmark numbers alone
    aren't enough for a researcher to trust one specific prediction. */
function ResultsTable({ result }: { result: TargetPredictionV2Result }) {
  const [expanded, setExpanded] = useState<Set<string>>(new Set());
  const toggle = (t: string) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      if (next.has(t)) next.delete(t);
      else next.add(t);
      return next;
    });

  if (!result.results.length) {
    return <EmptyState title="No targets matched" hint="No similar compounds were found for this query in the index." />;
  }

  return (
    <div>
      <Notice>
        Ligand-based prediction: {result.regime === "pooled"
          ? `${result.density} similar compounds were found (≥ threshold of ${result.density_threshold}), so results pool weighted evidence across the 10 nearest compounds, potency- and orthologue-weighted.`
          : `Only ${result.density} similar compounds were found (below the threshold of ${result.density_threshold}), so results fall back to a single best-matching compound per target.`}{" "}
        The score shown is a raw, un-calibrated evidence signal (see each row's label) — not a probability that the
        compound binds this target. See target_prediction_v2/METHODS_AND_VALIDATION.md for the full validation
        account.
      </Notice>
      <div className="mt-2 text-[12px] text-inkmut">
        {result.results.length} target{result.results.length === 1 ? "" : "s"} ranked, out of{" "}
        {result.n_targets_indexed.toLocaleString()} indexed.
      </div>
      <div className="mt-1 divide-y divide-line/70">
        {result.results.map((r) => {
          const isOpen = expanded.has(r.target_chembl);
          const nativeN = r.evidence.native_neighbours ?? [];
          const orthoN = r.evidence.orthologue_neighbours ?? [];
          const hasEvidence = nativeN.length > 0 || orthoN.length > 0 || !!r.evidence.best_compound_smiles;
          return (
            <div key={r.target_chembl} className="py-3">
              <div className="flex items-center justify-between gap-3">
                <div className="min-w-0">
                  <span className="text-[11px] text-inkmut">#{r.rank}</span>{" "}
                  <span className="font-semibold text-ink">{r.target_pref_name || r.target_chembl}</span>
                  {r.target_pref_name && <span className="ml-1.5 text-[11px] text-inkmut">({r.target_chembl})</span>}
                </div>
                <span className="shrink-0 badge bg-brand-500/15 text-brand-800 cursor-help" title={r.confidence_label}>
                  score {r.score}
                </span>
              </div>
              <div className="mt-1 text-[11px] text-inkmut">{r.confidence_label}</div>
              {hasEvidence && (
                <button className="mt-1.5 text-[11.5px] font-medium text-brand-800 hover:underline" onClick={() => toggle(r.target_chembl)}>
                  {isOpen ? "Hide supporting evidence" : "Show supporting evidence"}
                </button>
              )}
              {isOpen && (
                <div className="mt-1.5 space-y-1 rounded-lg bg-surface2/40 p-2">
                  {r.evidence.best_compound_smiles && (
                    <div className="flex items-center justify-between gap-2 text-[11.5px]">
                      <span className="smi-mono truncate text-inkmut" title={r.evidence.best_compound_smiles}>
                        {r.evidence.best_compound_smiles}
                      </span>
                      <span className="shrink-0 text-inkmut">
                        best match{r.evidence.best_compound_pchembl != null ? ` · pChEMBL ${r.evidence.best_compound_pchembl}` : ""}
                      </span>
                    </div>
                  )}
                  {nativeN.map((n, i) => (
                    <div key={`n${i}`} className="flex items-center justify-between gap-2 text-[11.5px]">
                      <span className="smi-mono truncate text-inkmut" title={n.smiles}>
                        {n.smiles}
                      </span>
                      <span className="shrink-0 text-inkmut">
                        tanimoto {n.tanimoto}
                        {n.pchembl_value != null ? ` · pChEMBL ${n.pchembl_value}` : ""} · weight {n.potency_weight}
                      </span>
                    </div>
                  ))}
                  {orthoN.map((n, i) => (
                    <div key={`o${i}`} className="flex items-center justify-between gap-2 text-[11.5px]">
                      <span className="smi-mono truncate text-inkmut" title={n.smiles}>
                        {n.smiles}
                      </span>
                      <span className="shrink-0 text-inkmut">
                        tanimoto {n.tanimoto} · orthologue{n.species_provenance ? ` (${n.species_provenance})` : ""}
                      </span>
                    </div>
                  ))}
                </div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
