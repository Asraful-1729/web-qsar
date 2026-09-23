import { useEffect, useRef, useState } from "react";
import * as api from "../lib/api";
import { SectionIntro, SidebarLayout, useSidebarCollapsed } from "../components/Shell";
import { EmptyState, Notice, Spinner } from "../components/Feedback";
import type { CuratedCompoundSuggestion, TargetFishingResult, TargetPredictionV2Result } from "../lib/types";

/** Two independently-built, independently-validated methods live behind
    this one tab: v1 (target_fishing.py, best-similarity + a heuristic
    evidence score) and v2 (target_prediction_v2.py, density-adaptive
    pooling + potency weighting + orthologue evidence) — v2 significantly
    beats v1 on a genuinely held-out benchmark
    (target_prediction_v2/METHODS_AND_VALIDATION.md §5), so it's the
    default, but v1 stays available rather than being silently replaced. */
type Method = "v2" | "v1";

export function TargetFishingTab() {
  const [method, setMethod] = useState<Method>("v2");
  const [v1Available, setV1Available] = useState<"loading" | "no" | "yes">("loading");
  const [v2Available, setV2Available] = useState<"loading" | "no" | "yes">("loading");
  const [smiles, setSmiles] = useState("");
  const [threshold, setThreshold] = useState("0.4");
  const [state, setState] = useState<"idle" | "loading" | "error" | "done">("idle");
  const [error, setError] = useState("");
  const [result, setResult] = useState<TargetFishingResult | null>(null);
  const [resultV2, setResultV2] = useState<TargetPredictionV2Result | null>(null);

  const [suggestions, setSuggestions] = useState<CuratedCompoundSuggestion[]>([]);
  const [suggestOpen, setSuggestOpen] = useState(false);
  const suggestTimer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const boxRef = useRef<HTMLDivElement>(null);
  const [sidebarCollapsed, toggleSidebar] = useSidebarCollapsed();

  useEffect(() => {
    api
      .targetFishingStatus()
      .then((s) => setV1Available(s.available ? "yes" : "no"))
      .catch(() => setV1Available("no"));
    api
      .targetPredictionV2Status()
      .then((s) => setV2Available(s.available ? "yes" : "no"))
      .catch(() => setV2Available("no"));
  }, []);

  const available = method === "v1" ? v1Available : v2Available;

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
      if (method === "v1") {
        const r = await api.targetFishingSearch(smiles.trim(), parseFloat(threshold) || 0.4);
        setResult(r);
      } else {
        const r = await api.targetPredictionV2Predict(smiles.trim());
        setResultV2(r);
      }
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
        <label className="field-label">Method</label>
        <div className="flex gap-1.5 rounded-lg bg-surface2/60 p-1">
          <button
            className={`flex-1 rounded-md px-2 py-1.5 text-[12.5px] font-medium transition ${
              method === "v2" ? "bg-surface shadow-sm text-ink" : "text-inkmut hover:text-ink"
            }`}
            onClick={() => {
              setMethod("v2");
              setState("idle");
            }}
          >
            v2 (recommended)
          </button>
          <button
            className={`flex-1 rounded-md px-2 py-1.5 text-[12.5px] font-medium transition ${
              method === "v1" ? "bg-surface shadow-sm text-ink" : "text-inkmut hover:text-ink"
            }`}
            onClick={() => {
              setMethod("v1");
              setState("idle");
            }}
          >
            v1 (legacy)
          </button>
        </div>
        <div className="field-hint">
          {method === "v2"
            ? "Density-adaptive retrieval — significantly beats v1 and a published baseline (SEA) on held-out data. See target_prediction_v2/METHODS_AND_VALIDATION.md."
            : "Best-similarity ranking, the original method — kept available for comparison, not the current recommendation."}
        </div>
        {available === "loading" && <Spinner />}
        {available === "no" && <Notice>{method === "v1" ? "Target-fishing" : "Target Prediction v2"} data isn't available in this build.</Notice>}
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
            <div className="field-hint">Suggestions come from the target-fishing bioactivity index as you type (2+ characters).</div>
            {method === "v1" ? (
              <>
                <label className="field-label" style={{ marginTop: 12 }}>
                  Minimum Tanimoto similarity
                </label>
                <input
                  type="number"
                  min={0}
                  max={1}
                  step={0.05}
                  className="field-input"
                  value={threshold}
                  onChange={(e) => setThreshold(e.target.value)}
                />
              </>
            ) : (
              <div className="field-hint" style={{ marginTop: 12 }}>
                v2 automatically decides whether to pool multiple similar compounds' evidence or fall back to a
                single best match, based on how many similar compounds exist for your query — no threshold to set.
              </div>
            )}
            <button className="btn-primary mt-[18px]" onClick={run} disabled={state === "loading"}>
              Search
            </button>
            {state === "loading" && <div className="field-hint">Searching curated bioactivity data…</div>}
            {state === "error" && <div className="field-hint text-clay">{error}</div>}
            {state === "done" && method === "v1" && result && (
              <div className="field-hint">
                {result.n_targets_matched} target(s) matched above {threshold} similarity, out of {result.n_targets_searched} searched (
                {result.n_indexed_compound_target_pairs.toLocaleString()} compound–target pairs indexed).
              </div>
            )}
            {state === "done" && method === "v2" && resultV2 && (
              <div className="field-hint">
                {resultV2.results.length} target(s) ranked, out of {resultV2.n_targets_indexed.toLocaleString()} indexed. Regime:{" "}
                {resultV2.regime === "pooled" ? "pooled evidence" : "best-match fallback"} ({resultV2.density} similar compound
                {resultV2.density === 1 ? "" : "s"} found, threshold {resultV2.density_threshold}).
              </div>
            )}
          </>
        )}
      </aside>
      }
    >
      <main className="card min-h-[60vh] overflow-hidden p-[18px]">
        {method === "v1" ? (
          state !== "done" || !result ? (
            <EmptyState title="Target prediction" hint="Enter a compound in the sidebar to see its likely targets." />
          ) : (
            <ResultsTable result={result} />
          )
        ) : state !== "done" || !resultV2 ? (
          <EmptyState title="Target prediction v2" hint="Enter a compound in the sidebar to see its likely targets, with supporting evidence." />
        ) : (
          <ResultsTableV2 result={resultV2} />
        )}
      </main>
    </SidebarLayout>
  );
}

/** Reference evidence depth — NOT "reliability"/"confidence": that would
    imply overall confidence the target is correct, but the scaffold-split
    benchmark (TARGET_FISHING_BENCHMARK.md, 2026-09-20) only measured
    target-RECOVERY performance conditional on how many similar actives
    supported a target. Naming and bucket boundaries match that benchmark's
    own stratification exactly (1-2 / 3-10 / 11-50 / 51+), so the label is
    always traceable to a real measured number, never a vague heuristic.
    Colors reuse ConfidenceDot's existing high/medium/low/none scale
    (Feedback.tsx) for visual consistency with the rest of the app rather
    than inventing a new one. */
function referenceEvidenceInfo(n: number): { label: string; cls: string; recovery: string } {
  if (n <= 2) return { label: "Limited", cls: "bg-slateout/15 text-slateout", recovery: "0%" };
  if (n <= 10) return { label: "Low", cls: "bg-clay/15 text-clay", recovery: "21%" };
  if (n <= 50) return { label: "Moderate", cls: "bg-amber/15 text-amber", recovery: "72%" };
  return { label: "High", cls: "bg-brand-500/15 text-brand-800", recovery: "92%" };
}

type SortKey = "evidence_score" | "best_similarity" | "n_similar_actives" | "n_scaffolds" | "best_pchembl";
// Ordered to match the production default: Best similarity first — the
// scaffold-split validation benchmark (TARGET_FISHING_BENCHMARK.md,
// 2026-09-20) found it consistently beats Evidence score at actual target
// recovery, so that's now both the backend's default result order AND
// this dropdown's default selection (see useState below) — keeping both
// in sync matters because this component re-sorts client-side regardless
// of what order search() returns, so a mismatched default here would
// silently undo the backend's own ordering change on first render.
const SORT_OPTIONS: { value: SortKey; label: string }[] = [
  { value: "best_similarity", label: "Best similarity" },
  { value: "evidence_score", label: "Evidence score" },
  { value: "n_similar_actives", label: "Similar active count" },
  { value: "n_scaffolds", label: "Scaffold count" },
  { value: "best_pchembl", label: "Best pChEMBL" },
];

/** Every row is "similar known actives found for this target," not a
    calibrated probability — the disclaimer stays visible above the
    table rather than being a one-time notice, since it materially
    changes how a result should be read.

    Filtering/sorting here is deliberately ALL client-side: search() already
    returns every matched target's full evidence (n_similar_actives,
    n_scaffolds, best_pchembl, evidence_score, best_similarity) in one
    response, so narrowing or re-sorting what's already in hand needs no
    new request — unlike the sidebar's "Minimum Tanimoto similarity", which
    changes what search() itself aggregates over and genuinely needs a
    fresh search. Keeping that distinction is why these controls live here,
    in the results panel, instead of next to Tanimoto in the sidebar. */
function ResultsTable({ result }: { result: TargetFishingResult }) {
  const [minActives, setMinActives] = useState(1);
  const [minScaffolds, setMinScaffolds] = useState(1);
  const [minPchembl, setMinPchembl] = useState(0);
  const [sortBy, setSortBy] = useState<SortKey>("best_similarity");
  const [registryOnly, setRegistryOnly] = useState(false);

  if (!result.results.length) {
    return <EmptyState title="No targets matched" hint="Try lowering the similarity threshold." />;
  }

  const filtered = result.results
    .filter((r) => r.n_similar_actives >= minActives)
    .filter((r) => r.n_scaffolds >= minScaffolds)
    .filter((r) => (minPchembl > 0 ? r.best_pchembl != null && r.best_pchembl >= minPchembl : true))
    .filter((r) => (registryOnly ? !!r.target_id : true));
  const sorted = [...filtered].sort((a, b) => {
    const av = a[sortBy], bv = b[sortBy];
    if (av == null && bv == null) return 0;
    if (av == null) return 1;    // unknown (e.g. no pChEMBL data) sorts last regardless of direction
    if (bv == null) return -1;
    return bv - av;              // every sort field: higher = more/stronger evidence
  });

  return (
    <div>
      <Notice>
        Ligand-based prediction: results are targets with structurally similar known active compounds in a broad
        ChEMBL bioactivity pull, not a trained multi-target classifier or a calibrated probability. Ranked by best
        similarity — a scaffold-split validation benchmark (TARGET_FISHING_BENCHMARK.md) found this recovers true
        targets more reliably than the Evidence score shown alongside it; Evidence score is a heuristic summary of
        supporting chemical evidence, not a validated target-prediction score. Treat as a starting hypothesis to
        confirm with Docking or Screen, not a final answer.
      </Notice>
      <div className="mt-3 flex flex-wrap items-end gap-x-4 gap-y-2 border-b border-line pb-3">
        <div>
          <label className="field-label">Sort by</label>
          <select className="field-input h-[30px] py-0 text-[12.5px]" value={sortBy} onChange={(e) => setSortBy(e.target.value as SortKey)}>
            {SORT_OPTIONS.map((o) => (
              <option key={o.value} value={o.value}>
                {o.label}
              </option>
            ))}
          </select>
        </div>
        <div>
          <label className="field-label">Min. similar actives</label>
          <input
            type="number"
            min={1}
            className="field-input h-[30px] w-[70px] py-0 text-[12.5px]"
            value={minActives}
            onChange={(e) => setMinActives(Math.max(1, parseInt(e.target.value) || 1))}
          />
        </div>
        <div>
          <label className="field-label">Min. scaffolds</label>
          <input
            type="number"
            min={1}
            className="field-input h-[30px] w-[70px] py-0 text-[12.5px]"
            value={minScaffolds}
            onChange={(e) => setMinScaffolds(Math.max(1, parseInt(e.target.value) || 1))}
          />
        </div>
        <div>
          <label className="field-label">Min. best pChEMBL</label>
          <input
            type="number"
            min={0}
            step={0.5}
            className="field-input h-[30px] w-[70px] py-0 text-[12.5px]"
            value={minPchembl}
            onChange={(e) => setMinPchembl(Math.max(0, parseFloat(e.target.value) || 0))}
          />
        </div>
        <label className="mb-1.5 flex cursor-pointer items-center gap-1.5 text-[12.5px] font-medium text-ink">
          <input type="checkbox" checked={registryOnly} onChange={(e) => setRegistryOnly(e.target.checked)} />
          Only targets in app registry
        </label>
      </div>
      <div className="mt-2 text-[12px] text-inkmut">
        Showing {sorted.length} of {result.results.length} matched target{result.results.length === 1 ? "" : "s"}.
      </div>
      <div className="mt-1 divide-y divide-line/70">
        {!sorted.length && <div className="py-3 text-[12.5px] text-inkmut">No targets match these filters — try loosening one.</div>}
        {sorted.map((r) => (
          <div key={r.target_chembl} className="py-3">
            <div className="flex items-center justify-between gap-3">
              <div className="min-w-0">
                <span className="font-semibold text-ink">{r.target_pref_name || r.target_chembl}</span>
                {!r.target_id && (
                  <span className="ml-1.5 text-[11px] text-inkmut" title="This ChEMBL target isn't one of PhytoScreen's own docking/QSAR targets yet.">
                    ({r.target_chembl})
                  </span>
                )}
              </div>
              <span
                className="shrink-0 badge bg-brand-500/15 text-brand-800 cursor-help"
                title="Heuristic summary of supporting chemical evidence; not a probability of target binding and not validated as a target-prediction score. See the individual similarity/scaffold/pChEMBL numbers below for the actual supporting evidence."
              >
                Evidence score {r.evidence_score}
              </span>
            </div>
            <div className="mt-0.5 flex flex-wrap items-center gap-x-3 gap-y-0.5 text-[12px] text-inkmut">
              <span>
                {r.n_similar_actives} similar known active{r.n_similar_actives === 1 ? "" : "s"}
              </span>
              {(() => {
                const info = referenceEvidenceInfo(r.n_similar_actives);
                return (
                  <span
                    className={`badge cursor-help ${info.cls}`}
                    title={`Reference evidence depth: this target is supported by ${r.n_similar_actives} similar active compound${r.n_similar_actives === 1 ? "" : "s"} in the ChEMBL reference set. In our scaffold-held-out benchmark, Top-10 target recovery was 0% for 1-2 reference actives, 21% for 3-10, 72% for 11-50, and 92% for 51+ (this target: ${info.recovery} band). These are benchmark results, not probabilities that this target binds your query compound.`}
                  >
                    {info.label} reference evidence
                  </span>
                );
              })()}
              <span>{r.n_scaffolds} distinct scaffold{r.n_scaffolds === 1 ? "" : "s"}</span>
              <span>
                best/mean similarity {r.best_similarity}/{r.mean_similarity}
              </span>
              {r.best_pchembl != null && (
                <span>
                  best/mean pChEMBL {r.best_pchembl}
                  {r.mean_pchembl != null ? `/${r.mean_pchembl}` : ""}
                </span>
              )}
            </div>
            <div className="mt-1.5 space-y-1">
              {r.compounds.map((c, i) => (
                <div key={i} className="flex items-center justify-between gap-2 text-[11.5px]">
                  <span className="smi-mono truncate text-inkmut" title={c.smiles}>
                    {c.smiles}
                  </span>
                  <span className="shrink-0 text-inkmut">
                    tanimoto {c.tanimoto}
                    {c.pchembl_value != null ? ` · pChEMBL ${c.pchembl_value}` : ""}
                  </span>
                </div>
              ))}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

/** Target Prediction v2's results — density-adaptive retrieval, not the
    best-similarity ranking above. Every row carries an explicit
    confidence_label (never presented as a probability — Phase 4 measured
    calibration and it failed its own accuracy gate, see
    target_prediction_v2/METHODS_AND_VALIDATION.md §4) and full supporting
    evidence, since aggregate benchmark numbers alone aren't enough for a
    researcher to trust one specific prediction. */
function ResultsTableV2({ result }: { result: TargetPredictionV2Result }) {
  const [expanded, setExpanded] = useState<Set<string>>(new Set());
  const toggle = (t: string) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      if (next.has(t)) next.delete(t);
      else next.add(t);
      return next;
    });

  if (!result.results.length) {
    return <EmptyState title="No targets matched" hint="No similar compounds were found for this query in the v2 index." />;
  }

  return (
    <div>
      <Notice>
        Ligand-based prediction (v2): {result.regime === "pooled"
          ? `${result.density} similar compounds were found (≥ threshold of ${result.density_threshold}), so results pool weighted evidence across the 10 nearest compounds, potency- and orthologue-weighted.`
          : `Only ${result.density} similar compounds were found (below the threshold of ${result.density_threshold}), so results fall back to a single best-matching compound per target.`}{" "}
        The score shown is a raw, un-calibrated evidence signal (see each row's label) — not a probability that the
        compound binds this target. Validated to significantly outperform both the legacy (v1) method and a
        published baseline (SEA) on held-out data — see target_prediction_v2/METHODS_AND_VALIDATION.md.
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
