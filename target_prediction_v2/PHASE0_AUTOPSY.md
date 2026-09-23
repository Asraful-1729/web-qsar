# Target Prediction v2 — Phase 0 Autopsy

**Status:** G0 gate — complete (first pass, sample-based; see "What Phase 0 did NOT cover" for the honest scope limits).
**Ran against:** `target_fishing_v1_freeze/` (the frozen, checksummed v1 snapshot) — v1 itself, `backend/target_fishing.py` and its live index, was never touched.
**Date:** 2026-09-21.
**Procedure:** `target_prediction_v2_procedure_rev3.md`, Section 4.

---

## 0. What this is

Section 4 of the procedure doc asks for a measurement pass, not a rebuild: run v1 as-shipped against drug-realistic queries it was never validated on, find out *why* it might underperform published work, and fill in the decision table before writing a line of v2 architecture. This document is that pass, plus the code and data that produced it (all under `target_prediction_v2/phase0/`, entirely separate from the frozen snapshot and from production).

**Everything below is real, measured output** from the code in this directory — no numbers were estimated or reasoned out in place of running the pipeline.

---

## 1. Infrastructure built

| Artifact | Purpose |
|---|---|
| `phase0/fetch_chembl_reference.py` | Pulls ChEMBL `molecule.json` (max_phase 1–4, 12,785 rows) and `mechanism.json` (6,984 usable rows) — external data v1's own index doesn't carry (no `molecule_chembl_id`, no drug-mechanism annotation). |
| `phase0/build_eval_sets.py` | Standardizes the ChEMBL drug pull the *same way* v1's index build does, matches by SMILES against the frozen index, and builds Sets A, B, D with the Section 3.1 three-level ground truth (annotated / mechanistic / primary target). |
| `phase0/build_leakage_fp2.py` | Precomputes a second, independent fingerprint (RDKit topological, 2048-bit) for all 857,232 distinct indexed compounds — used only for leakage detection, never for ranking, so a near-duplicate that Morgan folds/collides similarly on doesn't slip through. |
| `phase0/leakage.py` | Per-query leakage control: removes the query and any compound at Tanimoto ≥ cutoff under *either* fingerprint from the reference before searching. |
| `phase0/metrics.py` | Section 3.2 metrics: any-annotated and primary-target Top-k/MRR/median-rank, precision/recall/MCC/targets-per-query, precision-coverage curve. |
| `phase0/run_phase0.py` | Orchestrates all of the above into the actual measurement + experiment runs below. |

---

## 2. Evaluation sets (Section 4, point 1)

Matched by standardized SMILES against v1's 857,232-compound index (v1 stores no `molecule_chembl_id`, so this is the only way to identify "is this drug even in the reference pool").

| Set | Definition | Pulled from ChEMBL | Matched into v1 index | Match rate |
|---|---:|---:|---:|---:|
| **A — Approved** | max_phase = 4 | 4,225 | **1,965** | 46.5% |
| **B — Clinical** | max_phase 1–3 | 8,559 | **2,750** | 32.1% |
| **D — Sparse-target** | subset of A∪B where the drug's *weakest-covered* annotated target has ≤10 supporting compounds | — | **334** | — |
| **C — Scaffold-held-out** | novel chemistry, not re-built — reused verbatim from the already-validated `target_fishing_v1_freeze/target_fishing_benchmark_report.json` (7,295 instances) | — | — | — |

**Honest coverage gap, stated directly:** more than half of ChEMBL's approved drugs (53.5%) and two-thirds of clinical-phase compounds (67.9%) have **no bioactivity row at all** passing v1's own index filter (human single-protein target, IC50/Ki/Kd/EC50, assay_confidence≥8, pChEMBL present). That's either a real data-coverage limit or a salt/tautomer standardization mismatch between the ChEMBL drug record and however v1 indexed the same structure — Phase 0 did not separate these two causes (flagged as a Phase 1 follow-up, not papered over here).

**Primary/mechanism-target ground truth (Section 3.1, "Decisions needed" #2):** only **1,469 of 4,715** matched drugs (31.2%) have a ChEMBL-curated mechanism-of-action target. Primary-target metrics below are computed on that 31.2% subset — real, but partial; the other ~69% only have any-annotated ground truth available.

**Set B naming caveat:** Section 4 calls this set "clinical *and preclinical*." ChEMBL's `max_phase` field has no defined value for genuine preclinical candidates — the vast majority of the database is just `max_phase` null, which is ordinary screening chemistry, not preclinical drug candidates. Set B here is **clinical-only** (phases 1–3); a real preclinical set needs a different, curated source, out of scope for this pull.

---

## 3. H1 — is the benchmark too easy? (tested directly)

V1's own frozen scaffold-split benchmark landed 98% of its evaluation instances in the "51+ known actives" (High) reference-evidence bucket. Measuring the same bucket distribution on drug-realistic queries:

| Reference-evidence bucket | Set A (approved, n=1,965) | Set B (clinical, n=2,750) |
|---|---:|---:|
| 1–2 (Limited) | 2.0% | 2.0% |
| 3–10 (Low) | 6.7% | 3.9% |
| 11–50 (Moderate) | 19.5% | 10.5% |
| 51+ (High) | **71.8%** | **83.5%** |

**Finding: H1 partially confirmed, moderately.** Real drugs land in the "easy" bucket less often than the original benchmark's 98% (72–84% vs. 98%), so the original benchmark was measuring a somewhat easier population — but it's a moderate gap, not an extreme one. 16–28% of drug queries genuinely fall into harder buckets. This matches the review's Section 1 decision to soften "the algorithm class is not the problem" rather than accept or fully reject it.

---

## 4. Baseline measurement: v1-as-shipped on drug-realistic queries

Ranked by `best_similarity` (v1's frozen production default), leakage-controlled at Tanimoto ≥0.95, one search per query at threshold=0.0 (so rank is always defined), n=200 (A), 200 (B), 334 (D, full matched set).

| | Set A (approved) | Set B (clinical) | Set D (sparse-target*) | Set C (scaffold-holdout, reused) |
|---|---:|---:|---:|---:|
| Any-annotated Top-1 | 46.0% [39–53%] | 70.0% [63–76%] | 66.5% [61–71%] | 56.4% |
| Any-annotated Top-10 | 75.0% [69–80%] | 90.5% [86–94%] | 86.8% [83–90%] | 92.2% |
| Any-annotated MRR | 0.559 | 0.768 | 0.737 | 0.705 |
| Primary-target Top-1 (n with MoA target) | 48.4% (n=64) | 64.9% (n=57) | 58.2% (n=141) | — (no primary-target tier in C) |
| Primary-target Top-10 | 79.7% | 94.7% | 89.4% | — |
| Primary-target MRR | 0.600 | 0.758 | 0.695 | — |

*Set D caveat (methodological limitation, flagged not hidden): "any-annotated" counts a hit if **any** of a drug's annotated targets is recovered — for a drug with one sparse target and several well-populated ones, a well-populated target can carry the "hit," masking whether the sparse one specifically was found. Set D's numbers above should **not** be read as "sparse targets are well-recovered" — that requires a target-specific metric Phase 0 didn't build. Flagged for Phase 1.

**Notable, unplanned finding:** Set B (clinical-phase) recovers *better* than Set A (approved) on every metric, not worse. Plausible explanation: approved drugs skew toward older chemotypes (including natural products) with fewer modern ChEMBL analogs, while active clinical programs generate dense series of structurally similar compounds against the same target — the opposite of what "drugs are harder" (H1/[2]) would predict on its own. Worth carrying into Phase 1's benchmark design rather than assuming approved-drug difficulty is uniform.

---

## 5. Reference comparison ([1]: 10-NN Morgan, clinical drugs, 0.348 precision / 0.423 recall / ~11.7 targets/query)

Set A's own precision-coverage curve, at the operating point closest to [1]'s breadth (~11.7 targets/query):

| Threshold | Precision | Recall | Targets/query | MCC |
|---|---:|---:|---:|---:|
| 0.4 (v1 default) | 0.099 | 0.545 | 39.6 | 0.230 |
| **0.5** | **0.233** | **0.421** | **13.1** | 0.312 |
| 0.6 | 0.354 | 0.306 | 6.2 | 0.328 |

**Finding: at matched breadth (~13 targets/query), v1's recall (0.421) is essentially equal to [1]'s reference (0.423) — but precision (0.233) is well below it (0.348), about a third lower.** This is informative and specific: v1 is not casting a narrower or less-recalling net than published work at the same breadth; it is admitting more false positives per true positive. That points toward a ranking/filtering problem (consistent with H3/H4 below) rather than a fundamental ceiling problem (the strong form of H1 the review already downgraded).

---

## 6. H3 — potency floor (tested directly)

**Confirmed, concretely, before running any experiment:** v1's index has no potency floor beyond "pChEMBL present." Of 1,312,849 indexed pairs: **110,289 (8.4%) are weaker than 10 µM**, **383,808 (29.2%) weaker than 1 µM**, minimum pChEMBL 1.05.

Potency-filter experiment (Set A, n=200, reduced reference restricted to pChEMBL≥floor):

| | Baseline (no floor) | ≥10 µM floor | ≥1 µM floor |
|---|---:|---:|---:|
| Any-annotated Top-10 | 75.0% | 71.5% (−3.5) | 68.0% (−7.0) |
| Any-annotated MRR | 0.559 | 0.560 (≈0) | 0.544 (−0.015) |
| Primary-target Top-1 (n=64) | 48.4% | 54.7% (+6.3) | 51.6% (+3.2) |
| Primary-target MRR | 0.600 | 0.642 (+0.042) | 0.627 (+0.027) |

**Finding: H3 is real but the effect is modest and two-directional, not a dominant fix.** Filtering weak potency *helps* primary-target ranking (removing noisy weak evidence sharpens the top of the list for the real mechanism target) but *hurts* any-annotated breadth-recall (some real, weaker-potency annotated off-targets get filtered out along with the noise). This matches Section 4's decision-table framing only partially — supports making potency threshold a tunable hyperparameter (Phase 2 already plans this) rather than picking one fixed cutoff as a silver bullet.

---

## 7. H4 — weighted k-NN vs. best_similarity (tested directly)

`score(t) = Σ over k nearest neighbour compounds i whose annotation set contains t of tanimoto_i^α` (the review's formula-fix, Section 1 point 5), swept over (k, α), same leakage-controlled queries as the baseline.

**Set A, primary-target metrics (n=64):**

| | best_similarity (frozen) | k=10, α=1.0 | k=25, α=1.0 | k=50, α=1.0 |
|---|---:|---:|---:|---:|
| Top-1 | 48.4% | 60.9% (+12.5) | 62.5% (+14.1) | **65.6% (+17.2)** |
| Top-5 | 73.4% | 78.1% | 79.7% | 79.7% |
| MRR | 0.600 | 0.707 | 0.709 | **0.719** |
| Median rank when found | 2 | 1 | 1 | 1 |
| % never found | 0.0% | 7.8% | 7.8% | 6.2% |

**Set B, primary-target metrics (n=57):**

| | best_similarity | k=10, α=0.5/1.0 | k=25/50, α=1.0 |
|---|---:|---:|---:|
| Top-1 | 64.9% | **80.7% (+15.8)** | 79.0% (+14.1) |
| MRR | 0.758 | 0.852 | 0.846 |
| % never found | 0.0% | 7.0% | 3.5–5.3% |

**Finding: H4 is clearly confirmed and is the strongest, cleanest result in this Phase 0 pass.** Weighted k-NN beats best_similarity by a large, consistent margin on primary-target Top-1 and MRR across both drug sets (+12 to +17 points absolute on Top-1, ~+15–20% relative on MRR) — at the cost of a real but modest coverage gap (6–8% of primary-target queries return no hit at all under kNN, vs. 0% under best_similarity, because kNN only considers targets annotated to the top-k neighbour *compounds*, while best_similarity ranks every target with any nonzero hit). Any-annotated metrics move much less — the gain is concentrated exactly where the review predicted (the harder, more specific "what is THE target" task), not the easier "any known target" task. This is directly actionable for Phase 1/3A: weighted k-NN, not best_similarity, should be the new baseline ranking to beat.

---

## 8. Leakage sensitivity (Section 4, point 2 — review point 4)

First pass compared non-identical query subsets and produced a confusing (backwards-looking) result — re-run as a proper **paired** comparison, same exact 100 queries at both cutoffs:

| | Tanimoto ≥0.95 (looser removal) | Tanimoto ≥0.90 (stricter removal) |
|---|---:|---:|
| Any-annotated Top-1 | 49.0% | 46.0% |
| Any-annotated Top-10 | 78.0% | 78.0% |
| Primary-target Top-1 (n=31) | 58.1% | 64.5% |
| Primary-target MRR | 0.695 | 0.733 |

**Finding: the leakage-cutoff choice makes essentially no difference** — every gap here is well within sampling noise for n=31–100 (binomial SE ≈9 points at this n). v1's numbers are not an artifact of near-duplicate leakage in the 0.90–0.95 similarity band. This is a genuine negative result, reported as one rather than dropped.

---

## 9. Popularity baseline (sanity floor)

Rank every target by raw global popularity, ignore the query molecule entirely (n=734, A+B+D combined):

| | Any-annotated | Primary-target (n=262) |
|---|---:|---:|
| Top-1 | 6.5% | 3.4% |
| Top-10 | 30.5% | 14.9% |
| MRR | 0.163 | 0.086 |

v1 (any configuration tested) beats this floor by a wide margin on every metric — confirms the method is doing real structure-based work, not just recovering popular targets.

---

## 10. Decision table (Section 4, point 9)

| Row | Applies? |
|---|---|
| Potency filter or weighted k-NN closes **most** of the gap to [1] | **Partially.** Weighted k-NN closes a large fraction of the primary-target gap (H4, clean win) but not the whole gap, and trades away some coverage. Potency filtering is real but small and two-directional (H3). Neither alone is a full fix. |
| v1 near [1] but users are unhappy (expectation/presentation problem) | **Partially.** Recall matches [1] almost exactly at matched breadth; precision lags — this reads as a real ranking-quality gap, not purely a presentation gap, though H2 (ground-truth looseness/breadth) likely also contributes to the precision shortfall and needs Phase 1 investigation. |
| Most queries below 0.4 similarity (coverage limit) | **Not tested this pass** — would need the actual per-query achieved-similarity distribution, not built here. Flagged for Phase 1. |
| Servers win mainly on particular classes (data gap) | **Not tested** — no SwissTargetPrediction/SEA comparison was run (see below). |

**Net read:** this is not a single-row outcome. The evidence points to a **combination of a fixable ranking method (H4, strong, adopt weighted k-NN) and a partially-fixable evidence-quality issue (H3, real but modest)**, on top of a benchmark that was moderately (not drastically) easier than drug-realistic use (H1). Recommend proceeding into Phase 1 with weighted k-NN as the new baseline to beat, potency threshold as an explicit hyperparameter, and the ground-truth completeness question (H2) as a named early Phase 1 task rather than an assumption.

---

## 11. What Phase 0 did NOT cover (explicit, not silently skipped)

- **Sample, not census.** A(200/1965), B(200/2750) sampled at seed=42; D and the leakage-sensitivity pairing used full/matched subsets. Scaling to the full matched population is mechanical (raise the `SAMPLE_*` constants in `run_phase0.py`, re-run) — each query costs ~4–9s, so a full census is a multi-hour background job, not run in this pass.
- **H2 (ground-truth completeness/looseness)** — flagged as a likely contributor to the precision gap (Section 5) but not directly measured; needs its own investigation (e.g., how many of v1's "false positives" at matched breadth are actually true, unannotated polypharmacology vs. genuine noise).
- **H5 (0.4 floor / novel-compound coverage)** — Set D as built conflates "has a sparse target" with "is recovered via that sparse target" (see caveat in Section 4); needs a target-specific metric before drawing a conclusion.
- **H6 (fingerprint/branch ablation ladder)**, **H7 (reliability layer)**, **H8 (data-scope stratification, e.g. bacterial targets)** — none tested; out of scope for Phase 0's autopsy focus, belong to Phase 1/3/4 per the procedure doc's own phase boundaries.
- **External server comparison (SwissTargetPrediction, SEA)** — no public bulk API for either; needs a small manual/scripted-per-compound run on a common query set, not attempted here (this is Gate G5's job in a later phase, not a Phase 0 requirement per Section 4 vs. Section 7).
- **Family-level and within-family-discrimination metrics** — ChEMBL's `target_component`/`protein_classification` data exists and is pullable, but the target→family mapping wasn't built this pass. Explicitly diagnostic-only per Section 7 (not a gate), deferred without blocking G0.
- **Calibration metrics (Brier/ECE)** — no calibrated scorer exists yet to evaluate; correctly Phase 4's job, not Phase 0's.

---

## 12. Reproducing / extending

All code is in `target_prediction_v2/phase0/`, read-only against `target_fishing_v1_freeze/`:

```bash
PY=/home/storage/jonaid/projects/qsar-desktop/backend/.venv/bin/python
cd target_prediction_v2/phase0
$PY fetch_chembl_reference.py      # ~1 min, pulls ChEMBL drug + mechanism data
$PY build_eval_sets.py             # ~20s, builds Sets A/B/D + ground truth
$PY build_leakage_fp2.py           # ~25 min, one-time second-fingerprint precompute
$PY run_phase0.py                  # ~70-90 min at current SAMPLE_* sizes
```

Raw results: `phase0/results/phase0_report.json`, `phase0/results/leakage_sensitivity_paired_095.json`.
