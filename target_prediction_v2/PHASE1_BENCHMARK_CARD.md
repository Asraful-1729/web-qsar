# Target Prediction v2 — Benchmark Card

**`BUILD_PLAN.md` Phase 1 item 4.** Consolidates every real, measured fact about the benchmark this program has built and validated so far — nothing here is aspirational or a description of intent; every number traces to a specific prior document or a script in `phase0/` / `phase0b/` / `phase1/`. Read alongside `PHASE1_HOLDOUTS_AND_SPLITS.md` for the partition mechanics.

---

## 1. Data source

ChEMBL, pulled via the public REST API (`www.ebi.ac.uk/chembl/api/data`) at the dates recorded in each pull's own manifest — no single "ChEMBL version number" was pinned at the start of this program (a real gap, noted in §8). The frozen v1 index (`target_fishing_v1_freeze/`) is the fixed reference point every measurement in this program runs against; it was never modified.

## 2. Ground-truth recipe — as currently implemented (not yet Phase 2's rebuilt version)

This is the recipe the **existing v1 index** was built with, which every Phase 0/0b/1 measurement in this program has used. It differs from rev5 §4.2's target recipe in ways already measured and disclosed — Phase 2 owns closing these gaps, not this card:

| Recipe element | Current state | Source |
|---|---|---|
| Bioactivity types | IC50, Ki, Kd, EC50 only — no PubChem `Potency` | Confirmed by construction (index-build filter), `PHASE0_AUTOPSY.md` §6 (L3) |
| Relations | **Only `=`** — zero censored (`<`/`<=`) records exist anywhere in the index | Confirmed empirically: 0 ChEMBL records anywhere have `pchembl_value` set with a non-`=` relation (a hard ChEMBL platform convention, not a filter choice) — `PHASE0_AUTOPSY.md` §6 (L4) |
| Assay confidence | `assay_confidence_score ≥ 8` (hard filter, not a weight) | `V1_FREEZE.md` |
| Potency floor | **None beyond "pChEMBL present."** 8.4% of indexed pairs weaker than 10µM, 29.2% weaker than 1µM, min pChEMBL 1.05 | `PHASE0_AUTOPSY.md` §6 (H3) |
| Salt/parent collapse | Not explicit — standardization strips salts to the largest fragment at index-build time (same effect for search purposes), but molecule-level ChEMBL identity (which parent a salt form maps to) isn't separately tracked | `build_target_fishing_index.py` |
| Species scope | 100% human, by construction | `V1_FREEZE.md` |
| Document/assay metadata per pair | **Not retained** for the full ~1.3M-pair reference (only for the small, query-scoped mechanism-support pulls) | `PHASE0B_ADDENDUM.md` §13, this card §5 |

**Explicit statement, per rev5 §6.5**: because this recipe is stricter in some ways (zero censored records, hard confidence filter) and looser in others (no potency floor) than [1]'s own recipe, **our figures are not directly comparable to published numbers without accounting for these differences** — this is why the program's internally-reproduced baselines (§6) are the trustworthy comparators, not the external literature figures.

## 3. |T| (supported targets per query), by set and tier

Full matched population, computed directly from `eval_sets.json` (not sample-limited):

| Set | n | Tier | Mean \|T\| | Median \|T\| |
|---|---:|---|---:|---:|
| A (approved) | 1,965 | annotated | 8.81 | 3 |
| A | 1,965 | primary/intended | 0.46 | 0 |
| B (clinical) | 2,750 | annotated | 6.54 | 2 |
| B | 2,750 | primary/intended | 0.38 | 0 |
| D (sparse-target) | 334 | annotated | 33.84 | 12 |
| D | 334 | primary/intended | 0.77 | 0 |

(Full histograms in `PHASE0_AUTOPSY.md`/Phase 0b working data.)

## 4. Target universe

**4,658 human single-protein targets** (v1's frozen scope, unchanged). Orthologue scoping (not yet built into the index): 2,600 candidate orthologue targets across mouse/rat/cow/pig/rabbit/dog/guinea pig/macaque/chicken, covering 32.6% of the 4,658 human targets by name-match proxy (`PHASE0B_ADDENDUM.md` §10). No target-family/class breakdown built yet — diagnostic-only per rev5 §7, not a gate.

## 5. Ground-truth coverage — primary vs. mechanistically-supported

| Tier | Scope | Compound-level coverage |
|---|---|---:|
| Primary/intended | **Full matched population** (4,715 drugs, A+B combined) | **31.2%** (1,469/4,715) — set-specific: A 34.6%, B 28.7%, D 42.2% |
| Mechanistically supported (2-of-3 rule) | **Full matched population** for A/B (4,715 drugs; D was already full-population, n=334, from Phase 0b) | A **50.4%** (990/1,965), B **57.2%** (1,573/2,750), D 74.3% |

**Item 7 — CLOSED.** Full-population recompute done: `phase1/fetch_mechanism_support_full.py` (batched, recipe-filtered ChEMBL pull, 4,715/4,715 molecules, 0 failures, 2,012s) → `phase1/compute_mechanism_coverage_full.py` (applies the identical 2-of-3 rule from `PHASE0B_ADDENDUM.md` §13) → `phase1/data/mechanism_coverage_full.json`.

**Correction to the earlier sample-based recommendation, disclosed plainly**: Set A's full-population figure (50.4%) essentially matches its 200-drug sample estimate (50.5%) — the sample was representative. **Set B's full-population figure (57.2%) is meaningfully lower than its sample estimate (63.5%, a 6.3pt drop)** — B no longer clearly clears the ~60% adequacy bar at full population scale; it sits just under it. Pair-level coverage (of pairs with pulled data, matching §13's denominator convention): A 33.7% (4,112/12,217), B 34.2% (5,511/16,126) — closely matched between sets, unlike the compound-level figures, meaning the earlier B-vs-A gap was concentrated in *how many drugs have at least one supported target*, not in the underlying per-pair support rate.

**Adequacy read (rev5 Decision 3's ~60% bar), revised**: primary/intended coverage (31.2% overall) is **below** the adequacy bar for all sets, unchanged. Mechanistically-supported coverage clears it clearly for **D only** (74.3%); **A (50.4%) and B (57.2%) both now sit below or at the margin of the bar** at full-population scale — this is a real downgrade from the earlier sample-based read, not a data-quality artifact (data-pull coverage was 70.6% of A's pairs and 89.6% of B's, both high). **Revised recommendation**: treat mechanistically-supported as the powered ground-truth tier for D only; report A's and B's primary/mechanistic metrics with an explicit under-coverage caveat rather than treating B as adequately powered.

## 6. Averaging convention

**Macro (per-query, then averaged) is primary**, matching [1] and [3] — rev5's own mandate, confirmed necessary by measurement: macro vs. pooled differ by up to ~500% relative depending on threshold/tier (`PHASE0B_ADDENDUM.md` §2). Pooled is reported alongside as a secondary, always labeled. Every table in this program that reports precision/recall states which convention it's using — this is now house style, not an occasional caveat.

**One number requiring a standing caveat, not a headline claim**: at [1]'s own operating point (top-10, no cutoff), Set A's macro recall (0.474) numerically exceeds [1]'s published figure (0.423). Per `PHASE0B_REV6_RESPONSE.md` §2.4, this is very likely a **ground-truth-scope effect** (our broadened relations/assay handling converting previously-uncounted predictions into true positives), **not a claimed modelling win** — must be stated this way wherever this number is used externally, never as "beats the published state of the art."

## 7. Leakage protocol — all four levels, current status

| Level | Status | Finding |
|---|---|---|
| Near-duplicate (Tanimoto ≥0.90/0.95, two independent fingerprints) | **Closed** | Negligible sensitivity to the 0.90-vs-0.95 cutoff choice, confirmed via a proper paired comparison (`PHASE0_AUTOPSY.md` §8) |
| Scaffold-level (whole analogue series removed) | **Closed, now mandatory** | Material effect: −4pts Top-1, −7pts Top-10 vs. near-duplicate-only control (`PHASE0_AUTOPSY.md` §8). **This is the control every H4a/H4b/density result in this program has used since Rev 6.** |
| Document/assay-campaign | **Open, correctly deferred to Phase 2** | v1's aggregated index has no document_chembl_id for the full reference pool; not fakeable without a real data pull, which Phase 2's rebuild already plans |
| Temporal | **Open, correctly deferred** | No release/first-seen date retained per pair currently |

## 8. Evaluation sets — final sizes, this card's authoritative numbers

| Set | Full matched population | Sample used for density/H4 work (n=800 each, A/B) |
|---|---:|---:|
| A (approved) | 1,965 | 800 |
| B (clinical, phases 1–3) | 2,750 | 800 |
| D (sparse-target subset) | 334 | 334 (all) |
| C (scaffold-held-out, reused from v1's own benchmark) | 7,295 evaluation instances, 5,001 distinct compounds | — |

**Set C drug-realism check (P1-C), run this pass**: Set C's held-out compounds have **near-zero direct overlap with Sets A/B** (21 of 5,001 compounds, 0.42%) — it is drawn from the general ChEMBL screening pool, not restricted to real drugs. However, a random 2,000-compound sample shows a genuinely **drug-like physicochemical profile** (mean MW 477, median MW 447, **84.7% Lipinski pass rate**) — comparable to typical medicinal-chemistry compound sets, not arbitrary molecules. **Verdict: Set C is drug-*like*, not composed of real drugs.** Acceptable to keep for its stated purpose (novel-scaffold generalization), but the benchmark card and any downstream reporting must describe it precisely as "drug-like ChEMBL screening compounds," never as "real drugs" or "approved drugs."

**New this pass — the locked scaffold test and tuning/validation partition** (`PHASE1_HOLDOUTS_AND_SPLITS.md`): built from the full 3,975-compound de-duplicated A+B population, seeded (42), scaffold-group-safe.

- **Locked scaffold test: 795 compounds, 598 scaffold groups — sealed, not touched again until G2.**
- Tuning: 2,385 compounds (scaffold-split) / 2,385 (random-split, same size, different assignment).
- Validation: 795 / 795.

## 9. Baselines — see `PHASE1_BASELINES.md` for the full reference table.

## 10. What this card does not yet cover

- Full-population mechanistically-supported coverage (§5, in progress).
- The renewable temporal holdout and document/assay-campaign split type (§7 — correctly deferred to Phase 2).
- The P1-A adjudication study (§ see `PHASE1_ADJUDICATION_STUDY_STATUS.md` — requires a human domain reviewer, not something this program can generate itself).
- A pinned ChEMBL release version number for the whole program (a real, disclosed gap — different pulls across this program's history were not all timestamped against one fixed release).
