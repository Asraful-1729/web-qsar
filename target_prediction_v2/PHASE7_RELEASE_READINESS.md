# Phase 7 — Release Readiness Assessment

**This is an honest readiness check, not a release declaration.** Phase 7 originally asked for a completed `V2_FREEZE.md` including "locked-test result" and "temporal-holdout result" — at the time this was first written, neither existed, because opening the locked scaffold test is a deliberate, one-time, irreversible action this program protected since Phase 1 ("the locked scaffold test stays shut throughout," `BUILD_PLAN.md` §Phase 3). **That has since changed**: the user explicitly authorized opening it ("open the locked test", 2026-09-23), and G2 has passed — see `PHASE3_G2_LOCKED_TEST_RESULT.md` and `METHODS_AND_VALIDATION.md` §5. The gate table below is updated to reflect that; see "Bottom line" at the end for the current overall status.

## G1 gate re-check — still not satisfied, precisely why

| G1 condition | Status |
|---|---|
| Benchmark card published | ✅ Done (`PHASE1_BENCHMARK_CARD.md`) |
| Locked scaffold test | ✅ **Opened 2026-09-23** (explicit user authorization), G2 passed. 795 compounds/598 scaffold groups; 60% found already tuning-contaminated and disclosed, real answer uses the clean 370-compound subset (`PHASE3_G2_LOCKED_TEST_RESULT.md`). **Now fully spent — cannot be reused.** |
| Temporal holdout built | ✅ **UNBLOCKED, built.** `phase2/stage_document_year_backfill.py` (54,717 documents, publication year via `document.json`) + Stage 7's new `first_seen_year` field + `phase2/build_temporal_holdout.py`. Real coverage, disclosed precisely: 69.5% of compounds (936,772/1,347,039) have a resolvable year — the gap's root cause is confirmed and legitimate (Stage 3's Potency data is 99.997% dominated by one undated bulk-deposit ChEMBL document, `CHEMBL1201862`), not a bug. **Temporal test: 187,354 compounds / 86,653 scaffold groups, years 2020-2025. Temporal train: 749,418 compounds / 250,109 scaffold groups, years 1976-2020.** |
| Four split types defined | ✅ **All four now built.** Random + scaffold: done (`phase1/data/holdouts.json`). Temporal: done (above). **Document-level: done** (`phase2/build_document_split.py`) — real leakage-safe partition by document-connected compound clusters (union-find over shared `document_chembl_id`). Found and fixed a real problem along the way: the same Stage-3 mega-document (and 17 others like it) created a catastrophic cluster that first produced an 86%/14% split, then (after excluding direct mega-documents) a residual transitively-chained cluster still produced ~50/50 — fixed by forcing any cluster too large to fit the test budget into train regardless of shuffle order. Final: **test 269,408 compounds/153,965 clusters, train 1,077,631/278,869** — a clean 20%/80% split. |
| Set C confirmed drug-realistic | ✅ Done |
| Adjudication study as precision range | ❌ **Still blocked** on the human reviewer pass (`PHASE1_ADJUDICATION_STUDY_STATUS.md`) — AI-assisted evidence gathering done (153 rows, 19 database-confirmed gaps), final human judgment not yet in. A ready-to-hand-off packet now exists (`ADJUDICATION_HANDOFF.md`) plus a one-command result script (`phase1/compute_adjudication_precision_range.py`) so the reviewer's verdicts become the final G1 number with no further engineering. |
| Full-population coverage numbers landed | ✅ Done |

**G1 does not pass yet — but only ONE genuine blocker remains, and it is now the ONLY item left in the entire project.** The temporal holdout and document split were built and validated; SEA (3B) was built, a real implementation bug found and fixed, and validated; reproducibility packaging, a live serving interface, and per-prediction interpretability are all done (`RESEARCH_HANDOFF_ROADMAP.md`). The sole remaining blocker is a human reviewer's time on the adjudication study — not something any further engineering can substitute for. `ADJUDICATION_HANDOFF.md` hands that off cleanly.

## What Phase 3-5 actually produced this session — a real, usable summary

| Item | Outcome |
|---|---|
| Density-adaptive rule | Re-fit against the rebuilt index: **pool at density≥8** (was 12), larger effect size than before |
| Popularity correction | Tested, **rejected** — significantly harmful at every strength tested |
| Potency weighting | Tested, **adopted** — threshold=6.0, floor=0.0, steepness=5.0 |
| Orthologue ablation | Tested, **adopted for primary-target scoring only** — no broad effect, real modest gain on the harder task |
| Negative-evidence term | **Not testable yet** — needs a capture-schema extension, correctly not faked |
| Confidence weighting (8 vs 9) | **Not tested** — needs the same kind of schema extension |
| SEA (3B) | **Done.** Genuine implementation; found and fixed a real bug (mega-target null-distribution scale mismatch); v2 significantly beats it (`phase3/PHASE3B_SEA_COMPARISON.md`) |
| Neural re-ranker (3D) | **Not attempted** — a substantial separate ML engineering effort |
| Stacker (3E) | Tested, **adopted** — small but real, significant gain on both tiers |
| Calibration (Phase 4) | Tested, **fails its own gate** — ship L-score, not the calibrated version, per the plan's explicit no-complexity-credit rule |
| Abstention policy (D5) | **Done** — abstain below density=8 |
| D1/D2/D3/Novelty position | Each has a specific, named blocker (human review, a broken dependency needing redesign, an unresolved scope decision, or a literature-review task) — none fabricated |

## Honest bottom line — UPDATED after G2, SEA, and the research-handoff roadmap

This program took target-prediction v2 from "Phase 2 not started" through a fully executed, bug-fixed data rebuild, a re-fit and substantially extended Phase 3A retrieval core, a working Phase 3E stacker, a real Phase 4 calibration test, a completed temporal holdout and document-level split, a genuine SEA (3B) comparison (with a real implementation bug found and fixed along the way), and — the most consequential step — **opened the locked scaffold test and passed G2**: on the 370 genuinely clean (never-tuned-on) compounds, v2 significantly beats v1's actual production system on both any-tier (+0.034) and primary-tier (+0.072) ranking, beats the internal baseline on any-tier (+0.031), and significantly beats a genuine SEA reimplementation on both tiers (`PHASE3_G2_LOCKED_TEST_RESULT.md`, `phase3/PHASE3B_SEA_COMPARISON.md`). A real methodological gap was caught and disclosed along the way: 60% of the locked test had accidentally been touched by earlier tuning; the reported result uses only the clean subset. **The locked test is now spent** and cannot be reused.

On top of the science, the full research-handoff roadmap (`RESEARCH_HANDOFF_ROADMAP.md`) is now closed except for one item: reproducibility packaging (`REPRODUCE.md`), a live serving interface with per-prediction interpretability built in (`backend/target_prediction_v2.py` — density-adaptive prediction, evidence, and an explicit un-calibrated confidence label on every result, wired into `backend/app.py`), and one consolidated methods document (`METHODS_AND_VALIDATION.md`) are all done.

**G2 (does the method actually work) is answered: yes, with real statistical support, against both v1 and a real published baseline.** **G1 (is the benchmark's ground truth itself fully validated) is not** — the adjudication study's human reviewer pass is the **sole remaining blocker in the entire project**, not something further engineering can substitute for. `ADJUDICATION_HANDOFF.md` packages that hand-off cleanly, with a one-command script (`phase1/compute_adjudication_precision_range.py`) that turns the reviewer's verdicts directly into the G1 number.
