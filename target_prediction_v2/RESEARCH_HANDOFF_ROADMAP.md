# Research Handoff Roadmap — Making v2 High-Quality for Professor Handoff

**Context**: the goal shifted from "ship a product feature" to "produce a defensible research tool for academic use." This changes priorities — statistical rigor and reproducibility now matter more than UI polish. This document tracks the plan; update status inline as items close.

---

## 1. Adjudication study — human review (PRIORITY, longest lead time)

**Status: NOT STARTED — needs the professor specifically.**
The one task requiring genuine domain expertise, not more engineering. `phase1/data/adjudication_pilot_output.csv`, 153 rows, 19 already database-confirmed by AI-assisted evidence gathering (`PHASE1_ADJUDICATION_STUDY_STATUS.md`). **Reframe this as bringing the professor in to help validate the tool, not handing them a finished black box** — get this in front of them as early as possible given the lead time.

## 2. Comparison against a real published method (SEA)

**Status: DONE.**
`BUILD_PLAN.md` 3B was previously marked "not attempted" (disclosed, not faked) — built a genuine SEA implementation (Gumbel/EVD null-distribution fitting, not a shortcut). Found and fixed a real bug along the way: the first null model broke down for mega-promiscuous targets (up to 195,809 ligands), producing implausible z-scores (up to 641) and a false "SEA catastrophically worse than v1" result. Fixed with a continuous power-law null model (log-log R²≈1.0), re-validated via true-positive z-score sanity checks. Final, trustworthy result: v2 significantly beats both SEA and v1 (`phase3/PHASE3B_SEA_COMPARISON.md`); SEA itself still significantly underperforms v1's simple nearest-neighbour search on this benchmark — judged a genuine finding (near-neighbour-rich data favors best-match retrieval over SEA's ensemble design), not an artifact.

## 3. One consolidated methods document

**Status: DONE** (`METHODS_AND_VALIDATION.md`, repo root).
Everything in this program is real and rigorous but scattered across 50+ phase-specific markdown files accumulated over the full build. One authoritative document now exists: method description, what's validated (with the real numbers, including the completed SEA comparison), exact known limitations, citation/provenance — written for a reader who wasn't in this conversation. Should be kept in sync as later roadmap items close.

## 4. Reproducibility packaging

**Status: DONE** (`REPRODUCE.md`, repo root).
The pipeline is still a sequence of individually-invoked scripts (real `_run2`/`_v2`/`_fixed` retry artifacts from bugs hit and fixed live during this build, kept rather than deleted for transparency), but a single ordered, dependency-accurate table now exists — a stranger can reproduce the index, every Phase 3 result, and the G2 result without archaeology through 50+ markdown files. No `run_all.sh` wrapper (disclosed as a known gap in `REPRODUCE.md` §7) — each stage genuinely needed live inspection of its output the first time, which is how the real bugs (salt-collapse, mean-vs-max pChEMBL, 18.99% contamination, the locked-test overlap, the SEA null-model bug) were actually caught.

## 5. Usable interface — CLI/Python-API first

**Status: DONE** (`backend/target_prediction_v2.py`).
A live serving module: `predict(smiles, top_k)` in Python, plus a CLI (`python3 target_prediction_v2.py "<SMILES>" [top_k]`) and two FastAPI routes (`GET /api/target_prediction_v2/status`, `POST /api/target_prediction_v2/predict`) wired into `backend/app.py` alongside (not replacing) the existing v1 `target_fishing` endpoints. Verified end-to-end: aspirin correctly recovers Prostaglandin G/H synthase 1 (COX-1) at rank 2 with real supporting evidence; sparse queries correctly trigger the best-similarity fallback regime; invalid SMILES and out-of-range parameters return proper 400/422 errors (tested via FastAPI's TestClient, both regimes and error paths). Along the way, fixed a disclosed gap in the v2 index itself: `target_pref_name` was never fetched (`build_v2_index_files.py`'s own docstring flagged this) — built `phase2/stage_target_name_backfill.py` (4,887/4,887 ChEMBL target names resolved) so results show real names, not bare ChEMBL ids.

## 6. Per-prediction interpretability

**Status: DONE** — satisfied by item #5's design, not a separate build.
Every `predict()` result includes *why*, not just a rank: supporting neighbour compounds with their Tanimoto similarity and potency weight (`evidence.native_neighbours`), orthologue provenance where used (`evidence.orthologue_neighbours`, including source species), which retrieval regime fired (`"pooled"` vs `"best_similarity_fallback"`) and why (explicit density count vs. the density≥8 threshold), and an explicit un-calibrated confidence label on every single result (never presented as a probability, consistent with Phase 4's calibration test failing its own gate).

---

## Already strong, don't redo
- Every methodological choice has real statistics (scaffold-clustered BCa CIs throughout, disclosed negative results — popularity correction rejected, calibration failed its own gate, rather than hidden).
- Data provenance: pinned ChEMBL release (`ChEMBL_37`, 2026-05-01 — first time this program pinned one), checksummed freeze artifacts (`V2_FREEZE.md`), real bugs found and fixed with disclosure rather than silently patched (salt-collapse, mean-vs-max pChEMBL, 18.99% invalid-target contamination, the locked-test tuning-overlap issue).
- G2 passed on genuinely clean, never-tuned-on data (`PHASE3_G2_LOCKED_TEST_RESULT.md`) — v2 significantly beats v1 on both any-tier and primary-tier ranking.
