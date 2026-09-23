# Reproducing Target Prediction v2 — from nothing to a validated, servable model

**Purpose**: a single, ordered, copy-pasteable script sequence for a reader who was not present for this program's development (research-handoff roadmap item #4). Every script referenced here already exists and was actually run to produce every result in `METHODS_AND_VALIDATION.md` — this document does not describe new work, it documents the real dependency order of work already done, so it can be reproduced (or re-run against a fresh ChEMBL release) without archaeology through 50+ phase-specific markdown files.

## 0. Environment

Everything below runs under a dedicated conda environment with a working RDKit build (the system default `python3` has a broken/incompatible numpy2+pandas combination and no RDKit):

```
/home/storage/jonaid/.conda_envs/vsdock/bin/python3
```

It needs: `rdkit`, `numpy`, `pandas`, `scipy` (already present) and `scikit-learn` (installed mid-program for the Phase 4 calibration test and the Phase 3E stacker — `pip install scikit-learn` into vsdock if starting fresh). The live serving module (`backend/target_prediction_v2.py`) needs the same environment (or any Python with rdkit/numpy/pandas) plus `fastapi`/`uvicorn` if serving over HTTP via `backend/app.py`.

All ChEMBL API calls use the pinned release recorded by `stage0_pin_release.py` (`ChEMBL_37`, 2026-05-01) — re-running against a live/different release will reproduce the *method*, not the exact numbers in `METHODS_AND_VALIDATION.md`.

Every staged script below is **checkpointed and resumable** — real ChEMBL API transient failures (HTTP 500s) were hit repeatedly during the original run and every stage handles being re-invoked after a partial failure by skipping already-fetched records. If a stage script exits partway through, just re-run it.

## 1. Phase 2 — data rebuild (`phase2/`)

Run in this exact order — each stage's output is a real input to the next, not just a numbering convention:

| # | Script | Produces | Depends on |
|---|---|---|---|
| 1 | `stage0_pin_release.py` | pinned release manifest, `data/single_protein_target_ids.json` (verified allowlist) | — |
| 2 | `stage1_fetch_base_evidence.py` | `data/stage1_activities.jsonl` | 1 |
| 3 | `stage2_fetch_censored.py` | `data/stage2_censored.jsonl` | 1 |
| 4 | `stage3_fetch_potency.py` | `data/stage3_potency.jsonl` | 1 |
| 5 | `stage4_fetch_inactives.py` | `data/stage4_inactives.jsonl` | 1 |
| 6 | `stage5a_build_orthologue_map.py` | `data/orthologue_target_map.json` | 1 |
| 7 | `stage5b_fetch_orthologue_activities.py` | `data/stage5b_orthologue_activities.jsonl` | 6 |
| 8 | `stage6_standardize.py` | `data/stage6_std_cache.json`, `data/stage6_fingerprints.npz`, `data/stage6_fp_meta.json` (salt-stripped via `rdMolStandardize.FragmentParent()`, not RDKit's `SaltRemover`) | 2-5, 7 |
| 9 | `stage_confidence_backfill.py` | `data/assay_confidence_map.json` (270,225 assay lookups — `confidence_score` is not in `activity.json` itself) | 2-5 |
| 10 | `stage_document_year_backfill.py` | `data/document_year_map.json` (54,717 document lookups — enables the temporal holdout) | 2-5, 7 |
| 11 | `stage7_aggregate.py` | `data/stage7_native_human_pairs.jsonl`, `data/stage7_orthologue_pairs.jsonl` (max, not mean, pChEMBL per pair; single-protein-target cross-referenced against #1's allowlist — fixes an 18.99% contamination bug if skipped) | 1, 8, 9, 10 |
| 12 | `build_v2_index_files.py` | `v2_index/compounds.csv.gz`, `v2_index/fingerprints.npz` (v1-compatible format) | 8, 11 |
| 13 | `build_v2_leakage_fp2.py` | `data/v2_leakage_fp2.npz` → filter to `data/v2_leakage_fp2_filtered.npz` (RDKit topological fingerprints for leakage-controlled evaluation only, not production serving) | 8 |
| 14 | `build_temporal_holdout.py` | `data/temporal_holdout.json` (scaffold-group-safe, year-ordered split) | 11 |
| 15 | `build_document_split.py` | `data/document_split.json` (document-connected compound clusters, union-find) | 11 |
| 16 | `stage_target_name_backfill.py` | `data/target_name_map.json` (4,882 ChEMBL target names — needed for the serving layer's display names, not the ranking math) | 12 |

**Phase 1** (`phase1/build_holdouts.py`, not re-run above) already produced `phase1/data/holdouts.json`, which contains the **locked scaffold test** — 795 compounds sealed at the start of this program. **Do not open it casually**: it is a one-time, deliberately irreversible evaluation, already spent (§3 below) and must not be reused.

## 2. Phase 3 — retrieval core, levers, and comparators (`phase3/`)

| # | Script | Produces | Depends on |
|---|---|---|---|
| 17 | `capture_scaffold_strict_v2.py` (run for Set A and Set B) | `results/capture_scaffold_strict_v2_{A,B}.jsonl` — the tuning sample every lever below was tested against | Phase 2 §1 |
| 18 | `fit_density_adaptive_rule_v2.py` | `PHASE3A_DENSITY_REFIT.md` — the density≥8 pooling threshold | 17 |
| 19 | `test_potency_weighting.py` / `test_potency_weighting_round2.py` | `PHASE3A_POTENCY_WEIGHTING.md` — adopted (threshold=6.0, floor=0.0, steepness=5.0) | 17 |
| 20 | `build_orthologue_index.py` | `orthologue_index/compounds.csv.gz`, `fingerprints.npz` | Phase 2 §1 |
| 21 | `test_orthologue_ablation.py` | `PHASE3A_ORTHOLOGUE_ABLATION.md` — adopted for primary-tier scoring only | 17, 20 |
| 22 | `test_popularity_correction.py` | `PHASE3A_POPULARITY_CORRECTION.md` — **rejected** (harmful at every strength tested); kept only as a disclosed negative result | 17 |
| 23 | `test_stacker.py` | `PHASE3A_STACKER.md` — real small gain, **not** serialized into production (see `backend/target_prediction_v2.py`'s docstring) | 17, 19, 21 |
| 24 | `capture_fresh_holdout.py` | `results/capture_fresh_holdout_{A,B}.jsonl` — a genuinely new sample, explicitly excluding every compound already in #17 | Phase 2 §1 |
| 25 | `test_fresh_holdout_validation.py` | `results/fresh_holdout_validation.json` — first assembled-pipeline check before opening the locked test | 18, 19, 21, 24 |
| 26 | `build_sea_null_distributions_v2.py` | `results/sea_null_distributions_v2.json` — continuous power-law null model (fixes a real bug in the superseded `build_sea_null_distributions.py`/discrete-bin version — see `PHASE3B_SEA_COMPARISON.md`) | Phase 2 §1 |
| 27 | `test_sea_vs_v2_fixed.py` | `results/sea_vs_v2_result_fixed.json` — the trustworthy SEA comparison (`PHASE3B_SEA_COMPARISON.md`) | 24, 26 |

### 3. The locked scaffold test — one-time, already spent

```
capture_locked_test.py   ->  results/capture_locked_test.jsonl
test_locked_test_g2.py   ->  results/locked_test_g2.json  (PHASE3_G2_LOCKED_TEST_RESULT.md)
```

**Do not re-run these against a new sample of the same locked compounds and call it a fresh G2 check** — the locked test (`phase1/data/holdouts.json`) was opened once (2026-09-23, explicit user authorization) specifically because it is meant to be spent exactly once. A 60% tuning-overlap contamination was found and disclosed (`METHODS_AND_VALIDATION.md` §5.1); the real result uses only the 370 genuinely clean compounds. If a truly new held-out answer is ever needed again, it requires a **new** locked set, built the same way `phase1/build_holdouts.py` did, sealed before any further tuning touches it.

## 4. Phase 4 — calibration (`phase4/`)

```
test_calibration.py  ->  PHASE4_CALIBRATION.md
```
Depends on Phase 3 §2 (uses the same tuning captures + 5-fold scaffold-grouped CV). **Result: failed its own pre-registered gate** — ships the raw L-score, not this calibrated version. Re-running this only matters if re-attempting calibration with a different method.

## 5. Phase 6 — pathogen-module scoping (`phase6/`, optional)

```
build_bacterial_target_universe.py    -> data (665 bacterial single-protein targets)
build_bacterial_homology_check.py     -> data (12/665 share a human gene symbol)
```
Independent of everything else — a scoping exercise, not a required dependency for the core human-target model.

## 6. Serving — no build step

`backend/target_prediction_v2.py` reads directly from `phase2/v2_index/`, `phase3/orthologue_index/`, and `phase2/data/target_name_map.json` (steps 12, 16, 20 above) — nothing further to build. See its own module docstring for the exact method, and `METHODS_AND_VALIDATION.md` for what is/isn't validated. CLI: `python3 backend/target_prediction_v2.py "<SMILES>" [top_k]`. HTTP: `POST /api/target_prediction_v2/predict` once `backend/app.py` is running (`uvicorn app:app --app-dir backend --host 127.0.0.1 --port 8000`, per the repo root `README.md`).

## 7. Known gaps in this reproducibility pass (disclosed, not hidden)

- Script names still carry real `_run2`/`_v2`/`_fixed` suffixes from bugs found and fixed live during this build (e.g. `build_sea_null_distributions.py` is superseded by `_v2`, `test_sea_vs_v2.py` by `_fixed`) — the superseded files are kept for transparency, not deleted, but this table is the authoritative "what to actually run" reference.
- No single `run_all.sh` wrapper exists — the table above is accurate and ordered, but each script is still invoked individually. Building a wrapper is a reasonable follow-up if this needs to be re-run unattended, not done here since every stage already needed live inspection of its output when it was originally run (that is how the real bugs in §1 were caught).
- `requirements.txt`/environment pin: `vsdock`'s exact package versions were not exported to a `environment.yml` — re-create by installing `rdkit`, `numpy`, `pandas`, `scipy`, `scikit-learn` (any reasonably recent mutually-compatible set; this program hit real numpy2/pandas incompatibility with the system default Python, which is why a dedicated env was used at all).
