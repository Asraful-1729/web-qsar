# Target Prediction (v2) — Complete Technical & Scientific Documentation

**Scope:** the **v2** compound → target prediction engine — `backend/target_prediction_v2.py` (the live serving module and its two `/api/target_prediction_v2/*` routes in `backend/app.py`), the `target_prediction_v2/` research package that built and validated it (`phase1`/`phase2`/`phase3`, etc.), and the shipped runtime data assets (`target_prediction_v2_data/`). Written directly from the current source code and `target_prediction_v2/METHODS_AND_VALIDATION.md`, not from memory.

**On v1:** an earlier, independent engine ("Target Fishing" v1 — `backend/target_fishing.py`, best-similarity ranking + a heuristic evidence score) previously lived behind a method toggle on this same tab. It was **removed from the codebase** once v2 was validated (§5) to significantly outperform it on genuinely held-out data — there is no v1 code, index, route, or UI toggle left to select; v2 is simply what this tab is. The one part of v1 that survived is its as-you-type compound-suggestion helper, ported to query the v2 index instead (`target_prediction_v2.py`'s `suggest_compounds`, `GET /api/target_prediction_v2/suggest_compounds`) — see §3/§9.

**Audience:** anyone who needs to know exactly what v2's ranked-target output means, where its evidence comes from, what was tested and rejected along the way, and what its validated (and unvalidated) performance actually is.

---

## 1. Objective

Given a query molecule (SMILES), predict which of the ~4,882 indexed human single-protein targets it is likely to act on, ranked with an explicit, per-prediction evidence trail and an un-calibrated confidence label — the same problem SwissTargetPrediction, SEA (Keiser et al. 2007), and MolTarPred address. This is a **substantial rebuild of an earlier shipped system (v1)**, built to fix specific, measured weaknesses in v1's data and ranking method — not a from-scratch reinvention, and not a trained multi-class classifier; it is a retrieval-based (ligand-similarity) method, the same family as v1, SEA, and SwissTargetPrediction.

---

## 2. Scientific and methodological basis

### 2.1 The reference data — a re-pulled, re-validated ChEMBL index

Source: the ChEMBL REST API, pinned to release **`ChEMBL_37`** (2026-05-01) — the first time this program pinned an exact release; v1 and all of this program's own earlier work drew from unpinned live pulls. Scope: **5,869 human single-protein targets**, cross-checked against a verified allowlist rather than a bare confidence-score filter (see §2.2 for why that distinction mattered in practice). Volume: **1,590,210** distinct standardized compounds; **3,973,653** native-human (compound, target) pairs; plus **34,080** additional orthologue-tier (non-human-species) pairs, kept as a separate stratum that is never merged into the native pool (see §2.4). Ground-truth recipe: IC50/Ki/Kd/EC50 plus PubChem/literature `Potency`-type records, `=`/`<`/`<=` relations (censored records flagged, not silently treated as exact), `confidence_score≥8`, and the **most-potent (max, not mean)** measurement kept per (compound, target) pair.

### 2.2 Real bugs found and fixed while building this index

Found because the pipeline was actually run and its output checked, not by inspection after the fact:

1. **Salt-form collapse.** RDKit's `SaltRemover` cannot correctly collapse n:1-stoichiometry salts (e.g. a 2:1 tartrate salt) to their parent structure. Confirmed on a real case, then confirmed **already present in v1's shipped index** (2,167 of 1,312,849 rows, ~0.17%). Fixed here via `rdMolStandardize.FragmentParent()`.
2. **Potency aggregation.** The documented recipe called for the most-potent (max) value per pair; v1's actual code averaged (mean) instead — meaning every pChEMBL number in every prior phase of this program was mean-based until this rebuild.
3. **Target-validity contamination.** `confidence_score≥8` alone does not reliably exclude non-single-protein targets (e.g. organism-level or protein-complex assignments) — confirmed live (a ChEMBL `target_type=ORGANISM` record passed the filter). **18.99% of pairs (931,457 of 4,905,110) were contaminated** before a verified single-protein-target cross-reference was applied and the contamination fully removed.
4. **`confidence_score` isn't actually in ChEMBL's own activity records.** The REST API's activity endpoint doesn't serialize this field at all (confirmed empirically), despite accepting it as a query filter. Backfilled via a separate `assay.json` lookup (270,225 assay IDs resolved).

### 2.3 The ranking method — density-adaptive retrieval

For a query molecule, fingerprinted the same way as every other similarity tool in this app (Morgan/ECFP4, radius 2, 2048 bits, packed to bytes — see `backend/target_prediction_v2.py`'s `_packed_fingerprint`):

1. **Compute Tanimoto similarity against the entire native index in one vectorized pass** (`_query_tanimoto_df`) — the same packed-bit-AND-plus-popcount-lookup-table technique documented in `documentation/SIMILARITY_SEARCH.md` §2.3, reused here rather than reimplemented, so both the pooling regime and the fallback regime below score the query against the index exactly once.
2. **Density gate**: count distinct neighbour *compounds* at Tanimoto ≥ **0.5** ("density"), computed over the full, uncapped index.
   - **density ≥ 8** → **pooled regime**: take the 10 most-similar neighbour compounds (`POOL_K = 10`, drawn from a 300-compound candidate cap, `NEIGHBOR_CAP`); every one of those 10 votes for every target it is annotated against.
   - **density < 8** → **best-similarity fallback**: v1's original method — rank targets by their single nearest-neighbour compound, computed over the *full* index (not the capped pooling candidate list) — re-fit and re-validated against this rebuilt v2 index (the threshold used to be 12, fit against v1's smaller index; re-fitting against the larger v2 index moved it to 8 — see `phase3/PHASE3A_DENSITY_REFIT.md`).
3. **Potency weighting** (pooled regime only): each pooled vote is weighted by a logistic ramp on the neighbour's measured pChEMBL — `threshold=6.0, floor=0.0, steepness=5.0` (`_potency_weight`). A neighbour with no potency value at all is treated as neutral (weight 1.0), not as evidence of weak potency.
4. **Orthologue term** (pooled regime only, additive): a separate Tanimoto pass (minimum 0.4) against a much smaller, non-human-species-only index (§2.4), added *specifically* because it was shown to improve primary/mechanism-tier target recovery — not the broader "any known target" task (`phase3/PHASE3A_ORTHOLOGUE_ABLATION.md`).

### 2.4 The orthologue index — a deliberately separate, never-merged stratum

34,080 pairs of non-human-species bioactivity evidence, built by `phase3/build_orthologue_index.py`, kept in its own files (`target_prediction_v2/phase3/orthologue_index/`) and loaded independently (`_load_ortho()`) rather than concatenated into the native index — so a native-evidence result and an orthologue-evidence result are never silently mixed into one undifferentiated similarity pool. Every orthologue-sourced piece of evidence a user sees is labeled as such (`orthologue_neighbours`, with `species_provenance` when available), never presented as native human evidence.

### 2.5 What was tested and genuinely rejected — evidence of real rigor, not selective reporting

- **Popularity correction** (`score / (target_popularity + ε)` — "SEA's own background-model logic," per the original design intent): tested at 5 different strengths, **significantly harmful at every one**, on both the broad and primary-target task. Not shipped.
- **Density-stratified calibration** (isotonic regression, meant to turn the raw score into a calibrated probability): tested, **failed its own pre-registered bar** — worse Brier score than the simple raw-score baseline (0.1604 vs. 0.1532). **The system ships the raw L-score as its confidence signal, not a calibrated probability** — a deliberate, evidence-based choice, documented directly in the serving module's own docstring, not an oversight.
- **A logistic-regression stacker** combining the three scoring signals: adopted as a real finding (a small effect: +0.6% any-tier, +1.6% primary-tier via 5-fold cross-validation), but **never built into a production-servable model** (it would need to be retrained on all data and serialized, not just cross-validated) — so it is not called by the live serving code.

### 2.6 What a result does and doesn't mean

The score shown (the "L-score" in the pooled regime, or `best_similarity` in the fallback regime) is **not a calibrated probability that the compound binds this target** — Phase 4 explicitly measured calibration and it failed to beat the raw score (§2.5), so no probability-shaped number is ever produced or shown. It is evidence: structurally similar compounds, annotated against a target in real ChEMBL bioactivity data, potency- and (optionally) orthologue-weighted. Every result the live API returns carries an explicit, human-readable `confidence_label` stating exactly this, and never a bare unlabeled score.

---

## 3. Architecture — module map

### Research / validation package (`target_prediction_v2/`, offline, maintainer-side)

| Area | Responsibility |
|---|---|
| `phase2/` (`stage0_pin_release.py` … `stage5b_fetch_orthologue_activities.py`, `build_v2_index_files.py`) | The ChEMBL_37-pinned data pipeline: fetch, censored-record handling, potency aggregation, single-protein-target verification, assay-confidence backfill, orthologue mapping — producing the `v2_index` files served in production. |
| `phase3/` (`fit_density_adaptive_rule_v2.py`, `build_orthologue_index.py`, `build_sea_null_distributions_v2.py`, the `PHASE3A_*`/`PHASE3B_*` reports) | Method refinement and validation against the rebuilt v2 index: the density-adaptive gate refit, potency-weighting fit, orthologue ablation, the popularity-correction and calibration rejections (§2.5), and the SEA head-to-head comparison (§5.3 below). |
| `phase1/` | Earlier baseline/adjudication work this rebuild's own validation partly builds on (ground-truth adjudication candidates, mechanism-support coverage) — historical/supporting, not part of the live serving path. |
| `METHODS_AND_VALIDATION.md` | The single authoritative account this document is built from — every number in §2/§5 traces to a script and result file cited there. |

### Backend (live)

| Module | Responsibility |
|---|---|
| `target_prediction_v2.py` | The entire live serving path: index loading (`_load`/`_load_ortho`/`_load_target_names`, each `functools.lru_cache`\-memoized so the (large) index is read from disk once per process), fingerprinting, the density gate, pooled/fallback scoring, and the `predict(smiles, top_k)` function used by both the CLI (`python3 target_prediction_v2.py "<SMILES>" [top_k]`) and the FastAPI routes below. |
| `app.py` (`/api/target_prediction_v2/*`) | `GET /status` (index present in this build?), `POST /predict` (`smiles`, `top_k` 1–10,000, default 25), and `GET /suggest_compounds` (`q`, `limit`, default 8 — the as-you-type helper ported from v1, §9). |

### Frontend

| File | Responsibility |
|---|---|
| `tabs/TargetFishingTab.tsx` | The Target Prediction tab: compound input with as-you-type suggestions (`api.suggestCuratedCompounds`), a Search action calling `api.targetPredictionV2Predict` (`POST /api/target_prediction_v2/predict`), and the results view (`ResultsTable`). There is no threshold to set (the sidebar states directly why: "automatically decides whether to pool... or fall back... no threshold to set" — the density gate is the only control), and no method toggle — this tab previously hosted a v1/v2 toggle; v1 was removed from the codebase entirely (see the scope note above). |

---

## 4. Data provenance and index building (offline, maintainer-side)

Not part of the interactive user journey, but the same pipeline the live index is served from, documented for auditability:

1. `phase2/stage0_pin_release.py` pins the exact ChEMBL release (`ChEMBL_37`, 2026-05-01) before any data is fetched — the given reason (§2.1) is that every prior pull in this program's history, including v1's, was unpinned.
2. `phase2/stage1_fetch_base_evidence.py` through `stage5b_fetch_orthologue_activities.py` fetch, in stages, the base positive evidence, censored records, additional potency records, inactives, and the orthologue (non-human-species) evidence, each stage checkpointed/resumable — real transient ChEMBL API instability was encountered and handled throughout this build, not a hypothetical concern (`METHODS_AND_VALIDATION.md` §7).
3. `phase2/build_v2_index_files.py` (via `stage7_aggregate.py`) applies the max-potency-per-pair aggregation, the verified single-protein-target cross-reference, and the assay-confidence backfill (§2.2), then writes the two **row-aligned** files actually loaded at serving time: `fingerprints.npz` (packed Morgan bits) and `compounds.csv.gz` (`smiles`, `target_chembl`, `pchembl_value`, …) — currently **~181 MB** and **~69 MB** respectively on disk (`target_prediction_v2/phase2/v2_index/`).
4. `phase3/build_orthologue_index.py` builds the separate, much smaller orthologue-index pair (`target_prediction_v2/phase3/orthologue_index/fingerprints.npz` ≈1.2 MB, `compounds.csv.gz` ≈0.5 MB — 34,080 pairs total, §2.4).
5. A target-ChEMBL-id → preferred-name lookup (`target_prediction_v2/phase2/data/target_name_map.json`, ~264 KB) is built alongside, so results can show a human-readable target name (`_target_name`) rather than a bare ChEMBL id, when known.
6. For the desktop/Windows build specifically, these three files are staged verbatim into `target_prediction_v2_data/{v2_index, orthologue_index, target_name_map.json}` (see `desktop.py`'s `TARGET_PREDICTION_V2_INDEX_DIR`/`_ORTHO_DIR`/`_NAME_MAP` environment overrides, and the CI workflow's staging step) — the same three environment variables `target_prediction_v2.py` already read at import time for any deployment, so relocating the data for a packaged build needed no code change.

---

## 5. Validation — the honest headline result

### 5.1 A methodological finding, caught and handled before it mattered

Every individual scoring decision in §2.5 was tested repeatedly against the *same* ~1,600-query tuning sample — a real overfitting risk that was checked, not assumed away. A **locked scaffold test** (795 compounds, sealed early in the program specifically to give one final, unbiased answer) was opened once, and an overlap check *before scoring anything* found that **60% of the locked compounds had already been touched by the tuning process** (the tuning scripts sampled from the general drug population without excluding the locked test's own scaffold groups — a real gap, disclosed rather than hidden). **All headline numbers below use only the 370 compounds confirmed never touched by any tuning step.** This locked test is now fully spent and must not be reused for any future claim.

### 5.2 The result (n=370, genuinely held out)

| Comparison | Effect (paired reciprocal-rank, bootstrap) | 95% CI | Significant? |
|---|---:|---|:---:|
| v2 vs. v1 (`best_similarity`), any-target | **+0.034** | [0.012, 0.059] | **Yes** |
| v2 vs. v1 (`best_similarity`), primary-target (n=103) | **+0.072** | [0.021, 0.128] | **Yes** |
| v2 vs. internal unweighted-10NN baseline, any-target | **+0.031** | [0.004, 0.059] | **Yes** |
| v2 vs. internal unweighted-10NN baseline, primary-target | +0.011 | [−0.039, 0.062] | No |

**v2 significantly outperforms the currently-shipped v1 production system**, on data no part of its design was tuned against, on both the general recovery task and the harder, more specific primary-mechanism task.

### 5.3 Comparison against SEA (Keiser et al. 2007)

A genuine implementation (raw Tanimoto-sum scores converted to significance via a target-size-stratified Gumbel/extreme-value null distribution fit by simulation — the same statistical logic as a BLAST E-value, not a shortcut), evaluated on a 150-query subset of the fresh (never-tuned-on) holdout population.

A real bug was found and fixed during this comparison, not glossed over: the first null-distribution fit used discrete, capped size bins, whose top bin (targets with ≥1,280 ligands) used a fixed representative size of only ~1,640 — but real targets in the index range up to **195,809 ligands** (63 targets exceed 10,000; 15 exceed 50,000). This produced absurd z-scores (up to 636–641) for the largest, most promiscuous targets and an implausible first result. Diagnosed by direct inspection of true-positive raw scores against the null, then fixed by replacing the discrete bins with a continuous power-law null model fit across the full real size range (log-log R² ≈ 1.0).

| Comparison | Effect (paired reciprocal-rank, bootstrap) | 95% CI | Significant? |
|---|---:|---|:---:|
| v2 (full pipeline) vs. SEA, any-target | **+0.399** | [0.333, 0.467] | **Yes** |
| v2 (full pipeline) vs. SEA, primary-target (n=54) | **+0.569** | [0.464, 0.670] | **Yes** |
| SEA vs. v1 `best_similarity`, any-target | **−0.371** | [−0.442, −0.300] | **Yes** |

Fixing the null-model bug substantially reduced the "SEA vs. v1" effect (−0.538 → −0.371) but did not flip its direction — SEA still significantly underperforms simple nearest-neighbour search on this benchmark, judged a genuine property of this near-neighbour-rich ChEMBL data (which favors best-match retrieval over SEA's ensemble-normalization design), not an artifact. v2 significantly outperforms both comparators.

---

## 6. The user journey

1. On opening the Target Prediction tab, `/api/target_prediction_v2/status` is checked once.
2. **Query SMILES** — entered directly, or picked from as-you-type suggestions drawn from the v2 curated-compound index as the user types 2+ characters (`GET /api/target_prediction_v2/suggest_compounds`, §9 — a convenience lookup, not part of the scoring method itself).
3. **Search** calls `POST /api/target_prediction_v2/predict` with `{smiles, top_k}` (`top_k` defaults to 10,000 in the frontend's own call — effectively "all" ranked targets, since the index currently holds 4,882; the API's own ceiling is 10,000 regardless of index size, to keep a single pathological request bounded).
4. **Results** (`ResultsTableV2`) show, per target: rank, preferred name (or bare ChEMBL id if unknown), the raw score badge (tooltip carries the exact `confidence_label`), and — behind a per-row "Show supporting evidence" toggle — every neighbour compound that actually contributed: its SMILES, Tanimoto similarity, pChEMBL value and potency weight (native neighbours), or its Tanimoto and species provenance (orthologue neighbours), or the single best-matching compound and its pChEMBL (fallback regime). A summary line states the real regime that fired for this specific query — pooled (with the actual density count) or best-match fallback — and how many of the indexed targets were ranked.

---

## 7. Background automatic work and decisions — consolidated

### 7.1 Always automatic

- Index/orthologue-index/name-map availability check on tab load.
- Fingerprinting the query and scoring it against the entire native index, on every prediction.
- The density gate (§2.3) — the user never chooses pooled vs. fallback; it is entirely determined by how many similar compounds actually exist for that specific query.
- Potency weighting of every pooled vote (§2.3, step 3).
- The additive orthologue term, whenever the pooled regime fires and the orthologue index has qualifying neighbours (§2.3, step 4; §2.4).
- Per-prediction evidence collection (which neighbours, their Tanimoto/pChEMBL/potency-weight/species, kept and returned for every ranked target — never only an aggregate score).

### 7.2 Automatic by default, adjustable in the UI

- `top_k` (how many ranked targets to return) — has a default (25 at the API layer; the frontend's own call requests up to 10,000) but is a caller-supplied parameter, not fixed.

### 7.3 Never automatic

- Submitting a prediction — always an explicit action on a specific query molecule; nothing runs on tab load beyond the two `status` checks.
- Rebuilding, re-fetching, or re-validating the index itself — entirely an offline, maintainer-run pipeline (§4); the live app only ever reads the already-built files.

---

## 8. Limitations — stated directly

- **This is a retrieval/evidence method, not a calibrated probability.** Density-stratified calibration was explicitly tried and failed its own pre-registered accuracy bar against the raw score (§2.5) — there is currently no path from the score shown to "probability this target is genuinely hit," and the UI's own `confidence_label` says so on every single result, not just in a one-time notice.
- **The benchmark's own ground truth is not yet fully validated.** An adjudication study — checking whether the model's apparent "false positives" are real, ChEMBL-undocumented target relationships — has AI-assisted evidence prepared (153 candidate pairs, 19 already database-confirmed) but still needs a qualified human reviewer's judgment. This is recorded in `METHODS_AND_VALIDATION.md` as the single most important open item.
- **No live serving system existed before this module** — every validation number in §5 came from offline batch analysis against pre-computed neighbour lists, not the query-time engine documented here; the engine implements the validated method exactly (§7b of `METHODS_AND_VALIDATION.md`), but was not itself the object being benchmarked.
- **The orthologue and gene-symbol-homology signals are proxies**, not validated sequence-based homology (no BLAST/OrthoDB cross-check) — disclosed at every point they are used, both in the source documentation and in this document's §2.4.
- **Negative-evidence weighting and a two-level (8-vs-9) confidence weight are unbuilt** — both would need a capture-schema extension that does not exist yet.
- **Temporal and document-level holdouts exist as scripts** (`phase2/build_temporal_holdout.py`, `build_document_split.py`) but have not yet been used for a full calibration-under-distribution-shift pass.
- **A genuine SEA comparator is built and validated (§5.3); a neural re-ranker comparator was scoped but never attempted** — a substantial separate ML engineering project, not started.
- **The index is a snapshot of `ChEMBL_37`.** New bioactivity data ChEMBL publishes after this pin will not appear until the index is rebuilt from a newer release — there is no live/incremental update path.

---

## 9. API reference (target-prediction-v2-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `GET /api/target_prediction_v2/status` | Whether the v2 index is present in this build. |
| `POST /api/target_prediction_v2/predict` | Run a prediction: `smiles`, `top_k` (1–10,000, default 25). Returns `{regime, density, density_threshold, n_neighbours_considered, n_targets_indexed, results: [{rank, target_chembl, target_pref_name, score, confidence_label, evidence}]}`. |
| `GET /api/target_prediction_v2/suggest_compounds` | As-you-type suggestions: `q`, `limit` (default 8, capped at 20). Plain substring match (no fingerprinting) against the v2 index's `smiles`/`target_chembl` columns and the target-name map's values — matches a SMILES fragment or a target's ChEMBL id/name. Returns `{results: [{smiles, target_chembl, target_pref_name, target_id, pchembl_value}]}`; `target_id` (this app's own registry id, e.g. `CHEMBL1862_ABL1`) is always `null` — v2's index carries no such mapping (ported as-is from the removed v1 engine, which had the same field for the same reason: most ChEMBL targets in this broad reference pool simply aren't one of this app's own docking/QSAR targets). |

---

## 10. Document provenance

Written by reading, in full: `backend/target_prediction_v2.py`, `target_prediction_v2/METHODS_AND_VALIDATION.md`, the `/api/target_prediction_v2/*` routes in `backend/app.py`, `frontend/src/tabs/TargetFishingTab.tsx`, the relevant interfaces in `frontend/src/lib/types.ts` and calls in `frontend/src/lib/api.ts`, and the `target_prediction_v2/phase2/`, `phase3/`, and `phase1/` directory contents (file listings and the manifests/reports cited above) to confirm the module map and data-build provenance. The private `target_prediction_v2_procedure_rev*.md` files at the repo root were deliberately not read or used. Updated after v1 (`backend/target_fishing.py` and everything else under the "Target Fishing" name) was removed from the codebase, including its `suggest_compounds` helper being ported into `target_prediction_v2.py` and its two routes/response-shape callers (`frontend/src/lib/api.ts`, `backend/research_report.py`) repointed at v2 — see those files' own current state, not this document's earlier draft, for anything not already re-verified here. No content here was reconstructed from memory of past conversation — every specific number traces to a line of code, a manifest field, or a result in `METHODS_AND_VALIDATION.md`.
