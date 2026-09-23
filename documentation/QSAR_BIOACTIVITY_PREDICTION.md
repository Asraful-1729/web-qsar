# QSAR Bioactivity Prediction — Complete Technical & Scientific Documentation

**Scope:** the ligand-based potency-prediction feature — the Predict tab (single target), the Compare tab (multi-target analysis), and the Target Info tab (full model transparency/audit browser) — plus the serving pipeline underneath all three (`backend/serving/`, `backend/analysis.py`, `backend/factory_browser.py`). Written directly from the current source code and from the real, shipped model-bucket artifacts (`models/<target_id>/`), not from memory. Docking (a separate, physics-based method) has its own document (`DOCKING.md`); ADMET and Target Fishing are separate features, mentioned here only where they intersect.

**Audience:** anyone who needs to know exactly what a predicted pIC50 number means, where it came from, and where its warranted trust ends.

---

## 1. Objective

Given a target protein with a **pre-trained QSAR model** (built offline, shipped as a "bucket" of files) and a set of candidate small molecules (SMILES), predict each molecule's potency against that target as a **pIC50** value (−log₁₀ of the IC50 in molar units — higher means more potent), rank candidates by it, and — just as importantly — say honestly which predictions should be trusted and which shouldn't. This is a **ligand-based, machine-learned** method: it reasons from statistical structure–activity patterns learned from real measured ChEMBL bioactivity data, in contrast to Docking's physics-based 3D pose simulation. The app treats the two as independent, complementary lines of evidence, never substitutes one for the other, and is explicit throughout that **a predicted number is a prioritization aid, not a validated potency measurement.**

---

## 2. Scientific and methodological basis

### 2.1 Scope boundary, stated up front

This app's backend is **prediction-only — it never trains a model.** The phrase appears literally at the top of `app.py`: *"Prediction-first serving app. NO training here."* Every target's model was produced by a separate, external training pipeline (internally called "the factory," evidenced by `config_used.yaml`/`run_metadata.json`/`chosen_model/metadata.json` inside each shipped bucket) that is **not part of this repository**. This document describes:
- **In full, from the actual code**: how a shipped model is loaded, how a molecule is turned into the exact feature vector the model expects, how a prediction is produced, and how applicability/confidence are computed — all of this *is* in this repo (`backend/serving/`).
- **From the shipped artifacts, not from training code**: what training methodology produced each model — this is reconstructed honestly from the real config and metadata files every bucket ships with (§2.7), since those files are real, on-disk facts, not something outside this repo's visibility.

### 2.2 The feature space

Every molecule is converted to the **same** feature space the factory trained on (`serving/featurize.py`):

1. **Standardization** — RDKit's `rdMolStandardize.Cleanup`, then `LargestFragmentChooser` (keeps only the largest covalently-bonded fragment — strips counterions/salts), then `Uncharger` (neutralizes formal charges where chemically reasonable). A molecule that fails to parse at all is marked unparsable and carried through the pipeline as such, never silently dropped.
2. **RDKit 2D descriptors** — the full `Descriptors.descList` (every descriptor RDKit ships, by name — molecular weight, LogP, topological/electronic descriptors, etc.).
3. **MACCS keys** — 167 structural-key bits (`MACCS_0`…`MACCS_166`).
4. **Morgan/ECFP4 fingerprint** — radius 2, 2048 bits (`Morgan_0`…`Morgan_2047`), the same circular-fingerprint convention used elsewhere in this app (target fishing, decoy generation).

**Per-descriptor fault isolation**: if any single RDKit descriptor throws (a known, occasional RDKit behavior on unusual structures), only that one value falls back to `0.0` — one bad descriptor never fails the whole molecule. This mirrors the same robustness the training pipeline used, so a molecule doesn't get differently-shaped features at prediction time than it would have at training time.

**Why this matters more than almost anything else in the pipeline** (stated directly in the module's own docstring): *"This is the single most important correctness boundary in the app: a naming/order mismatch here silently corrupts every prediction."* A `self_check()` runs once per target at load time (§2.6) specifically to catch this class of bug loudly instead of shipping quietly-wrong numbers.

### 2.3 Two-stage inference — Chemprop, then AutoGluon

This is **not** a single model. Every shipped target uses a two-stage pipeline (confirmed from the actual bucket contents and `model_adapter.py`, not just a high-level description):

1. **Stage 1 — Chemprop D-MPNN.** A directed message-passing graph neural network (the [Chemprop](https://github.com/chemprop/chemprop) package) runs directly on the molecular graph (built from the standardized SMILES, no hand-engineered features) and produces one scalar prediction, `chemprop_pred`.
2. **Stage 2 — AutoGluon stacked ensemble.** `chemprop_pred` is appended as one more column to the RDKit-descriptor/MACCS/Morgan feature vector from §2.2, and the **full** vector is fed to an [AutoGluon](https://auto.gluon.ai/) `TabularPredictor` — a stacked/bagged ensemble of classical tabular models (the shipped example for CHEMBL203_EGFR includes CatBoost, LightGBM/LightGBMXT, RandomForest, ExtraTrees, each at multiple bagging/stacking levels, combined by AutoGluon's own `WeightedEnsemble`). This ensemble's output is the **final predicted pIC50**.

In other words: the graph neural network's own opinion becomes one input feature among thousands to a second, classical-ML ensemble model, rather than being used on its own. Both stages run purely in inference mode on pre-trained, frozen artifacts — nothing is fit or updated at request time.

**Every target's `selected_features.csv` must end with a `chemprop_pred` column** — `Target.__init__` raises a hard error at load time if it doesn't, refusing to serve a bucket that doesn't match this assumed two-stage format rather than silently mis-featurizing it.

### 2.4 Applicability domain (AD) — is this molecule "in scope" for this model?

A QSAR model's predictions are only trustworthy for molecules chemically similar to what it was trained on. This is checked, not assumed (`serving/applicability.py`):

- At load time, the mean and standard deviation of every feature column are computed from the target's own **training** feature matrix (`Data/fit.csv` — the compounds actually used to fit the shipped model, not the full curated dataset).
- For a new molecule, a per-feature z-score is computed, and the **mean absolute z-score across all features** is the molecule's domain distance.
- **Threshold: mean |z| ≤ 3.0** (`DEFAULT_Z_THRESHOLD`) is "in-domain." Above that, the molecule is flagged out-of-domain.

**Where this gate is actually enforced — an important, precise distinction:** the underlying `Target.predict_smiles()` (in `model_adapter.py`) computes a numeric prediction for *every* molecule regardless of domain status — the raw pIC50 value exists in that internal DataFrame either way. **The gate is enforced one layer up**, in `app.py::_rows()`: *"Out-of-domain rows NEVER carry a potency number — AD gating is enforced here, not just in the UI."* A row that is out-of-domain has `predicted_pIC50` set to `None` before it is ever serialized to JSON and sent to the client. The same rule is applied independently in `analysis.py` for the Compare tab's multi-target matrix. This means: **the number never reaches the user for an out-of-domain molecule**, by construction at the API boundary — not merely hidden by a UI convention that a client could bypass.

### 2.5 Confidence tiers — an honest, deliberately modest claim

`serving/confidence.py` opens with a direct statement worth reproducing verbatim, because it explains a choice that could easily look like a missing feature rather than a deliberate one:

> *"The shipped models are AutoGluon TabularPredictors with a single fit/test split — there is no per-compound calibration set, so no genuine per-compound conformal interval exists for them. Rather than fabricate one, the tier below is out-of-domain status + the target's own held-out TEST RMSE, labelled honestly as a test-set error, not a prediction interval."*

The tiering rule:

| Condition | Tier | Label |
|---|---|---|
| Out of applicability domain | **out** | "Outside training chemistry" |
| In domain, target's test RMSE unknown | **med** | "Medium confidence (target test error unknown)" |
| In domain, test RMSE ≤ 0.5 pIC50 | **high** | "High confidence (held-out test RMSE X.XX pIC50)" |
| In domain, 0.5 < test RMSE ≤ 1.0 | **med** | "Medium confidence (held-out test RMSE X.XX pIC50)" |
| In domain, test RMSE > 1.0 | **low** | "Low confidence (held-out test RMSE X.XX pIC50 — wide expected error)" |

This tier is a property of **the target's own model quality** (its held-out test error), combined with **this specific molecule's** domain status — not a per-prediction statistical interval. The label always states which of the two is being reported ("test_rmse" vs. "applicability_domain" as the `confidence_basis`), so a user never has to guess what kind of confidence claim they're looking at.

### 2.6 The startup self-check — a real, hard-failing guard against silent drift

Every time a target bucket is loaded (`Target.__init__`), two independent checks run before any prediction is ever served:

1. **`F.self_check(feature_columns)`** — featurizes a known molecule (aspirin) and asserts every column name the bucket's `selected_features.csv` expects (other than `chemprop_pred`) is actually produced by the current featurizer. If the featurizer's descriptor set has drifted from what this model was trained on (an RDKit version change adding/removing/renaming a descriptor is the realistic failure mode), this raises immediately with the exact missing column names, rather than silently filling them with zeros and serving confidently-wrong predictions.
2. **A feature-parity check against AutoGluon itself** — `self.predictor.features()` (what AutoGluon's loaded model actually asks for at prediction time) must be a *subset* of what `selected_features.csv` declares. AutoGluon legitimately wanting *fewer* columns is harmless (its own internal preprocessing may drop some); AutoGluon wanting a column nobody told the featurizer to produce is treated as **possible corruption**, and the bucket is refused rather than served.

Both checks exist specifically because a naming/version mismatch here is invisible in normal operation — the pipeline would run, produce a number, and that number would simply be wrong, with nothing in the response to indicate it. This is the app's primary defense against that failure mode.

### 2.7 What the shipped training artifacts reveal about model-building methodology

Every target bucket carries real, non-fabricated evidence of how it was built (read directly from `config_used.yaml` / `run_metadata.json` / `chosen_model/metadata.json` — genuine files this app ships, not a description of intent):

- **Training framework**: AutoGluon 1.5.0 `TabularPredictor`, `preset: best_quality`, `time_limit: 7200` seconds, **10-fold bagging**, **2 stacking levels** — a substantial, resource-intensive AutoML search across many base learner families (observed candidate pool for one real target: CatBoost, LightGBM/LightGBMXT, RandomForest, ExtraTrees, each at multiple bag/stack levels, combined into `WeightedEnsemble_L2/L3/L4`).
- **Chemprop D-MPNN settings**: 5-fold CV, up to 200 epochs with 25-epoch early-stopping patience, hidden dimension 300, depth 3, dropout 0.1.
- **Data split**: a fixed random state (42), 15% test fraction; one real target's counts: 7,433 total compounds → 6,480 fit / 953 test, spanning 2,720 distinct Murcko scaffolds (i.e., a real, scaffold-diverse compound set, not a narrow analog series).
- **Activity label**: `pIC50`, consistent with the rest of this app's convention.
- **Model-quality gates actually applied, with real pass/fail results recorded per target**:
  - **Tropsha's applicability-domain/QSAR-validation criteria** — three conditions checked and recorded individually (`R²>0.6`; `(R²−R₀²)/R²<0.1`; `0.85≤k≤1.15`), plus an overall `Tropsha_Pass` boolean.
  - **Y-randomization** — the target label is shuffled and the model re-evaluated; `Y_Random_DeltaR2` records how much worse the shuffled-label model performs, which is the standard check that a model is learning real structure–activity signal rather than overfitting noise. (One real target: `Y_Random_DeltaR2 = 0.8785` — a large drop, consistent with genuine signal.)
  - **Applicability-domain coverage** of the test set itself (`AD_Coverage_pct`).
- **Headline metrics recorded per target, all real, all traceable to `run_metadata.json`**: `R2_Test`, `Q2_Test`, `RMSE_Test`, `MAE_Test`, `Pearson_r`, `Bias`, `SDEP`, `Chemprop_Alone_R2` (the graph network's own standalone performance, for comparison against the full two-stage ensemble), plus `N_Features`, fit/test/total compound counts, and total training runtime.

**Nine standard diagnostic plots** ship per target (`Plots/`), each independently downloadable and annotated (`factory_browser.py::FILE_ANNOTATIONS`):

| # | Plot | What it shows |
|---|---|---|
| 01 | Actual vs Predicted | Test-set scatter — the core calibration check |
| 02 | Residuals vs Predicted | Whether error grows/shrinks with predicted potency |
| 03 | Residual Distribution | Shape of the error distribution |
| 04 | Residual Q-Q | Normality check on residuals |
| 05 | Model Comparison | Cross-validated performance of every base model AutoGluon tried |
| 06 | Y-Randomization | The shuffled-label control described above |
| 07 | Feature Importance | Top contributing features |
| 08 | Applicability Domain | AD coverage of the test set |
| 09 | Target Distribution | Distribution of the real training activity values |

### 2.8 CPU-only, by design

Both inference stages are forced to run on CPU regardless of what hardware is available (`accelerator="cpu"` for the Chemprop/Lightning trainer; thread counts for OMP/MKL/OpenBLAS explicitly capped). The app is designed to never require or depend on a GPU being present at serving time, even though the original training run (per `chosen_model/metadata.json`) used one (`cuda_available: true`, an RTX 4090) — training and serving have deliberately different hardware requirements.

---

## 3. Architecture — module map

### Backend

| Module | Responsibility |
|---|---|
| `serving/model_adapter.py` | **The only module that knows the shipped bucket format.** Loads a target bucket (AutoGluon predictor + Chemprop checkpoint + feature list + AD parameters), runs the self-checks (§2.6), and exposes `Target.predict_smiles()` — the single entry point everything else in the app uses. Bounded LRU cache (default 2 targets resident at once — buckets are large) so repeated predictions against the same target don't reload from disk every time, while memory stays bounded across a session that touches many targets. |
| `serving/featurize.py` | SMILES → standardized SMILES → the full descriptor/MACCS/Morgan feature space (§2.2). Owns `self_check()`. |
| `serving/applicability.py` | Applicability-domain mean-|z| computation (§2.4). |
| `serving/confidence.py` | Confidence tiering from AD status + test RMSE (§2.5). |
| `analysis.py` | Multi-target analysis for the Compare tab: potency matrix, selectivity, polypharmacology, consensus ranking, coverage, per-compound ADMET — everything downstream of calling `predict_smiles()` once per selected target. |
| `factory_browser.py` | Read-only browser/downloader for a target's entire model bucket — metrics, the 9 standard plots (annotated), the full curated dataset, the raw AutoGluon/Chemprop artifacts, or the whole bucket as one ZIP. This is the deliberate "make the rigor visible and auditable" feature — nothing about a model's construction is hidden from a user who wants to check it. |
| `downloads.py` | On-demand fetch of a target's model bucket (and/or docking data) from remote storage, since the installed app does not ship all ~101GB of model buckets locally — a target becomes usable the moment its bucket finishes downloading. Shared infrastructure, not QSAR-specific logic. |

### Frontend

| File | Responsibility |
|---|---|
| `tabs/PredictTab.tsx` | Single-target prediction: target picker, molecule input, ranked results table, CSV export. |
| `tabs/CompareTab.tsx` | Multi-target analysis: target checklist, molecule input, six analytical views over `analysis.py`'s output. |
| `tabs/TargetInfoTab.tsx` | Per-target metrics dashboard + the 9 QSAR plots + bucket file browser, backed by `factory_browser.py`. |
| `components/TargetPicker.tsx` (`PlainTargetSelect`) | Single-target dropdown with download-gating, shared by Predict and Target Info (distinct from `TargetBrowser.tsx`'s disease-first browsing used by Docking/Screen). |
| `components/MoleculeInputPanel.tsx` + `lib/useMoleculeInput.ts` | The shared three-mode molecule input (paste / CSV / SDF) used across Predict, Compare, and other tabs. |
| `components/DownloadGateBar.tsx` | Shared download-progress UI (same component Docking uses). |

---

## 4. The user journey

### 4.1 Predict tab (single target)

**Target** — a plain dropdown of every target whose bucket is complete enough to list (`MA.list_target_ids()`: has `chosen_model/`, `selected_features.csv`, and `Data/fit.csv`). Picking a target not yet downloaded to this machine triggers the same download-gate progress bar used elsewhere in the app; the first ready target is auto-selected on load if none is picked yet. Below the dropdown: real numbers for the currently-selected target (compound count, test R², test RMSE) — not just its name.

**Input** — three modes (`MoleculeInputPanel`, shared everywhere):
- **Paste**: one SMILES per line, split/trimmed client-side.
- **CSV file**: parsed **entirely client-side** — looks for a `smiles`/`SMILES` header (case-insensitive) and reads that column; falls back to the first column if no such header exists. No server round-trip for CSV.
- **SDF file**: uploaded to `/api/parse_sdf`, parsed server-side with RDKit (deliberately not client-side — RDKit's bond/stereochemistry perception on real SDF files is far more reliable than any JS-side parser), returning a clean SMILES list.

**"Rank compounds"** submits synchronously (no job/polling — this is fast enough to run inline) to `/api/predict` (or `/api/predict_csv` for the CSV path).

**Results:**
- A header strip: target name, the actual best model AutoGluon selected (e.g. `LightGBMXT_BAG_L2`), test R² and test RMSE (each with an explanatory tooltip on hover), Tropsha pass/fail, and the in-domain/submitted count.
- A disclaimer line, always shown, verbatim from the backend (`DISCLAIMER = "Prioritisation aid, not a substitute for assays. Trust predictions only for in-domain molecules; treat the top of the list as a shortlist."`).
- **Only in-domain, successfully-parsed compounds appear in the ranked table** — ranked descending by predicted pIC50, with rank number, compound (SMILES, truncated with full text on hover), predicted pIC50, AD z-score, and a confidence chip (colored dot + label, hover for the full explanation).
- Below the table: plain-text counts of out-of-domain and skipped (unparsable) compounds — the current UI does not show these as a browsable list, only as a count, consistent with the "never given a trusted potency number" principle (there is nothing trustworthy to show per-row for them beyond the fact that they exist).
- **Download CSV** exports the ranked, in-domain table only (rank, SMILES, predicted pIC50, AD z, confidence label).

### 4.2 Compare tab (multi-target)

**Pick targets** — a scrollable checklist of every downloaded target; multiple may be checked.

**Input** — the same three-mode `MoleculeInputPanel`.

**"Compare"** calls `/api/predict_multi`, which runs the *full* single-target pipeline independently once per selected target (`analysis.predict_matrix`), then computes several derived views over the combined result (`analysis.analyse`, §2's module-level docstring is worth restating: potency matrix, selectivity, polypharmacology, consensus rank, coverage, ADMET).

**Six views**, switched by a segmented toggle, all reading from the one response:

- **Matrix** — compound × target grid. Each in-domain cell is a heat-colored pIC50 value (color scales between pIC50 4.5 and 8.0); an out-of-domain or unparsable cell shows a plain italic dash/value, never heat-colored (visually distinct from real evidence at a glance). A trailing "In-dom" column shows each compound's in-domain-target-count out of the total selected.
- **Ranking** — the consensus ranking: compounds sorted first by how many targets they're confidently active on (`n_active`, pIC50 ≥ 6.0 by default — `ACTIVE_CUT`), then by mean in-domain predicted pIC50. Coverage (fraction of selected targets the compound was in-domain for) is shown alongside every row, so a high rank achieved on only one target's evidence is never disguised as a rank achieved with broad support.
- **Selective** — compounds with a ≥1.0 log-unit gap (`SELECTIVE_GAP`) between their best in-domain target and the next-best — i.e., a real, quantified selectivity signal, not just "high on one target."
- **Multi-target** — compounds "active" (pIC50 ≥ 6.0) on **two or more** targets simultaneously — the polypharmacology view, useful for diseases where hitting multiple targets is desirable (the module docstring's own example: Alzheimer's AChE + BChE).
- **Best/target** — for each selected target independently, its top 5 in-domain compounds by predicted pIC50.
- **ADMET** — a compact per-compound drug-likeness table (MW, LogP, QED, Lipinski-violation badge, structural-alert count) computed by the app's deterministic ADMET layer (`admet.admet_profile`, a separate feature, called here for convenience so a Compare-tab user doesn't need to re-run the ADMET tab separately).

A single disclaimer (from `analysis.py`, more specific than the Predict tab's) is always shown, stating both the selectivity-gap and active-pIC50 thresholds used and reiterating that out-of-domain predictions are shown but never counted as evidence.

### 4.3 Target Info tab (transparency / audit browser)

**Target** — the same plain dropdown as Predict.

Once picked, shows, all pulled live from the real bucket on disk (nothing here is a cached summary computed once and stored — it is read directly from the same files a maintainer would open by hand):

- A header strip: best model, Tropsha pass/fail, and whether a docking setup also exists for this target (cross-referencing `docking_registry.json` via `/api/docking/status`).
- A metrics grid: test R², test RMSE, Pearson r, AD coverage %, Y-randomization ΔR², and fit/test compound counts.
- **All 9 standard QSAR plots**, rendered as images directly (each fetched from `/api/factory/download/{target}?path=Plots/...`), captioned with their real annotation from `factory_browser.py`'s table (§2.7) — not a generic "plot" label.
- A link to the full cleaned training dataset (`Data/full_cleaned.csv`) with its real file size, downloadable directly.
- (Not currently surfaced in this tab's UI, but available via the same backend: **every other file in the bucket** — the raw AutoGluon predictor/learner pickles, `run_metadata.json`, `selected_features.csv`, the AutoGluon version file, training-environment metadata — each individually downloadable via `/api/factory/download/{target}?path=...`, or the entire bucket at once via `/api/factory/download_all/{target}` as a single ZIP.)

This tab exists specifically so that a model's real-world quality evidence is never just a marketing claim — every number and plot shown anywhere else in the app (Predict tab's header stats, Compare tab's target list) traces back to files a user can open, inspect, and re-derive independently from this same tab.

---

## 5. Background automatic work and decisions — consolidated

### 5.1 Always automatic, no user action

- SMILES standardization (cleanup, largest-fragment, uncharge) before featurization — the user never sees or controls this step; it happens identically for every prediction so results are consistent with how the model was trained.
- The full descriptor/MACCS/Morgan feature computation.
- The Chemprop forward pass and its inclusion as an AutoGluon input feature (the user never sees `chemprop_pred` as a separate number — it is folded into the single final prediction).
- Applicability-domain scoring for every molecule.
- Confidence tiering for every in-domain molecule.
- Suppression of the predicted number for any out-of-domain molecule, enforced at the API layer (§2.4).
- Per-descriptor fault isolation (one bad RDKit descriptor never fails the whole molecule).
- The startup self-check (§2.6) on every target load — invisible unless it fails, in which case the target refuses to load rather than silently mis-predicting.
- LRU eviction of loaded target models once more than `PHYTO_MODEL_CACHE_SIZE` (default 2) targets have been used in a session — a memory-management decision, invisible to the user beyond a possibly-slightly-slower first prediction against a target that had to be reloaded.
- Numeric safety clipping (`±1e6`, NaN/inf → 0) before handing features to AutoGluon's underlying sklearn-family models, which would otherwise raise on an extreme out-of-range value from an unusual molecule rather than simply predicting (badly, but not crashing) — the intent is that an unusual molecule comes back correctly flagged out-of-domain, never as a server error.

### 5.2 Automatic by default, not currently user-adjustable in the UI

- The applicability-domain z-score threshold (fixed at 3.0, `AD.DEFAULT_Z_THRESHOLD`) — a code constant, not exposed as a setting.
- The confidence-tier RMSE cutoffs (0.5 / 1.0 pIC50, `HIGH_RMSE`/`MED_RMSE` in `confidence.py`).
- Compare tab's "active" threshold (pIC50 ≥ 6.0) and selectivity gap (≥1.0 log unit) — real, stated constants (`analysis.py::ACTIVE_CUT`/`SELECTIVE_GAP`), shown to the user in the disclaimer text but not adjustable from the UI.

### 5.3 Never automatic — always an explicit user choice

- Which target(s) to predict/compare against.
- Which molecules to submit, and by which input method.
- Downloading a target's model bucket in the first place (the download-gate is opt-in the moment a not-yet-downloaded target is picked, but always requires that explicit pick).
- Downloading any individual bucket file, plot, or the whole-bucket ZIP from Target Info.

---

## 6. Limitations — stated directly

- **A predicted pIC50 is a point estimate from a machine-learned model trained on historical ChEMBL data, not a measured value.** It reflects patterns in existing bioactivity data, which can be sparse, noisy, or biased toward certain chemical series for reasons unrelated to the target biology.
- **No genuine per-compound confidence interval exists.** The confidence tier is a coarse combination of (a) whether a molecule resembles the training chemistry at all, and (b) the model's own aggregate test-set error — not a calibrated, per-prediction statistical bound. This is stated as a deliberate, honest limitation in the code itself, not something this document is adding after the fact.
- **Out-of-domain does not mean "wrong" — it means "the model has no basis to claim accuracy here."** A molecule far from the training chemistry might still be correctly predicted by chance; the point of the gate is that this can't be claimed with any confidence, so the number is withheld rather than presented as if it could be trusted equally.
- **The training pipeline itself is external to this repository.** This document describes it faithfully from the real artifacts every model bucket ships with (config files, metrics, metadata), but the actual training code (feature selection process, exact AutoGluon/Chemprop invocation, data-cleaning pipeline that produced `Data/full_cleaned.csv`) cannot be audited from within this codebase — only its documented outputs can.
- **Two-stage inference means a Chemprop failure silently zeros that one feature, not the whole prediction.** If Chemprop cannot process a given (already RDKit-parseable) standardized SMILES, `chemprop_pred` for that row defaults to `0.0` rather than blocking the rest of the pipeline — the molecule still gets a prediction from the AutoGluon ensemble using its other ~thousands of features, but with one designed-to-be-informative feature silently missing. This is not currently surfaced to the user as a distinct warning state.
- **The applicability domain check is a single aggregate statistic (mean |z| across all features), not a rigorous domain-boundary method** (e.g. not a convex-hull or leverage-based approach). It is a computationally cheap, real, and genuinely informative proxy, but a coarser one than some published AD methods.
- **Compare tab's "active"/"selective" thresholds are fixed, sensible defaults, not validated cutoffs for every disease/target context** — they are disclosed explicitly in the UI's own disclaimer text so a user applies their own judgment about whether pIC50 ≥ 6.0 (1 µM) or a 1-log selectivity gap is the right bar for their specific question.
- **Model quality varies substantially by target**, exactly as `Test R²`/`Test RMSE`/`Tropsha_Pass` are designed to reveal per target — the app deliberately surfaces this per-target variation (Target Info tab) rather than implying uniform reliability across every shipped model.

---

## 7. API reference (QSAR-prediction-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `GET /api/targets` | List every loadable target with its headline metrics (best model, compound count, test R²/RMSE, AD coverage, Tropsha pass). |
| `GET /api/health` | How many target buckets are present, whether ADMET-AI and Docking are available. |
| `GET /api/health/ml_backends` | Whether AutoGluon's native base-learner libraries (xgboost/lightgbm/catboost) actually loaded their compute engine, not just their Python module — a frozen-build-specific smoke check. |
| `POST /api/predict` | Predict + rank a SMILES list against one target. |
| `POST /api/predict_csv` | Same, from an uploaded CSV. |
| `POST /api/parse_sdf` | SDF → SMILES list (server-side RDKit parsing), used by the SDF input mode across multiple tabs. |
| `POST /api/predict_multi` | Multi-target analysis (Compare tab) — potency matrix, selectivity, polypharmacology, consensus ranking, best-per-target, ADMET. |
| `GET /api/factory/targets` | Every target bucket present, with headline metrics and whether it has a real trained model file. |
| `GET /api/factory/target/{id}/metrics` | Just the metrics for one target. |
| `GET /api/factory/bucket/{id}` | Full file listing for one target's bucket, annotated. |
| `GET /api/factory/download/{id}?path=` | Download one specific file from a target's bucket. |
| `GET /api/factory/download_all/{id}` | The entire bucket (minus the internal `_cache/`) as one ZIP. |

---

## 8. Document provenance

Written by reading, in full: `backend/serving/model_adapter.py`, `featurize.py`, `applicability.py`, `confidence.py`, `backend/analysis.py`, `backend/factory_browser.py`, the relevant sections of `backend/app.py` (`/api/predict*`, `/api/parse_sdf`, `/api/targets`, `/api/health*`), the beginning of `backend/downloads.py`, and every frontend file in the Predict/Compare/Target-Info path (`PredictTab.tsx`, `CompareTab.tsx`, `TargetInfoTab.tsx`, `TargetPicker.tsx`, `MoleculeInputPanel.tsx`, `useMoleculeInput.ts`). The training-methodology claims in §2.7 were read directly from real shipped artifacts for an actual target bucket (`models/CHEMBL203_EGFR/config_used.yaml`, `run_metadata.json`, `chosen_model/metadata.json`) — genuine on-disk files, not a description of intended design. No content here was reconstructed from memory of past conversation.
