# ADMET Profiling — Complete Technical & Scientific Documentation

**Scope:** the ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) profiling feature — the ADMET tab, its two independent prediction layers, and the isolated worker service that runs the heavier of the two. Written directly from the current source code (`backend/admet.py`, `backend/admet_endpoints.py`, `backend/admet_service.py`, the `/api/admet*` routes in `backend/app.py`, `frontend/src/tabs/AdmetTab.tsx`) — not from memory.

**Audience:** anyone who needs to know exactly what an ADMET flag or predicted endpoint means, which of the two layers produced it, and how much to trust it.

---

## 1. Objective

Given a set of candidate molecules, produce a **drug-likeness and liability profile** — physicochemical properties, rule-of-thumb drug-likeness checks, known problematic-substructure alerts, and (when available) machine-learned predictions across 34 real pharmacokinetic and toxicity endpoints — **independent of any specific biological target**. This is deliberately the one prediction feature in the app that is *not* about "does this bind protein X" — it answers "is this molecule likely to behave like a usable drug at all," which matters regardless of which target a researcher is ultimately interested in. Like every other prediction in this app, results here are a triage aid, explicitly never used to silently filter or exclude a compound from consideration.

---

## 2. Scientific and methodological basis

### 2.1 Two independent layers, one always on

**Layer 1 — Deterministic (RDKit), always available, no external dependency.** Pure, fast, local cheminformatics — the same molecule always produces the same output, with no model weights, no network call, and no possibility of being "down."

**Layer 2 — Learned (ADMET-AI), available only if a separate worker process is running.** A real trained multi-task model making genuine per-endpoint predictions across 34 ADMET properties. If the worker isn't reachable, the app degrades to Layer 1 alone with a clearly-stated reason — it never fabricates a learned-layer number to fill the gap.

This two-layer design exists specifically so that a heavy, potentially slow model never blocks or risks the app's fast, always-on core feature.

### 2.2 Layer 1 — deterministic drug-likeness (`admet.py::deterministic_profile`)

For each molecule, after standardization (the same `serving.featurize.standardise_smiles` used throughout the app — RDKit cleanup, largest-fragment selection, uncharging):

**Physicochemical descriptors** — molecular weight, LogP (RDKit `Crippen.MolLogP`), TPSA, H-bond donor/acceptor counts, rotatable bonds, aromatic ring count, fraction of sp³ carbons, heavy-atom count, and **QED** (Quantitative Estimate of Drug-likeness — RDKit's `QED.qed`, a single 0–1 score combining several of the above into one common drug-likeness summary metric).

**Four published rule-of-thumb drug-likeness checks**, each a real, named literature rule, computed exactly as defined:
- **Lipinski's Rule of Five**: violations counted from {MW > 500, LogP > 5, HBD > 5, HBA > 10}; `lipinski_pass` = at most 1 violation (the standard "≤1 violation" convention, not "zero violations").
- **Veber**: rotatable bonds ≤ 10 **and** TPSA ≤ 140 — oral-bioavailability-oriented.
- **Egan**: TPSA ≤ 131.6 **and** −1 ≤ LogP ≤ 5.88 — the "egg" absorption model's boundary.
- **Ghose**: MW between 160–480, LogP between −0.4 and 5.6, heavy atoms between 20–70 — a drug-likeness filter derived from known oral drugs.

**Explicit, code-level statement of intent, reproduced because it governs how every downstream feature treats these flags**: *"Informational only; natural products commonly violate these while remaining bioactive. Never used to filter or penalise."* This matters specifically because PhytoScreen's domain is natural-product chemistry, where large, polar, multi-ring-system compounds (glycosides, alkaloids, terpenoids) routinely and legitimately violate these rules while still being real, potent, bioactive molecules — the app is built around not letting a rule-of-thumb silently disqualify real candidates.

**Structural alerts** — three independent, real RDKit `FilterCatalog`s checked against every molecule:
- **PAINS** (Pan-Assay Interference Compounds) — substructures known to cause frequent false positives across unrelated assays.
- **BRENK** — a broader medicinal-chemistry "undesirable substructure" catalog.
- **NIH** — the NIH's own structural-alert catalog.

Every match is real and named (catalog + RDKit's own description string), never invented. Alert count alone does not gate anything; it's informational, same as the rule-of-five checks.

**Cautions** — a short, human-readable list generated from real thresholds: a PAINS match ("may assay-interfere; verify experimentally"), LogP > 6.5 ("solubility/permeability risk"), MW > 800 ("oral absorption less likely — common for glycosides/natural products," an explicit acknowledgment that this specific flag will legitimately fire often in this app's actual compound domain rather than being an unexpected edge case).

### 2.3 Layer 2 — learned endpoints (ADMET-AI)

[ADMET-AI](https://github.com/swansonk14/admet_ai) is a real, published, Chemprop-based multi-task model trained on curated ADMET datasets (the [Therapeutics Data Commons](https://tdcommons.ai/) collection — the endpoint names visible in `admet_endpoints.py`, e.g. `HIA_Hou`, `CYP3A4_Veith`, `hERG`, `AMES`, `LD50_Zhu`, are TDC's own standard dataset identifiers, not this app's invention). It is a genuinely separate, real prediction model, not a relabeling of the deterministic layer.

**34 endpoints mapped to five display groups** (`admet_endpoints.py::ENDPOINTS`), each with:
- **task**: `class` (a classification probability, 0–1) or `reg` (a regression value in a real physical unit).
- **unit**: e.g. `log cm/s`, `%`, `L/kg`, `hr`, `mL/min/kg`, `kcal/mol`, or `prob` for classification tasks.
- **polarity**: `risk` (higher = more concerning — e.g. every CYP inhibition endpoint, hERG, AMES mutagenicity, DILI), `benefit` (higher = more favorable — e.g. human intestinal absorption, oral bioavailability), or `neutral` (a plain PK number with no inherent "good/bad" judgment attached — e.g. half-life, volume of distribution, plasma protein binding).

| Group | Representative endpoints |
|---|---|
| **Absorption** | Human intestinal absorption, oral bioavailability, Caco-2 permeability, P-gp substrate status, PAMPA permeability, aqueous solubility, lipophilicity (LogD), hydration free energy |
| **Distribution** | Blood-brain barrier penetration, plasma protein binding, volume of distribution |
| **Metabolism** | CYP1A2/2C19/2C9/2D6/3A4 inhibition (five separate, real isoform-specific endpoints — a single "metabolism risk" score is never fabricated by averaging them), plus CYP2C9/2D6/3A4 substrate-likelihood |
| **Excretion** | Half-life, hepatocyte clearance, microsomal clearance |
| **Toxicity** | Mutagenicity (AMES), hERG cardiotoxicity, drug-induced liver injury, carcinogenicity, clinical toxicity, skin sensitization, acute toxicity (LD50), **plus the full Tox21 nuclear-receptor and stress-response panel** (12 separate real assay-derived endpoints: androgen/estrogen receptor, aryl hydrocarbon receptor, aromatase, PPAR-gamma, oxidative-stress/genotoxicity/heat-shock/mitochondrial-toxicity/p53-stress-response readouts) |

**Color tone** for a classification endpoint is derived directly from its polarity and two fixed probability thresholds (`RISK_HIGH = 0.70`, `RISK_MED = 0.40`): for a `risk` endpoint, probability ≥ 0.70 is "bad" (red/clay), ≥ 0.40 is "warn" (amber), else "good"; for a `benefit` endpoint the same thresholds apply with the good/bad sense inverted. Regression-task and neutral-polarity endpoints are never color-judged — they display as a plain number in their real unit, since there is no universal "good" volume of distribution or half-life independent of the specific drug-development context.

**DrugBank-approved percentile** — where ADMET-AI provides it, each endpoint also carries a percentile rank against DrugBank-approved drugs (e.g. "this compound's predicted hERG risk is higher than 82% of approved drugs") — real context from the model's own reference distribution, not computed by this app.

**Endpoints not in the curated map** — any column ADMET-AI returns that isn't named in `admet_endpoints.py`'s table is **not currently rendered anywhere** in the grouped UI display (the module docstring's "shown under Other" description is aspirational; the actual `grouped_learned()` code only ever populates the groups explicitly named in `ENDPOINTS`, so an unmapped column is silently absent from the grouped view — though see the "raw" field below, which does carry it).

**Full raw row preserved regardless.** Every column ADMET-AI actually returned for a compound — every physicochemical descriptor, every task prediction, every `*_drugbank_approved_percentile` — is retained verbatim in a `raw` field alongside the curated grouped view, specifically so the **"Download all data (CSV)"** export doesn't need a second server round-trip and never loses information the grouped UI happens not to surface.

### 2.4 Why the learned layer runs in its own process

`admet_service.py` is a **separate FastAPI application**, run as its own process (default port 8100), loading the ADMET-AI model once at startup and keeping it resident in memory. The reasoning, stated directly in the module docstring: *"Runs ADMET-AI in its OWN process so its heavy CPU load is isolated from the main web app (potency/compare tabs stay responsive)."* The main app talks to it purely over HTTP (`ADMET_SERVICE_URL`, default `127.0.0.1:8100`) and never imports ADMET-AI directly — if the worker process isn't running at all, or crashes, or is simply slow, the main app (Predict, Compare, Docking, every other tab) is architecturally incapable of being blocked by it.

**Health check, cached and self-healing** (`admet.py::learned_status`): a `GET /health` probe against the worker, cached for 15 seconds on a positive result but only 3 seconds on a negative one — deliberately asymmetric, so the app re-probes quickly and "self-heals" the moment a worker that was down comes back up, without needing a restart of the main app.

**Sequential job queue inside the worker** (`admet_service.py`): submitted jobs are processed by exactly **one background thread**, one job at a time — *"so multiple big runs queue instead of thrashing all CPU cores at once."* This is a deliberate throughput-vs-fairness tradeoff: a large batch from one user will not be starved or interleaved unpredictably with another's, at the cost of queuing behind whatever's already running.

### 2.5 Synchronous vs. asynchronous dispatch

The main app (`app.py::admet`) picks the request mode based on batch size, transparently to the caller:

1. **Worker unavailable, or nothing in the batch parsed successfully** → deterministic-only result, returned immediately (`mode: "result"`).
2. **≤ 50 valid molecules** (`ADMET_SYNC_MAX`) → a single synchronous call to the worker's `/profile` endpoint, returned immediately once it responds.
3. **> 50 valid molecules** → an asynchronous job (`worker_submit` → `/jobs`), returning a `job_id` the frontend polls every 1.5 seconds until `done`, showing live progress (compounds processed / total) rather than one static wait.

The worker itself batches internally at 512 molecules per model call regardless of how the request arrived (`BATCH` in `admet_service.py`), and rejects a synchronous `/profile` call outright above 1,000 molecules (directing the caller to the async `/jobs` path instead) — a batch that large is assumed to be a job, not an inline request.

---

## 3. Architecture — module map

| Module | Responsibility |
|---|---|
| `admet.py` | Runs inside the main app. Owns the deterministic layer entirely (§2.2); talks to the worker over HTTP for the learned layer; merges the two into one profile per compound; health-check caching; sync/async client logic. |
| `admet_endpoints.py` | The declarative endpoint map (§2.3) — group, human label, task type, unit, and risk/benefit/neutral polarity for every named ADMET-AI column. Editing this file is the entire mechanism for adding, renaming, or regrouping an endpoint; no other code needs to change. |
| `admet_service.py` | The isolated worker process — loads ADMET-AI once, exposes `/health`, `/profile` (sync), `/jobs`+`/jobs/{id}` (async), runs jobs through one sequential background-thread queue. |
| `app.py` (`/api/admet`, `/api/admet/job/{jid}`) | The public API surface: decides deterministic-only vs. sync vs. async based on batch size and worker availability (§2.5); tracks in-flight async jobs. |
| `analysis.py` | Reuses `admet.admet_profile()` per-compound to populate the Compare tab's ADMET view (documented in `QSAR_BIOACTIVITY_PREDICTION.md` — the same deterministic+learned profile, just consumed by a different tab). |
| `research_report.py` | Also reuses `admet.admet_profile()` for the ADMET section of the "Research Report" evidence-chain document (documented in `DOCKING.md` §7). |

### Frontend

| File | Responsibility |
|---|---|
| `tabs/AdmetTab.tsx` | The tab itself: molecule input, submission (sync or job-polling), the results table, per-compound expandable endpoint detail, full-data CSV export. |
| `components/MoleculeInputPanel.tsx` + `lib/useMoleculeInput.ts` | The same shared three-mode input (paste/CSV/SDF) used by Predict and Compare. |

---

## 4. The user journey

### 4.1 Input

Identical to Predict/Compare: **Paste** (one SMILES per line), **CSV file** (client-side parsed, `smiles`/`SMILES` header preferred, first column as fallback), or **SDF file** (server-side RDKit parse via `/api/parse_sdf`). There is no target selection on this tab at all — ADMET profiling is target-independent by design, and the sidebar states this directly ("Target-independent").

### 4.2 "Profile compounds"

Submits to `/api/admet`. Depending on what comes back (§2.5):
- **Immediate result** — deterministic-only (worker unavailable) or a small synchronous batch.
- **A progress bar** — "Profiling N compounds in the ADMET-AI worker…" with a live percentage, for a batch large enough to be dispatched as an async job.

### 4.3 Results

- If the learned layer isn't available for this run, a notice states exactly why (the same real reason string from `learned_status()` — e.g. instructions to start the worker — never a generic "unavailable").
- A disclaimer, always shown: drug-likeness flags are informational, natural products often violate them while remaining bioactive, they are never used to filter compounds — plus a note on whether learned ADMET-AI endpoints are included in this particular run.
- **Download all data (CSV)** — exports every raw column ADMET-AI returned (not just the curated grouped subset), one row per submitted compound, including unparsed or worker-missing compounds (smiles only, columns blank) rather than silently dropping them from the export.

**Results table**, one row per compound: molecular weight, LogP, QED, a Lipinski-violations badge (color-coded pass/fail), a structural-alert count badge, and — only if the learned layer ran — a "Tox flags" badge (count of learned endpoints this specific compound tripped into the "bad" tone).

**Expanding a row** (available whenever the learned layer produced grouped data for that compound) reveals every endpoint, grouped under Absorption/Distribution/Metabolism/Excretion/Toxicity headers, each shown as a small colored chip: label, the real value in its real unit (or a rounded percentage for a classification probability), and — where available — its DrugBank-approved percentile. Hovering any chip explains exactly what it is (a learned Chemprop prediction, its unit, and the percentile's meaning) rather than leaving the number to be guessed at.

---

## 5. Background automatic work and decisions — consolidated

### 5.1 Always automatic

- SMILES standardization before any descriptor is computed (same standardizer as the rest of the app).
- All four rule-of-thumb drug-likeness checks and all three structural-alert catalogs, for every parseable molecule.
- The worker health probe before deciding sync/async/deterministic-only routing — the user never manually declares "the worker is up."
- Batch-size-based routing (≤50 sync / >50 async / worker down → deterministic-only) — invisible beyond the resulting wait experience.
- Internal 512-molecule batching inside the worker, regardless of how the request arrived.
- Color-tone assignment for every classification endpoint, from its declared polarity and the fixed 0.70/0.40 probability thresholds.
- Health-check caching with the asymmetric 15s-positive/3s-negative TTL, so a worker coming back online is detected quickly without any restart.

### 5.2 Automatic by default, not currently user-adjustable

- The risk color thresholds (`RISK_HIGH`/`RISK_MED`) — fixed constants, not exposed as settings.
- The sync/async cutover point (50 compounds) and the worker's own hard `/profile` limit (1,000).
- Which endpoints appear in which group, and their labels/units/polarity — entirely driven by `admet_endpoints.py`, a code file, not a runtime setting.

### 5.3 Never automatic

- Which molecules are submitted, and by which input method.
- Starting the ADMET-AI worker process at all — this is an operator/deployment decision (`uvicorn admet_service:app --port 8100`, requiring `pip install admet-ai`), not something the main app can do for itself. Its absence degrades the feature gracefully rather than failing.

---

## 6. Limitations — stated directly

- **Drug-likeness rules and structural alerts are heuristics from historical drug/screening-hit populations, not universal correctness criteria** — the app's own code explicitly rejects using them as filters for exactly this reason, especially given natural products routinely and legitimately violate several of them.
- **The learned layer's predictions are model outputs, not measurements** — each of the 34 endpoints is a real, separately-trained prediction with its own accuracy characteristics (not documented per-endpoint within this app; consult ADMET-AI's own published validation for per-task performance). This app does not currently surface a per-endpoint confidence or error metric the way the QSAR potency models do (Test R²/RMSE) — an ADMET-AI prediction is shown with its value and DrugBank-approved percentile context, but not a stated accuracy figure.
- **An endpoint's DrugBank-approved percentile describes where a compound falls relative to already-approved drugs on that one endpoint — it says nothing about the other 33 endpoints, and does not imply overall drug-likeness by itself.**
- **CYP-isoform predictions are five genuinely separate endpoints** (1A2/2C19/2C9/2D6/3A4) — there is no single combined "metabolism risk" score; a compound flagged on one isoform and clean on the other four is reported exactly that way, never collapsed into one number.
- **Structural-alert catalogs (PAINS/BRENK/NIH) catch known problematic substructure *patterns*, not all possible liabilities** — a clean alert count is not proof of safety, only the absence of a match against these three specific, real but non-exhaustive catalogs.
- **The worker's own model version/training-data snapshot is whatever `pip install admet-ai` resolves to at deployment time** — this app does not pin or independently verify a specific ADMET-AI model version; users concerned with exact reproducibility of learned-layer numbers should check the ADMET-AI package version in their own deployment.
- **Any ADMET-AI output column not named in `admet_endpoints.py` is invisible in the grouped UI view** — present only in the raw CSV export, not in the interactive table, so a newly-added ADMET-AI endpoint (from a package upgrade) needs a corresponding entry added to `admet_endpoints.py` before it will appear anywhere but the raw download.

---

## 7. API reference (ADMET-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `POST /api/admet` | Submit a SMILES list. Returns an immediate result (`mode: "result"`) or a job id to poll (`mode: "job"`), depending on batch size and worker availability. |
| `GET /api/admet/job/{jid}` | Poll an async ADMET job's status/progress/result. |
| `GET /health` *(on the worker, port 8100 by default — not part of the main app's API surface)* | Whether the ADMET-AI model loaded successfully in the worker process. |
| `POST /profile` *(worker)* | Synchronous prediction for ≤1000 molecules. |
| `POST /jobs` + `GET /jobs/{id}` *(worker)* | Asynchronous prediction job submission/polling, used internally by the main app for large batches. |

---

## 8. Document provenance

Written by reading, in full: `backend/admet.py`, `backend/admet_endpoints.py`, `backend/admet_service.py`, the `/api/admet` and `/api/admet/job/{jid}` routes in `backend/app.py`, and `frontend/src/tabs/AdmetTab.tsx`. No content here was reconstructed from memory of past conversation — every claim traces to a specific line of code read during this pass.
