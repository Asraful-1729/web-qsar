# Screen — Complete Technical & Scientific Documentation

**Scope:** the Screen tab — the combined pipeline that runs **QSAR potency prediction + ADMET profiling + Docking** against one target in one submission and fuses them into a single ranked shortlist. This document assumes familiarity with `DOCKING.md` and `QSAR_BIOACTIVITY_PREDICTION.md` and does **not** re-explain how Vina docking, PoseBusters, GNINA, redocking validation, or the Chemprop+AutoGluon QSAR pipeline work internally — only Screen's own orchestration, fusion, and integration logic, which is genuinely new and not documented elsewhere. Written directly from `backend/serving/screen.py`, the `/api/screen/*` routes in `backend/app.py`, and `frontend/src/tabs/ScreenTab.tsx`.

---

## 1. Objective

Run **one submission** — a target and a list of compounds — through every independent line of evidence the app can produce for that target (potency, drug-likeness/toxicity, physical docking pose and score), and return **one ranked shortlist** with every underlying signal still visible per row, rather than requiring a researcher to run three separate tabs and manually cross-reference the results. This is the app's headline feature — the module's own docstring calls it exactly that: *"the headline feature: an explicit, ordered 8-step pipeline the researcher can watch run, ending in a ranked, honestly-caveated shortlist with every signal attached."*

---

## 2. The 8-step pipeline

`serving/screen.py::run()` executes these in strict order, reporting each transition to the UI live:

| Step | What happens |
|---|---|
| 1 | Parse & standardize every submitted SMILES (RDKit cleanup, largest-fragment, uncharge — the same standardizer used everywhere in the app). Unparsable inputs are recorded as `skipped`, never silently dropped without a trace. |
| 2 | Featurize (RDKit 2D descriptors + MACCS + Morgan/ECFP4) — reported as its own visible step for transparency, but **actually executed inside step 4's `predict_smiles` call**, not as a separate pass; the code comment is explicit about this: *"steps 2-4 share one featurisation pass for correctness + speed."* |
| 3 | Applicability-domain check (mean |z| ≤ 3.0 against the target's training chemistry) — see `QSAR_BIOACTIVITY_PREDICTION.md` §2.4 for the method; not re-derived here. |
| 4 | QSAR potency prediction — the full two-stage Chemprop→AutoGluon pipeline, one call for the whole batch. |
| 5 | ADMET profiling — deterministic layer always; the learned (ADMET-AI) layer if its worker is reachable, exactly as in `ADMET.md`. |
| 6 | Docking — **only if a usable receptor profile exists for this target** (see §3). This is the slow step; progress is reported per compound, not just per pipeline-step. |
| 7 | Rank & fuse the evidence into the final shortlist (§4). |
| 8 | Finalize — attach a real, data-derived methods/caveats note (§6). |

The frontend mirrors this exactly (`SCREEN_STEPS` in `ScreenTab.tsx`), showing a live numbered checklist with the current step highlighted and completed steps checked off — not a generic spinner.

---

## 3. Whether docking runs at all — resolved once, per submission

Screen does not assume docking is available. Before step 6 can run, the pipeline resolves a docking profile using the **exact same logic** `/api/docking/submit` uses (Advanced Settings' `custom_profile` if the user picked a manual structure, otherwise the registry's automatic default for the target; Blind-mode box computed on demand if selected; an explicit box override applied if given). Docking is skipped, cleanly and with a stated reason (`docking_note`), in every one of these cases — never a crash, never a silently-empty docking column:

- The docking package itself isn't importable on this machine.
- The docking engine (Vina) isn't available/ready.
- No profile exists for this target at all (a `custom_profile` whose `target_id` doesn't match, or no registry entry and no manual pick).
- The target has no site-specific binding box on record and Blind mode wasn't selected.
- Blind mode was selected but no receptor file exists on disk to compute a whole-protein box from.

When docking *does* run, it uses whatever Advanced Settings values are in effect (exhaustiveness, pose count, GNINA on/off) — identical settings to what a plain Docking-tab submission with the same Advanced Settings would use, so a Screen result and a separate Docking-tab result for the same target/settings are directly comparable.

**Redocking validation runs automatically here too, exactly as in the plain Docking tab** (`DOCKING.md` §2.7) — once, after every compound has been docked, using the same structure/exhaustiveness/pose-count just used for the batch. It never blocks or fails the shortlist if it errors; the result (or the error) is carried in the response as `redocking_validation` and rendered by the same `RedockingValidationNote` component the Docking tab uses.

---

## 4. The fusion formula — how QSAR and docking scores become one ranked list

This is the part of Screen that doesn't exist anywhere else in the app, so it's documented in full.

### 4.1 Why raw scores are never averaged directly

A predicted pIC50 (typically 4–10) and a Vina score (typically −4 to −12 kcal/mol, more negative is better) are **on incompatible scales with opposite directionality**. Averaging them directly would be meaningless — and would let whichever value happens to have the larger numeric range silently dominate the ranking. Screen never does this.

### 4.2 Rank-based normalization (`_rank_score`)

Both QSAR and docking values are independently converted to a **0–1 rank-based score within the current batch only** — never a global or cross-run scale:

- Every compound with a real (non-`None`) value is ranked against the others in *this* submission.
- The best-ranked compound gets score **1.0**, the worst gets **0.0**, with every other compound's score its linear position between them (`1.0 − rank/(n−1)`).
- QSAR uses `higher_is_better=True` (higher pIC50 = more potent); docking uses `higher_is_better=False` (more negative Vina score = better binding).
- A compound with no value at all for a given axis (out-of-domain QSAR, or docking that didn't produce a pose) gets `score = None` on that axis, not a synthetic worst-case value — its absence is tracked explicitly, not encoded as a fabricated 0.0.

### 4.3 The fused score

```
fused = (1 − DOCK_WEIGHT) × qsar_rank_score + DOCK_WEIGHT × docking_rank_score      (both present)
fused = qsar_rank_score                                                              (docking didn't run, or no pose for this compound)
```

`DOCK_WEIGHT = 0.35` — a fixed, disclosed constant (35% docking / 65% QSAR when both are available). It is not learned, not per-target, and not currently exposed as a user setting. **When docking did not run for a compound at all** (the whole batch had no docking profile, or this specific compound simply produced no valid pose), the fused score is the **QSAR rank score alone** — docking's 35% weight is not redistributed, reduced, or otherwise renormalized; it simply isn't part of the calculation for that row, and the methods note (§6) states plainly whether docking contributed to ranking at all for the run.

Compounds with **no fused score at all** (QSAR itself was out-of-domain or unparsable, with no docking to fall back on) sort to the bottom, after every compound that has a real score — never interleaved as if a missing score were a low-but-real one.

### 4.4 Per-row caveats — generated, not templated boilerplate

Each shortlist row carries a real, specific caveat list, built from that row's own actual values, never a generic disclaimer repeated identically for everyone:

- `"Outside the QSAR model's training chemistry — potency not trusted."` — only when this compound's own AD check failed.
- `"Docking: {the actual reason}."` — only when this compound's own docking attempt didn't produce a usable result, quoting its real failure reason.
- `"{n} structural alert(s) (PAINS/Brenk/NIH) — informational, not a filter."` — only when this compound actually matched one, with its real count.
- The Blind-docking caveat — only appended to rows that actually have a docking result, when the run used Blind mode.

---

## 5. Gene-only targets — Screen degrades to plain docking, deliberately

A `GENE_<symbol>` target (a protein with disease-association evidence but no trained QSAR model — see `QSAR_BIOACTIVITY_PREDICTION.md`/`DOCKING.md` for how these arise) has nothing for steps 3/4 (AD/QSAR) to run against. Rather than error, the frontend detects this (`isGeneOnly(targetId)`) **before** submission and routes the entire request through the plain `/api/docking/submit` pathway instead of `/api/screen/submit` — the button itself relabels to **"Dock this target (no QSAR model)"** so this isn't a silent fallback the user has to infer. The result renders through `GeneOnlyDockResults`, which is **literally the same results-table shape the plain Docking tab uses** (reusing `DockDetailPanel`, `FreshDecoyButton`, `InteractionLogTable`/`InteractionTableToggle`, and the Docking-tab export endpoints), not a separate implementation — the code comment states this directly: *"reuses the Docking tab's own result shape/rendering rather than a hardcoded id, matching the original app's routing."*

---

## 6. The methods note — real numbers, not a fixed disclaimer

Every completed Screen run attaches a `methods_note` assembled from that run's own actual facts, not a static string:

- How many compounds were submitted and how many were dropped as unparsable.
- The applicability-domain threshold actually used (mean |z| ≤ 3.0) and that only in-domain compounds carry a trusted potency value.
- **Which docking scenario applied** — one of three real, mutually-exclusive sentences, chosen from what actually happened this run: docking ran in Blind mode (with the Blind caveat inlined) and contributed 35% of the fused score; docking ran site-specifically and contributed 35%; or docking did not run at all, with the real reason quoted from `docking_note`.
- That ADMET flags are informational, never used to exclude a compound.
- The standing disclaimer: *"This is a prioritisation aid, not a substitute for assays."*

---

## 7. The user journey

### 7.1 Sidebar

Identical setup to the Docking tab (`TargetBrowser` with `need: ["model", "docking"]` — a target must be reachable for *both* kinds of data for the download gate to consider it fully ready; `DockingModeSection`, `AdvancedSettingsPanel`), plus:

- **"Why this?"** (`WhyThisButton`, from `RecommendationPanel.tsx`) — shown only for real (non-gene-only) targets. Fetches and displays the same structural-evidence bundle documented in `DOCKING.md` §3 (`docking/recommend.py`'s crystallographic-quality context: resolution, RSCC/RSR, rank among qualifying structures) — explaining *why* this target's automatic default structure was chosen, one click away, without leaving the sidebar.
- Molecule input (paste/CSV/SDF, the same shared component as every other prediction tab) and the same optional Plant source field, threaded into export metadata identically to the Docking tab.

**"Run screen"** submits the full pipeline; for a gene-only target it reads **"Dock this target (no QSAR model)"** instead (§5).

### 7.2 Live progress

A vertical 8-step checklist (§2's table), the current step highlighted, completed steps checked. During step 6 specifically, the note switches to live per-compound docking progress (`"Docking 3/12…"`) rather than staying on the generic step label, since that's the step slow enough to need it. A **Stop** button is present throughout (cooperative cancellation, same pattern as every other long-running job in this app — checked once per compound during docking, the only per-item checkpoint the pipeline has).

### 7.3 Results

A header strip (submitted/parsed/skipped counts, whether docking was used at all), the real methods note (§6), and — if docking was skipped — the actual reason, shown as a notice rather than left to be inferred from an empty column. The redocking-validation note (identical component to the Docking tab) appears next, when applicable.

**The shortlist table**, one row per compound, ranked by fused score: rank, compound, QSAR pIC50 (or "out-of-domain," never a suppressed-but-implied number — consistent with the hard AD gate documented in `QSAR_BIOACTIVITY_PREDICTION.md` §2.4), QSAR confidence chip, **Vina score + docking confidence columns (only shown at all if docking was used for this run)**, **GNINA columns (only shown if any row actually has GNINA data)**, the fused score, a per-row caveats list (§4.4), and a Fresh Decoy Validation button (identical feature and methodology to `DOCKING.md` §2.10, linked back to this Screen job via `parentKind: "screen"` so a validated compound's result is included in the export package).

Expanding a row (available whenever it has a docking result worth showing) reveals the **exact same per-compound detail panel** the Docking tab uses (`DockDetailPanel` — interaction diagram with lightbox, 3D pose viewer, summary stat tiles, interaction table) — Screen does not maintain a separate, parallel results-detail implementation.

**Export row**: interaction-table toggle (cross-compound log, identical component to Docking); **CSV** (rank, both input and standardized SMILES, predicted pIC50 + in-domain flag, QSAR confidence, Vina score + docking confidence, fused score, and the semicolon-joined caveats list — Screen-specific, richer than a plain docking export since it carries the fusion columns); and, only when docking actually ran, **Failure log (.csv)** and **Full experiment package (.zip)** — the same export machinery documented in `DOCKING.md` §5.6, adapted here to read each compound's docking result from `shortlist[i].docking` rather than being the row itself.

---

## 8. Background automatic work and decisions — Screen-specific only

(For everything shared with plain Docking or plain QSAR prediction — ligand prep, PoseBusters, confidence tiering, AD gating, redocking validation mechanics — see the other two documents; not repeated here.)

### 8.1 Always automatic

- The docking-availability resolution in §3 — the user never explicitly declares "use docking" or "skip docking" for a Screen run; it's determined from what's actually available for the picked target and settings.
- Sharing one featurization pass across steps 2–4.
- Rank-based (never raw-unit) normalization before fusion.
- Per-row caveat generation from each row's own real values.
- The methods note's docking-scenario sentence selection (one of three, chosen from what actually happened).
- Redocking validation, exactly as in plain Docking.
- Gene-only detection and routing to the plain-docking code path, before submission.

### 8.2 Automatic by default, not currently user-adjustable

- `DOCK_WEIGHT = 0.35` — a fixed constant, not a per-run or per-target setting.

### 8.3 Never automatic

- Which target and compounds to submit; every Advanced Settings override (exhaustiveness, poses, GNINA, box, manual structure) — identical to plain Docking, since Screen reuses the exact same `useAdvancedDocking` state.
- Running Fresh Decoy Validation for any specific compound.
- Downloading any export artifact.

---

## 9. A6 reproducibility — Reproduce / alternate-ligand / research report, Screen's versions

Screen carries the **same three "present in the API, not currently wired into this tab's UI" features** documented in `DOCKING.md` §7, in Screen-shaped form:

- **`POST /api/screen/job/{jid}/reproduce`** — resubmits a finished Screen job with its exact saved parameters. If the original run used docking, the exact receptor file is pinned via `custom_profile` (identical contract to `docking_reproduce`); if the original run was QSAR-only (docking didn't run), reproduction stays QSAR-only too, faithfully.
- **`POST /api/screen/job/{jid}/alternate_ligand/build`** + **`/submit`** — re-center a finished Screen job's structure on a different real co-crystallized ligand from the same raw PDB, then resubmit the full pipeline with every other setting pinned. Reuses the same `_CUSTOM_RECEPTOR_JOBS` machinery and polling endpoint as the Docking tab's version — the receptor-build step doesn't care which tab originally created the job.
- **`POST /api/screen/job/{jid}/research_report`** — the same evidence-chain report as `DOCKING.md` §7, adapted for Screen's shortlist shape (`row["docking"]` nested rather than being the row itself).

All three are fully functional server-side; none currently has a rendered button in `ScreenTab.tsx`, for the same product-design-history reason documented in `DOCKING.md` §7 (not a sign of being broken or abandoned).

---

## 10. API reference (Screen-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `POST /api/screen/submit` | Submit a Screen run (target + compounds + optional Advanced Settings/plant source). |
| `GET /api/screen/job/{jid}` | Poll status — includes the current step (1–8), step label, and (during docking) per-compound progress. |
| `POST /api/screen/cancel/{jid}` | Request cancellation (checked once per compound during the docking step). |
| `GET /api/screen/job/{jid}/export.csv` | The fused shortlist as CSV (§7.3). |
| `GET /api/screen/job/{jid}/export_package` | Full experiment ZIP, only meaningful once docking ran. |
| `GET /api/screen/job/{jid}/interaction_diagram` | On-demand SVG/TIFF/PNG/PDF interaction-diagram export for one shortlist compound. |
| `GET /api/screen/job/{jid}/failure_log` | Docking-failure-only CSV for the run. |
| `POST /api/screen/job/{jid}/reproduce` | §9 — API-only. |
| `GET /api/screen/job/{jid}/alternate_ligands`, `POST .../alternate_ligand/build`, `POST .../alternate_ligand/submit` | §9 — API-only. |
| `POST /api/screen/job/{jid}/research_report` | §9 — API-only. |

---

## 11. Document provenance

Written by reading, in full: `backend/serving/screen.py`, every `/api/screen/*` route in `backend/app.py`, `frontend/src/tabs/ScreenTab.tsx`, and `frontend/src/components/RecommendationPanel.tsx`. Deliberately does not restate content already covered in `DOCKING.md` or `QSAR_BIOACTIVITY_PREDICTION.md` — see those documents for the mechanics of docking and QSAR prediction themselves. No content here was reconstructed from memory of past conversation; every claim, especially the fusion formula in §4 and the methods-note logic in §6, traces to a specific line of code read during this pass.
