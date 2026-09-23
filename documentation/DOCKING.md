# Docking — Complete Technical & Scientific Documentation

**Scope:** every part of PhytoScreen's structure-based docking feature — the scientific method behind it, every button and field in the UI, every automatic background decision, every API endpoint, and its known limitations. Written directly from the current source code (`backend/docking/`, `backend/app.py`'s docking routes, `frontend/src/tabs/DockingTab.tsx` and its supporting components) — not from memory or design intent, so it reflects what the software actually does today, including the parts that are dead code, unvalidated, or deliberately simplified.

**Audience:** anyone who needs to understand *why* a docking result looks the way it does — a researcher interpreting output, a developer maintaining the code, or a reviewer auditing the method.

---

## 1. Objective

Given a target protein (a binding pocket, defined by a 3D structure) and a set of candidate small molecules (SMILES strings), predict which molecules are likely to bind, in what pose, and with what relative strength — using **physics-based molecular docking**, not a machine-learned predictor. This complements the QSAR (ligand-based, trained-on-activity-data) prediction elsewhere in the app: docking reasons from 3D structure and physical interaction geometry, QSAR reasons from historical structure-activity patterns. Neither is a substitute for the other; the app surfaces both as independent lines of evidence.

**What this feature is not:** a replacement for experimental validation. Every score, pose, and confidence label produced here is a computational hypothesis. The interface and this document say so explicitly and repeatedly, rather than letting a number imply more certainty than it has.

---

## 2. Scientific and methodological basis

### 2.1 The core pipeline

```
Ligand (SMILES)
   → 3D conformer generation (RDKit ETKDG + MMFF)
   → PDBQT conversion (Meeko)
   → AutoDock Vina search against a prepared receptor + binding box
   → PoseBusters physical-validity gate (filters poses, does not rank them)
   → best valid pose selected; self-consistency measured against the other valid poses
   → GNINA CNN rescoring (optional second opinion — pose-quality confidence + predicted affinity)
   → confidence tier assigned (validity + self-consistency + CNN agreement)
   → protein–ligand interaction detection (PLIP if installed, else a built-in distance-based fallback)
   → 2D interaction diagram rendered (LigPlot+-style)
```

Every step above runs once per compound, automatically, with no separate user action required beyond submitting the batch.

### 2.2 Why AutoDock4 is gone

An earlier version of this module ran both AutoDock Vina and AutoDock4 and treated agreement between the two as a confidence signal ("cross-engine consensus"). AutoDock4 has been removed entirely: it is GPL-licensed (a distribution concern for this project), and a second physics-based scoring function adds little independent signal over Vina alone for the purpose of ranking. The single-engine pipeline's confidence signal is now built from three genuinely independent sources instead (§2.4).

### 2.3 PoseBusters — a gate, not a scorer

[PoseBusters](https://github.com/maabuu/posebusters) checks each candidate pose for physical plausibility: correct stereochemistry, sane bond lengths/angles, non-flat rings, no internal atomic clashes, and (when a receptor is supplied) no protein–ligand clashes. **This is a filter applied before ranking, never a score used to rank poses against each other.** A pose that fails any check is dropped; the best-scoring pose is then chosen only from among the poses that passed. If PoseBusters is not installed, the gate is bypassed — `check()` returns `(True, {"posebusters": "not installed — validity not checked"})` — this is a fail-*open* behavior on missing-tool, explicitly distinct from a fail-*closed* behavior on an actual check error (an exception during a real check returns `(False, {"error": ...})`, i.e. the pose is rejected). The distinction matters: "not checked" and "checked and failed" are not the same state, and the code keeps them separate rather than collapsing both to "invalid."

### 2.4 Confidence — three independent signals, no cross-engine agreement

With AutoDock4 gone, confidence for a single-engine (Vina) result is built from:

1. **Validity** — did any pose pass PoseBusters at all?
2. **Self-consistency** — do Vina's own top poses (from `--num_modes`) cluster together? Specifically: how many of the *other* physically-valid poses land within `rmsd_threshold` (default 2.0 Å) of the best one. This is Vina agreeing with itself across its own independent samples, not proof of correctness, but a real signal that the search converged rather than scattering.
3. **GNINA CNNscore** — a convolutional neural network's independent opinion on whether the pose looks like a real binding pose (0–1, ≥0.5 counted as "good").

Combination rule (`docking/consensus.py::assign_confidence`):

| Self-consistent? | Good CNN score? | GNINA available? | Confidence |
|---|---|---|---|
| — | — | no valid pose at all | **none** |
| yes | yes | yes | **high** |
| yes or | yes | yes | **medium** |
| no | no | yes | **low** |
| yes | — | **no (GNINA not installed/disabled)** | **medium** (capped) |
| no | — | no | **low** |

**GNINA capping is deliberate, not a missing feature**: without a second opinion, the pipeline refuses to claim "high" confidence no matter how self-consistent Vina's own search looks, because self-consistency alone cannot distinguish "converged on the right answer" from "converged on the same wrong answer." This is stated directly in the module's own design notes.

### 2.5 GNINA — CNN rescoring, an optional second opinion

[GNINA](https://github.com/gnina/gnina) runs a trained convolutional neural network directly on the receptor–ligand complex geometry, independent of Vina's own scoring function. It reports:
- **CNNscore** (0–1): the network's confidence that this looks like a genuine binding pose.
- **CNNaffinity**: a predicted, pKd-like binding-affinity value.
- **Affinity**: a Vina-style kcal/mol value computed by GNINA's own (also CNN-influenced) scoring.

GNINA is optional (Advanced Settings → "GNINA CNN rescoring" checkbox, on by default when installed). If unavailable or disabled, the pipeline proceeds without it and confidence is capped at `medium` as above — this is never silently substituted with a fabricated number.

### 2.6 Interaction detection — PLIP (gold standard) or a built-in fallback

Two independent detectors exist, selected automatically:

- **PLIP** (Protein–Ligand Interaction Profiler), if installed: real atom-level interaction typing — hydrogen bonds, π–π stacking (parallel or T-shaped), hydrophobic (alkyl / π-alkyl, distinguished by whether the contacting protein atom is actually part of an aromatic ring), π-cation, salt bridges, halogen bonds — every entry traces to a real detected atom pair, nothing invented.
- **Built-in distance-based fallback**, used only when PLIP is not installed: a much simpler rule (H-bond if both atoms are polar (N/O) within 3.6 Å; hydrophobic if both are carbon within 4.2 Å). This fallback only distinguishes two interaction types and carries no atom-level detail — its output is explicitly labeled `"(distance-based)"` in every interaction name so the results table never implies a precision the method doesn't have.

Every result records which detector produced it (`interaction_source`: `"PLIP"` or `"distance-based"`), and this is shown in the UI next to the interaction diagram.

### 2.7 Redocking validation — the built-in sanity check

Every docking or screening submission that has a real co-crystallized reference ligand automatically re-docks that reference ligand against the same receptor, structure, and search parameters just used for the user's own compounds, and computes the RMSD (root-mean-square deviation, atom-map-safe — see §2.8) between the redocked pose and the ligand's real experimental (crystal) pose. **< ~2 Å is the conventional "validated" threshold** used throughout: it means Vina's search, on this specific structure and box, can reproduce a pose already known to be physically correct. This is a check on the *docking setup* (receptor prep, box placement, search parameters), not a score for any compound the user submitted — the UI states this explicitly ("Sanity check on this structure's own docking setup … not a score for any compound you submitted").

This used to be a separate, opt-in "validate this target" action. It no longer is: it now runs automatically as part of every submission, using whatever structure and parameters the user already picked, because the validation is only meaningful for the *exact* setup being used for real compounds — a stale validation from different settings would be misleading. Redocking never blocks or delays the display of compound results; if it errors, the error is recorded and shown, but the compound results are unaffected.

Determinism: ordinary compound docking leaves Vina's search seed unset (a fresh random search each run — appropriate for exploring pose space). Redocking validation instead fixes both the Vina seed (`0xf00d`) **and** forces `--cpu 1`. The seed alone was tested and found insufficient — Vina's multi-threaded search still drifted by roughly 0.1 Å between identical runs even with a fixed seed; only single-threading removed that residual non-determinism. This matters because a "validated" badge should mean the same thing on every run of the same structure, not fluctuate with parallelism.

### 2.8 RMSD — atom-map-safe by construction

A naive RMSD calculation compares atoms by their index order in each file. This is unsafe across independently-generated structures: AutoDock Vina's own PDBQT output and a molecule reconstructed from crystal-structure PDB coordinates are not guaranteed to list the same atoms in the same order, even though they represent the same molecule. An index-wise RMSD on mismatched atoms produces a number that looks valid but is meaningless — a documented risk from the original docking-pipeline design review.

The fix (`docking/rmsd.py`): before computing any RMSD, the two poses' heavy-atom molecular graphs are compared (via canonical SMILES, stereochemistry stripped since a docked pose may not carry stereo information) to confirm they really are the same molecule. If they are not, `safe_rmsd()` **raises an error rather than returning a number** — it never silently produces a misleading value. When they do match, RDKit's `GetBestRMS` is used, which searches over molecular symmetry to find the true best atom correspondence, rather than trusting file order.

### 2.9 Receptor preparation — the full pipeline, once per structure

Before any docking can happen against a given PDB structure, it must go through preparation (`docking/receptor_prep.py`):

1. **Extract the reference ligand.** Every HETATM group in the structure is a candidate; a hand-curated blacklist (`scripts/select_receptor.py::ADDITIVE_BLACKLIST`, ~90+ entries) excludes crystallization additives, buffers, ions, cryoprotectants, glycosylation sugars, heavy-atom phasing derivatives, and non-hydrolyzable nucleotide analogs used as cofactor mimics. Among what remains, the largest (by heavy-atom count) group is picked automatically — unless the user (or an internal curated pick) specifies exactly which one.
2. **Strip to protein only.** Waters and all heteroatoms removed; if the structure carries multiple copies of the protein (a crystallographic dimer in the asymmetric unit — common), only one chain is kept, since docking against two overlapping copies is both physically wrong and can confuse bond perception at the chain–chain interface.
3. **Repair (PDBFixer).** Adds missing atoms, missing residues, and hydrogens at pH 7.0. A short vacuum energy minimization is then applied *only* to relieve local clashes PDBFixer's atom placement can introduce (confirmed case: two residues' rebuilt sidechains landing ~1.56 Å apart when they should be ~2.3–2.4 Å) — this is not a conformational refinement, just enough to stop a clash from being misread as a bonded atom downstream.
4. **Convert to PDBQT (Meeko).** Assigns AutoDock atom types and partial charges; the receptor is kept rigid (no torsions) for docking.
5. **Compute the binding site.** Pocket residues within 5 Å of the reference ligand are listed (for display and for the residue-based box definition), and a grid box is computed centered on the ligand (padding 8 Å by default, minimum size 20 Å in each dimension).

Every step's real, measured before/after facts (atom counts, residues removed, box dimensions) are recorded in a `prep_report` — never a generic label restating what the step does, but the actual numbers that step produced. This report is downloadable from the Manual Structure picker.

**Sanity checks, not silent rejections:** if the automatically-picked ligand is unusually small (<5 heavy atoms), unusually light (<180 Da) or heavy (>2000 Da), or its "pocket" has fewer than 3 nearby receptor residues, these are flagged as `ligand_sanity_warnings` and the whole profile is marked `build_status: "review_required"` — the receptor still builds (an automated pipeline shouldn't hard-fail on an edge case), but the flag is real and surfaced in the UI, never silently swallowed.

**No-ligand fallback:** if a structure genuinely has no qualifying small-molecule ligand (only buffers/ions/cryoprotectants), the receptor still builds, just as **Blind-only** (`site_source: "none_validated"`) — there is no pocket to center a site-specific box on, but whole-protein docking is still offered.

### 2.10 Fresh Decoy Validation — DUD-E-style enrichment, on demand

For any single compound worth deeper scrutiny, "Run Fresh Decoy Validation" generates a fresh set of decoys tailored to *that specific compound* and asks: does this compound's docking score beat those decoys?

Methodology (`scripts/generate_decoys.py`, adapted from Mysinger et al. 2012's DUD-E recipe):
1. **Property-matched**: a decoy must resemble the real compound's molecular weight (±30), LogP (±1.0), H-bond donors (±1), H-bond acceptors (±2), rotatable bonds (±2), and formal charge (exact) — so it's a physically plausible, similarly drug-like molecule, not a random reject.
2. **Topologically dissimilar**: a decoy must have Morgan/ECFP4 Tanimoto similarity < 0.25 to the query compound — otherwise it would trivially dock well by being nearly the same molecule, defeating the point.

Candidate pool: real ChEMBL compounds curated for the app's *other* ~64 QSAR targets (never the target under test), excluding anything that also appears in the compound's own training data. **Disclosed limitation, stated directly in the code and worth restating here:** this project has no external decoy database (e.g. ZINC) wired in. "Decoy" here means "a real compound curated for a different target's assay," not an independently-confirmed non-binder. This is a real, useful comparison (does this compound's score beat property-matched, structurally-distinct molecules from elsewhere), just not the same claim a true DUD-E/ZINC decoy set would make.

The compound and every decoy (default 50, adjustable 5–200) are docked with identical settings. The result is a **percentile rank** (what fraction of decoys this compound's score beats; ties count as half) and a discrimination label:

| Percentile | Label |
|---|---|
| ≥ 90 | Strong |
| ≥ 65 | Moderate |
| < 65 | Weak |

This is expensive (~n_decoys + 1 real Vina runs, so ~2–5 minutes for the default 50) and deliberately on-demand per compound, not run automatically for a whole batch. Decoys are docked concurrently (a thread pool, one Vina subprocess per worker, each capped to a fraction of the machine's cores via `--cpu`) rather than sequentially — the sequential version was measured at 10+ minutes wall-clock with 31 of 32 cores idle the whole time, since a single Vina run at default exhaustiveness only keeps ~8 cores busy.

### 2.11 Blind docking

An alternative to site-specific docking: instead of a small box centered on a known or assumed pocket, the search box spans the entire protein (computed the same way as a site-specific box, just from every atom in the receptor rather than just the ligand's neighborhood, with smaller padding — 4 Å vs. 8 Å, since the box is already large). This is offered whenever a receptor exists on disk, even for a target with no known ligand at all.

**Explicit, surfaced caveat** (`app.py::BLIND_CAVEAT`, shown in the UI whenever Blind mode is active): searching the entire protein surface is slower and substantially less reliable per-site than a validated pocket search. Any resulting hit should be treated as a candidate site to investigate, not a confirmed pose or affinity. Users are told to consider raising exhaustiveness well above the default when running a real blind search.

---

## 3. Architecture — module map

### Backend (`backend/docking/`)

| Module | Responsibility |
|---|---|
| `availability.py` | Detects which binaries/packages are present (Vina, GNINA, fpocket; rdkit, meeko, posebusters, pdbfixer, matplotlib) and whether the feature is usable at all. |
| `ligand_prep.py` | SMILES → 3D conformer (RDKit ETKDG, retried across 4 fixed seeds on failure, then a random-coordinates fallback) → MMFF optimization → PDBQT (Meeko). |
| `validity.py` | The PoseBusters gate (§2.3). |
| `rmsd.py` | Atom-map-safe RMSD (§2.8). |
| `consensus.py` | Best-pose selection + confidence tiering (§2.4). |
| `engines.py` | `VinaEngine` (the docking engine wrapper) and `GninaRescorer`/`NullRescorer` (the optional CNN second opinion). |
| `profile.py` | Registry (`docking_registry.json`) load/save, portable receptor-path resolution, cross-process file locking for concurrent registry writes, atomic (temp-file + rename) writes so a reader never sees a half-written file. |
| `receptor_prep.py` | The full one-time receptor-preparation pipeline (§2.9) — extraction, stripping, repair, PDBQT conversion, binding-site computation, the manual-override on-demand path. |
| `pipeline.py` | Per-compound orchestration (`dock_compound`) and per-target validation (`redock_reference`, `redock_reference_for_profile`). |
| `interactions.py` | The built-in, dependency-light distance-based interaction fallback and its own standalone 2D diagram renderer (superseded in practice by `interaction_diagram.py` below, which is what the live pipeline actually calls). |
| `interaction_diagram.py` | PLIP detection (if installed) + the LigPlot+-style 2D diagram renderer actually used by the pipeline; both the built-in PNG and downloadable SVG/TIFF/PDF export live here. |
| `enrichment.py` | Fresh Decoy Validation (§2.10). |
| `recommend.py` | Disease → ranked target → "why this structure" evidence bundle (crystallographic context from `panel_results_v2.csv`, not proof of docking accuracy). |
| `failure_diagnostics.py` | Classifies a failed result into a small set of concrete categories (invalid molecule, conformer generation failed, ligand conversion failed, engine unavailable, no valid pose, etc.) with a suggested action — built by pattern-matching real error strings the pipeline actually produces, not inventing categories nothing here triggers. |

### Supporting scripts (`backend/scripts/`), used by the docking pipeline

| Script | Role |
|---|---|
| `select_receptor.py` | Automated candidate-structure proposal for onboarding a *new* target (UniProt → best-resolution liganded PDB entry → best ligand). Also owns the canonical `ADDITIVE_BLACKLIST` and `MIN_LIGAND_MW` every other module imports. Admin/maintainer tool, not part of the interactive user journey. |
| `pdb_fetch.py` | Downloads a PDB entry from RCSB; falls back to mmCIF→legacy-PDB conversion (via `gemmi`) for depositions with no legacy `.pdb` file, preserving a name map for the rare case a 4–5 character modern ligand code gets truncated to fit legacy PDB's 3-character field. |
| `detect_chain.py` | Picks which protein chain to restrict receptor prep to, when a structure has multiple copies of the same protein — whichever chain actually carries the reference ligand. |
| `validate_target.py` | The full new-target validation pipeline (admin tool): fetch → prepare → extract crystal pose → redock → RMSD gate → optional enrichment test → registry update, with an accept-or-revert guard so a worse candidate can never silently overwrite an already-validated (better) entry. |
| `generate_decoys.py` | The DUD-E-style decoy generation shared by both the new-target enrichment test and the live Fresh Decoy Validation feature (§2.10). |

### Other backend modules that participate in a docking run

- `export_package.py` — builds the downloadable ZIP (§5.6).
- `research_report.py` — assembles the multi-stage "evidence chain" report (§7, dead-in-UI but live via API).
- `serving/screen.py` — the combined QSAR+docking Screen pipeline reuses `docking/pipeline.py` and `docking/profile.py` directly; not documented in full here (Screen has its own documentation scope) but shares the exact same docking core.

### Frontend (`frontend/src/`)

| File | Responsibility |
|---|---|
| `tabs/DockingTab.tsx` | The tab itself: readiness gate, sidebar (target/mode/SMILES/plant-source/advanced settings), submission, polling, and the results table. |
| `lib/useAdvancedDocking.ts` | All Advanced Settings + binding-site state for one target selection: automatic site loading, manual structure candidates, ligand picking, residue selection → box, drag-to-edit box, exhaustiveness/poses/GNINA overrides. Shared verbatim between Docking and Screen. |
| `components/TargetBrowser.tsx` | Disease-first target search (disease combobox → ranked targets → download gating → structure preview), shared by Docking and Screen. |
| `components/ManualStructurePicker.tsx` | The structure-override UI: candidate list, co-crystallized-ligand list, before/after receptor comparison, prep-report download. |
| `components/DockingModeSection.tsx` | Site-specific/Blind toggle + the current box summary + "View/edit binding site" entry point. |
| `components/BindingSiteModal.tsx` | The full 3D binding-site editor: residue checklist, drag-to-move/resize box handles, live 3D preview. |
| `components/AdvancedSettingsPanel.tsx` | Exhaustiveness, pose count, GNINA toggle, box-source explanation, reset-to-automatic. |
| `components/DockingPieces.tsx` | Shared results-rendering pieces: the per-compound detail panel, interaction table, Fresh Decoy button, complex-download button, and two currently-unused-in-UI components (`ResearchReportButton`, `AlternateLigandButton` — §7). |
| `components/RedockingValidationNote.tsx` | The automatic redocking-validation summary shown above the results table. |
| `components/PoseOverlayViewer.tsx` | The experimental-vs-redocked 3D pose overlay (part of redocking validation). |
| `components/PoseViewer.tsx` | The per-compound 3D pose modal (ligand + receptor context). |
| `components/InteractionLogTable.tsx` | The cross-compound consolidated interaction log. |
| `components/ReceptorPreview.tsx` / `ReceptorBeforeAfter.tsx` | Read-only 3D structure preview, and the before/after receptor-prep comparison view. |
| `components/DownloadGateBar.tsx` | Progress bar for on-demand model/docking-data downloads. |

---

## 4. The user journey, in order

### 4.1 Readiness gate

The Docking tab checks `/api/docking/status` on load. If the environment doesn't have what's needed (Vina on PATH; RDKit, Meeko, PoseBusters installed), the tab shows **"Docking — not yet enabled on this machine"** with a checklist of exactly which Python packages and engine binaries are present or missing, plus the planned pipeline description — never a bare "unavailable" message. GNINA (`has_gnina`) is tracked separately as optional; its absence does not block readiness, it only caps confidence (§2.4). Once ready, the tab activates automatically — there is no manual "enable" step.

### 4.2 Target selection (sidebar)

**Disease (optional)** — a searchable combobox over diseases from the curated panel CSV. Picking one narrows the Target list below to that disease's associated targets, ranked by disease-association score, and jumps focus straight to the target box. This is entirely optional; a user who already knows the target id can skip straight to the Target field.

**Target** — a searchable combobox with three visually distinct row types:
- ✓ **modeled** (brand color) — has a trained QSAR model.
- ⚙ **docking-only** (muted) — no QSAR model exists for this protein, but structure-based docking is still available (routed via a synthetic `GENE_<symbol>` id, not a real ChEMBL-prefixed target id).
- ⬇ **not-downloaded** — a target whose data bucket hasn't been fetched to this machine yet; picking it triggers the download gate.

Typing narrows across three pools simultaneously: already-downloaded targets, downloadable-but-not-yet-fetched targets, and every other registry target reachable by symbol or id (mostly `GENE_` docking-only entries) — a target not tied to any disease association is still findable this way. Only the first 30 matches are shown per query (`TARGET_PAGE`).

**Download gate**: if the picked target's data isn't on disk, a progress bar appears (`DownloadGateBar`) with live byte progress and a Stop button — the target only becomes usable once the download finishes.

**Manual structure (overrides the automatic recommendation)** — appears the instant a target is picked and its download (if any) completes. Lists every qualifying PDB structure on record for this target's gene (from the curated crystallographic-quality CSV), each annotated with its CSV quality rank, resolution, ligand-model quality (RSCC/RSR), and its suggested ligand — or, for a structure with no validated ligand at all, a `blind only` badge. The current automatic default is marked. Clicking a normal candidate lists its real co-crystallized ligands (fetched live from RCSB via `/api/docking/pdb_ligands`, an async job since it may need to download the structure first); **no ligand is pre-selected** — which one is biologically relevant is a judgment call the structure's biology determines, not a generic heuristic, so the user always picks explicitly. Clicking a blind-only candidate builds its receptor directly (there's nothing to choose). Picking a ligand triggers the full on-demand receptor-preparation pipeline (§2.9), with live step-by-step status text (`"Extracting reference ligand…"`, `"Repairing (PDBFixer)…"`, etc.) so a ~1–2 minute prep isn't one static, unmoving message. Manually-built structures are **never persisted to the shared registry** — they're a per-request profile the caller carries forward explicitly (`custom_profile`), so an individual's exploratory pick can never silently override the vetted default for everyone else, and can't race a concurrent batch validation run's registry writes. They *are* cached to disk per `(target_id, pdb_id, ligand_resname)` so re-picking the same combination doesn't rebuild from scratch — deliberately, since PDBFixer's hydrogen placement is not perfectly deterministic between independent runs, which would otherwise make redocking-validation RMSD drift run to run for what should be an identical pick.

Once a manual structure is built: **"Compare before/after receptor prep"** opens a side-by-side 3D view (original deposited structure with waters/heteroatoms visible, vs. the stripped/repaired/protonated version actually used for docking). **"Download prep report"** exports the real step-by-step before/after numbers (atom counts at each stage, pocket residue count, box dimensions) as a plain-text file, generated client-side from data already in the response — no extra server round-trip. If the pick triggered a sanity-check warning (unusually small/large ligand, too few pocket contacts), it's shown directly beneath the picker in an amber-bordered box, never suppressed.

Below the structure picker, once a target is active: a **"View 3D structure"** link opens a read-only 3D preview of whichever structure is actually in play (the automatic default, or the manual pick), labeled with its PDB id.

### 4.3 Docking mode

**Site-specific / Blind** segmented toggle.

- **Site-specific** (default): the search box is centered on the resolved binding site. The summary line states, in real numbers, how many residues are selected out of how many exist in the receptor, and the active box dimensions.
- **Blind**: the search box spans the entire protein. The summary states the whole-protein box's real dimensions. If no prepared receptor exists on disk for this target at all, this mode is unavailable and says so plainly rather than showing a guess.

**"View / edit binding site"** opens the full 3D editor (`BindingSiteModal`) whenever any box (site-specific or blind) is defined:

- A live 3D view of the receptor (cartoon by default; style switchable — cartoon/surface/other options via `StyleSelect`) with the currently-selected pocket residues highlighted as yellow sticks, and the search box drawn as a wireframe.
- **Every residue in the receptor is listed** (not just the automatically-detected pocket), with the automatic pocket pre-checked — filterable by name/number/chain once the list exceeds 15 entries. Checking or unchecking any residue recomputes the box from the checked set (debounced 400ms), via `/api/docking/box_from_residues`.
- **"Enable drag editing"** (site-specific mode only — the whole-protein box in Blind mode is fixed and explicitly not editable) turns on six colored resize handles (one pair per axis, anchored on the opposite face when resized) plus one orange move handle, positioned live via the 3D viewer's own screen-projection so they track the camera as it rotates. Dragging a handle recomputes center/size in real time and immediately redraws the box.
- **Precedence rule, and why it exists:** a manual drag "wins" over a later residue-checkbox change and stays put until the user explicitly resets — otherwise, toggling any single residue after a careful drag would silently discard that drag and recompute the box from the checked residues instead, which was an identified, real bug in an earlier version. The current box's exact source (residue-derived, drag-adjusted, or fully automatic) is stated in Advanced Settings' "Binding box" field.
- A **Current box** readout shows live center/size numbers regardless of source.

### 4.4 SMILES input and plant source

**SMILES (one per line)** — a plain textarea; each non-empty line after trimming becomes one compound submitted for docking.

**Plant source (optional)** — a free-text field ("e.g. *Curcuma longa*") purely for provenance: it's threaded through into the exported metadata so a batch of compounds stays traceable back to whatever natural source it came from, and is later reused (if provided) as the default literature-search query in the Research Report feature (§7). It has no effect on docking itself.

### 4.5 Advanced Settings (collapsed by default; auto-opens once a target is picked)

- **Exhaustiveness** (1–64, default 8 if left blank) — Vina's search-thoroughness parameter; higher values search more but take longer.
- **Number of poses** (1–20, default 9 if left blank) — how many candidate binding modes Vina returns per compound.
- **GNINA CNN rescoring** checkbox (on by default when GNINA is installed) — the only way to force the optional CNN second opinion off; unchecking it substitutes a `NullRescorer` that always reports "unavailable," which correctly caps confidence at `medium` per §2.4's table.
- **Binding box** — a read-only status line stating exactly why the current box is what it is (automatic / drag-adjusted / residue-derived / not yet defined).
- **Reset to Automatic** — discards every override (exhaustiveness, pose count, GNINA toggle, manual structure, box override) and reloads the target's registry default from scratch.

Every field here defaults to "Automatic" and is opt-in — nothing here is required, and nothing here is persisted between sessions; it exists purely for the current request.

### 4.6 Submission and background execution

Clicking **"Dock compounds"** validates first (a target must be picked; at least one non-empty SMILES line; a gene-only target with no automatic default must have a manual structure picked in Advanced Settings; site-specific mode must have *some* box available) and then submits. While running, the button reads **"Preparing structure…"** if a manual structure build is still in flight, otherwise it's simply disabled during submission/polling.

Server-side (`_run_docking_job`, a background thread), for each compound in order:
1. `ligand_prep.prepare_ligand` — SMILES → 3D conformer. Embedding is retried across four fixed seeds (`0xf00d, 1, 42, 7`) before falling back to a random-coordinates embed, specifically to rescue strained or unusual ring systems (natural products, the app's primary domain, are disproportionately likely to need this).
2. `VinaEngine.dock` — the actual search.
3. PoseBusters gate → best valid pose + self-consistency.
4. GNINA rescoring (if enabled and available).
5. Interaction detection (PLIP or fallback) + 2D diagram rendering.
6. Confidence assignment.
7. Failure classification, if applicable (§2.13/`failure_diagnostics.py`) — every non-`ok` result gets a category and a concrete suggested action, not just a raw error string.

After every compound in the batch finishes: **redocking validation runs automatically once**, against the same structure and parameters (§2.7) — this never blocks or fails the compound results even if it itself errors.

The UI polls every 2 seconds and shows **"Docking N/total… (minutes per compound)"** with a **Stop** button that requests cancellation — cancellation is checked *between* compounds (not mid-compound), so a stopped run still returns whatever finished before the stop request landed, never discarding completed work.

### 4.7 Results — the table

One row per submitted compound:

| Column | Content |
|---|---|
| (expander) | ▸/▾ — only present when there's something to show (an interaction diagram, an `ok` status, or a suggested action for a failure) |
| Compound | the SMILES, truncated with ellipsis, full string on hover |
| Confidence | a colored dot (brand=high, amber=medium, clay=low, gray=none) + the label |
| Vina (kcal/mol) | the best pose's raw Vina score, or `—` |
| GNINA CNN / GNINA affinity / GNINA (kcal/mol) | only shown as columns at all if *any* row in the batch has GNINA data — an empty column set for a run with GNINA disabled |
| Status | `"{n} valid pose(s)"` on success, or the failure reason |
| Fresh decoy check | the **"Run Fresh Decoy Validation"** button (§2.10), only offered for compounds that produced a usable pose |

Clicking a row (when it has detail to show) expands the **per-compound detail panel** below it.

**Redocking validation note** (above the table, once): states the RMSD number, validated/not-validated verdict, and offers **"Experimental vs. redocked pose overlay"** — a 3D view superimposing the real experimental (crystal) pose (gray sticks) against the freshly redocked pose (magenta sticks) of the *same* reference ligand only, with Reset-view and spin controls. This view is never offered for an arbitrary submitted compound — only for the specific reference ligand the validation used, since the overlay (and the RMSD it illustrates) only means something when both poses are of the identical molecule.

### 4.8 Results — the per-compound detail panel

Two-column layout:

**Left — "Binding site & interactions" card**, stretching to fill the row's height:
- If the compound's residue-contact overlap with the reference ligand was computed, and/or which interaction detector ran, both are stated as a small text line above the image.
- The 2D interaction diagram (§2.6), click-to-enlarge into a lightbox modal with the same image at full size plus **SVG** and **TIFF** download buttons (vector/high-resolution export, regenerated server-side on demand from the job's already-stored SMILES and interaction data — not cached per-format).
- If no diagram exists for this pose, a plain "No interaction diagram for this pose" message, not a broken image.

**Right column, stacked cards:**
- **Summary** — a compact stat-tile grid: Vina score, confidence, GNINA CNN score/affinity/Vina-style affinity (whichever are present), PoseBusters-valid pose count, pose self-consistency count, and — for a non-`ok` result — the status and reason, tinted amber.
- **3D pose** — **"View 3D pose"** opens a modal with the docked ligand rendered as ball-and-stick (all atoms including hydrogens — nonpolar hydrogen positions are geometry-completed for display, explicitly noted as an estimate rather than Vina's own optimized placement, §2 pipeline note in `pipeline.py::pose_pdb_with_hydrogens`), the receptor shown around it in a switchable style (cartoon/surface/etc.) with the actual contacting pocket residues highlighted, Reset-view and spin controls. Beside it, **"Download complex (PDB)"** (or "Download pose (PDB)" if no receptor is available) exports the combined receptor+ligand structure as a real file, in one horizontal row rather than stacked.
- **Protein-ligand nonbonding interactions** table — every detected contact (name, category, type, distance in Å), sorted by distance, in its own scrollable card.

### 4.9 Export

At the bottom of the results table, once a job has finished:
- **Interaction table (N interactions)** toggle — expands the cross-compound consolidated interaction log below the table: one row per interaction *event* (not aggregated by residue or compound), with a **Download CSV** button.
- **Failure log (.csv)** — every non-`ok` compound with its status/category/reason/suggested action, in one small CSV (a placeholder row if nothing failed).
- **Full experiment package (.zip)** — see §5.6 below.

---

## 5. Background automatic work and decisions — a consolidated reference

This section exists because the user-facing journey above deliberately hides most of this; here it is made explicit.

### 5.1 What happens automatically, with no user action

- Redocking validation, on every submission with a real reference ligand (§2.7).
- Confidence tiering, from validity + self-consistency + (if available) GNINA (§2.4).
- Interaction-detector selection (PLIP if installed, else the built-in fallback) — never user-chosen.
- Failure categorization for any non-`ok` result.
- The reference-ligand pick for a target's *automatic* default structure (largest non-additive HETATM group) — but never for a manually-picked structure, where the user always chooses the ligand explicitly (§4.2).
- The default docking box (ligand-centered, 8 Å padding, 20 Å minimum per side).
- The blind (whole-protein) box, computed for every structure with a receptor on disk regardless of whether it's ever used.
- Chain selection for a multi-copy structure (whichever chain actually carries the reference ligand).
- Ligand-embedding seed retries in `ligand_prep` (four fixed seeds, then a random fallback) — invisible unless it fails entirely.
- Reproducibility metadata capture (`run_metadata`) at submit time — target, box, structure source, engine settings, software versions, fixed embedding seeds — regardless of whether the user ever looks at it (it's what powers `/reproduce`, §7, and the export package's `metadata.json`).

### 5.2 What is automatic *by default* but user-overridable

- Docking mode (site-specific vs. blind).
- Which structure is used (Advanced Settings' manual override).
- Which co-crystallized ligand a manually-picked structure centers on.
- The binding box itself (residue selection or 3D drag).
- Exhaustiveness, pose count.
- Whether GNINA rescoring runs.

### 5.3 What is never automatic — always an explicit user choice

- Which ligand a manually-picked PDB structure centers on (§4.2) — deliberately no "main ligand" heuristic here, unlike the fully-automatic default path.
- Running Fresh Decoy Validation for any given compound (always opt-in per compound, never batch-triggered, because it is expensive).
- Downloading the full export package or the failure log.

### 5.4 The registry — `docking_registry.json`

Every target's automatic-default profile (receptor paths, box, reference ligand, validation status) lives in one shared JSON file, keyed by `target_id`. Notable engineering details:

- **Portable paths.** Receptor file paths are stored as basenames and resolved at load time relative to `docking_targets/<target_id>/` — not absolute paths baked in at prep time, which used to break the instant the project folder moved. A legacy absolute path is still honored if it happens to still exist, but the portable form is always preferred.
- **Cross-process locking.** A `registry_lock()` context manager (an `fcntl.flock` on a sibling `.lock` file) wraps every read-modify-write cycle, so two processes validating different targets concurrently (the batch validation script's parallel runner) can never have one's write silently clobber the other's — the second process blocks at its read until the first's write completes.
- **Atomic writes.** Every write goes to a same-directory temp file first, then `os.replace`s it over the real registry file — a concurrent reader (the live app, mid-request) can only ever see a fully-complete old or new file, never a half-written one.

### 5.5 How a new target gets onboarded (maintainer-side, not part of the interactive user journey)

This is background/administrative process, included here for completeness since it's part of the same scientific pipeline:

1. `select_receptor.py` proposes candidate (PDB id, ligand) pairs for a gene: UniProt lookup → best-resolution X-ray structures with a bound ligand (via UniProt's own cross-reference list, not a fragile RCSB schema query) → best non-additive ligand by formula weight.
2. `validate_target.py` runs the real pipeline against a chosen candidate: fetch → `receptor_prep.prepare_receptor` → fetch the ligand's real SMILES from RCSB's Chemical Component Dictionary → build a bond-order-correct crystal pose (`AssignBondOrdersFromTemplate`, since RDKit's PDB-coordinate-only bond perception frequently gets bond *orders* wrong even when atom connectivity is right) → **reference redocking** (the RMSD gate that actually flips `"validated": true`) → an optional enrichment test (this target's own known most-potent vs. least-potent compounds, or real DUD-E-style property-matched decoys — recorded either way as `enrichment_source` so nobody mistakes a lighter-weight proxy for the real thing; this step never gates `"validated"`, only RMSD does).
3. An accept-or-revert guard: if a new candidate's RMSD is *worse* than an already-validated entry's, both the registry entry **and** the on-disk receptor files (which the build step already overwrote unconditionally before this check could run) are restored to the prior, better entry — `--force` is required to intentionally accept a regression.

### 5.6 Export package contents (`export_package.py`)

One ZIP, built fresh on request (not cached), containing:
- `metadata.json` — the full reproducibility snapshot (§5.1).
- `receptor.pdb` — the exact receptor file used.
- `results.csv` — one row per compound: SMILES, status, Vina score, confidence, valid-pose count, failure category/reason/action, GNINA CNN score/affinity/Vina-affinity, Fresh Decoy percentile/discrimination/compound score (if run).
- `interactions.csv` — the same cross-compound interaction log the UI's "Interaction table" shows, one row per interaction event.
- `poses/<compound>.pdb` — the receptor+pose combined structure, per compound that produced a pose.
- `interactions/<compound>.png` — the 2D interaction diagram, per compound that has one.
- `fresh_decoy/<compound>.json` + `<compound>_decoys.csv` — only for compounds Fresh Decoy Validation was actually run for.
- `redocking_validation/summary.json` + `redocked_pose.pdb` + `experimental_pose.sdf` — the target-level validation result, when one exists for this run.

---

## 6. Limitations — stated directly, not buried

- **Docking scores are not binding affinities.** A Vina score is a relative, physics-approximation-based ranking signal, not a measured or even reliably-predicted ΔG. Two compounds with similar scores are not necessarily similarly potent.
- **Confidence tiers are heuristic, not calibrated probabilities.** "High" confidence means three specific, independent signals agree — it is not a statistically calibrated probability of correctness, and the code documents this distinction explicitly (§2.4).
- **The built-in interaction detector, when PLIP isn't installed, is coarse.** It distinguishes only H-bonds and generic hydrophobic contacts, with no atom-level detail — every such result is labeled `"(distance-based)"` so this limitation is visible in the output itself, not just in documentation.
- **Blind docking is markedly less reliable per-site than site-specific docking**, by construction (a much larger search volume for the same exhaustiveness) — stated as an explicit caveat every time Blind mode is used, not a silent caveat readers might miss.
- **Fresh Decoy Validation's "decoys" are not confirmed non-binders.** They are real compounds curated for other targets' assays, property-matched and structurally dissimilar to the query — a real, useful, but different comparison than a true experimentally-confirmed inactive set.
- **Redocking validation checks the docking *setup*, not any submitted compound.** A "validated" badge says Vina can reproduce a known-correct pose on this exact receptor/box/parameters; it says nothing about whether any of the user's own submitted compounds dock correctly.
- **Receptor preparation is a real engineering pipeline with real failure modes**, most already found and fixed (documented in `receptor_prep.py`'s own comments): PDBFixer-introduced clashes crashing downstream bond perception, OXT-atom valence errors at chain termini, multi-copy structures producing physically wrong receptors if not restricted to one chain, ligand HETATM records sometimes filed under a chain that isn't the one the ligand actually contacts. Each fix is real and tested against a specific confirmed case, not speculative.
- **PDBFixer's hydrogen placement is not perfectly deterministic** between independent builds of the same structure — this is why manually-built structures are cached rather than rebuilt on every pick (otherwise redocking-validation RMSD could drift between repeated runs of what should be an identical setup).
- **Module status labels in `docking/__init__.py` are stale** (they still describe an AutoDock4-era architecture) — the authoritative status is the README (`docking/README.md`) and this document, not that file's own docstring.
- **This pipeline was validated end-to-end on specific real structures** (e.g. `CHEMBL1862_ABL1`/1IEP) but large parts of `receptor_prep.py` and `engines.py` are marked in their own docstrings as written-but-not-independently-executed in the original development environment (PDBFixer/real-PDB testing was unavailable there) — the actual validation discipline is: never trust a target's docking output until its own redocking-validation RMSD is below threshold, which the UI always shows.

---

## 7. Features present in the code but not currently reachable from the UI

Documented explicitly so nobody mistakes "not in the UI" for "doesn't exist" or accidentally reintroduces it as if it were new:

- **"Reproduce this analysis"** (`POST /api/docking/job/{jid}/reproduce`) — resubmits a finished job with its exact saved parameters (the literal receptor file, box, exhaustiveness, GNINA setting, compound list) rather than whatever the registry/UI currently defaults to, so a result stays reproducible even if a target's default structure changes later. The backend endpoint is fully functional; the UI button that called it was removed at a later product-design pass (button-hierarchy cleanup) and has not been reintroduced.
- **"Dock again with a different ligand"** (`POST /api/docking/job/{jid}/alternate_ligand/build` + `/submit`) — lets a user re-center a finished job's structure on a *different* real co-crystallized ligand from the same raw PDB, keeping every other setting identical. Also fully functional server-side and in a still-defined-but-unused frontend component (`AlternateLigandButton` in `DockingPieces.tsx`); not currently rendered in `DockingTab.tsx`.
- **"Generate research report"** (`POST /api/docking/job/{jid}/research_report`) — assembles the full multi-stage evidence-chain report (natural source → chemical identity → reported literature activity → target prediction → QSAR → docking → interactions → ADMET → off-target analysis → a templated, non-fabricated plain-language summary and a Methods draft) for one already-docked compound. Fully functional; the frontend component (`ResearchReportButton`) is defined but was explicitly removed from the single-compound detail view per a later product decision, and not currently rendered anywhere in the Docking tab.

All three remain fully callable via the API directly, and their removal from the UI was a deliberate scope/clutter decision at a specific point in this project's history, not a sign of a broken or abandoned feature.

---

## 8. API reference (docking-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `GET /api/docking/status` | Readiness (tools/packages present), plus every registry target's id/name/site_source. |
| `GET /api/docking/targets` | Just the registry's target id/name list. |
| `GET /api/diseases`, `GET /api/diseases/{id}/targets` | Disease browsing (backs the disease combobox). |
| `GET /api/targets/{target_id}/recommendation` | The "why this structure" evidence bundle. |
| `GET /api/docking/receptor/{target_id}` | Raw receptor PDB text (registry default). |
| `GET /api/docking/receptor_file?path=` | Raw receptor PDB by explicit path (custom/manual structures) — path-restricted to `docking_targets/`. |
| `GET /api/targets/{target_id}/structure_candidates`, `GET /api/genes/{symbol}/structure_candidates` | Candidate PDB structures for the Manual Structure picker. |
| `GET /api/targets/{target_id}/binding_site` | Automatic-default binding-site evidence (pocket residues, box, all residues). |
| `GET /api/targets/{target_id}/reference_ligand.sdf` | The crystal reference pose, when one was built via the full validation pipeline. |
| `POST /api/docking/box_from_residues` | Recompute a box from a picked residue set. |
| `POST /api/docking/pdb_ligands` + `GET .../job/{jid}` | List every real co-crystallized ligand in a PDB entry (async job). |
| `POST /api/docking/receptor/custom` + `GET .../job/{jid}` | Build an on-demand receptor for a manually-chosen (PDB, ligand) pair (async job). |
| `POST /api/docking/submit` | Submit a docking batch. |
| `POST /api/docking/job/{jid}/reproduce` | Resubmit with a past job's exact parameters (§7 — API-only). |
| `GET /api/docking/job/{jid}/alternate_ligands`, `POST .../alternate_ligand/build`, `POST .../alternate_ligand/submit` | Re-center and redock against a different ligand from the same structure (§7 — API-only). |
| `POST /api/docking/cancel/{jid}` | Request cancellation of an in-progress batch. |
| `GET /api/docking/job/{jid}` | Poll job status/results. |
| `GET /api/docking/job/{jid}/export_package` | The full ZIP (§5.6). |
| `GET /api/docking/job/{jid}/interaction_diagram?smiles=&fmt=` | On-demand SVG/TIFF/PNG/PDF export of one compound's interaction diagram. |
| `GET /api/docking/job/{jid}/failure_log` | The failure-only CSV. |
| `POST /api/docking/enrichment/fresh` + `POST .../cancel/{jid}` + `GET .../job/{jid}` | Fresh Decoy Validation (§2.10). |
| `POST /api/docking/job/{jid}/research_report` | The evidence-chain report (§7 — API-only). |

---

## 9. Document provenance

Written by reading, in full: every file under `backend/docking/`, the docking-related routes in `backend/app.py`, `backend/export_package.py`, `backend/research_report.py`, the supporting scripts (`select_receptor.py`, `pdb_fetch.py`, `detect_chain.py`, `validate_target.py`, `generate_decoys.py`), and every frontend component reachable from `frontend/src/tabs/DockingTab.tsx` (`useAdvancedDocking.ts`, `TargetBrowser.tsx`, `ManualStructurePicker.tsx`, `DockingModeSection.tsx`, `BindingSiteModal.tsx`, `AdvancedSettingsPanel.tsx`, `DockingPieces.tsx`, `RedockingValidationNote.tsx`, `PoseOverlayViewer.tsx`, `PoseViewer.tsx`, `InteractionLogTable.tsx`, `ReceptorPreview.tsx`, `ReceptorBeforeAfter.tsx`, `DownloadGateBar.tsx`). No content here was reconstructed from memory of past conversation — every claim traces to a specific line of code read during this pass. Where the code itself documents a design decision or a fixed bug, that reasoning is reproduced here rather than paraphrased away.
