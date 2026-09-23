# Phase 2 Index Rebuild — Execution Plan

**Scope of this document**: how to actually *run* the rebuild decided in `BUILD_PLAN.md` Phase 2 items 1-7. Nothing here has been executed yet — this is the plan, staged and estimated, for review before committing real time/API/disk resources. Every stage below implements a specific, already-made decision; nothing new is decided here.

---

## Pre-flight checks (done, informing this plan)

| Check | Result |
|---|---|
| Disk space | `/home/storage`: 60GB free of 1.9TB (97% used). v1's *finished* index is 56MB — tiny. Raw intermediate pulls are the real space cost, estimated 3-5GB combined (see Stage sizing below). Feasible with margin, but should not be left uncompressed/undeleted after each stage. |
| `standard_relation__in` query support | Confirmed working, but combining it with `pchembl_value__isnull=false` is a no-op (that filter already implies `relation='='`, per L4) — censored/Potency/inactive pulls must filter on `standard_value__isnull=false` instead and compute pChEMBL locally from raw values, per items 1/3's decision. Three/four genuinely separate pulls, not one combined query. |
| ChEMBL release pinning | `status.json` → `chembl_db_version` + `chembl_release_date`, confirmed live (item 7). Call this **once, first**, before any pull starts, and record it — this is the closest this REST-based approach gets to a snapshot version (see Known Limitation below). |
| Single-protein-target allowlist | Not preserved from v1's build (confirmed absent from the freeze snapshot, only referenced in the build script's docstring) — rebuilt fresh: **5,869 ids**, exactly matching v1's own historical count, a reassuring consistency check that the human single-protein-target universe hasn't drifted. |
| **Sequential fetch throughput (MEASURED, not estimated)** | A first pilot (single-connection sequential pagination) measured only **~54 rec/s steady-state** — 5x worse than this plan's original v1-history-based ~2.5-3hr estimate, extrapolating to **~13.3 hours for Stage 1 alone**. Root cause: each individual page request takes ~13-18s round-trip regardless of page size — network/server latency-bound, not this machine's bandwidth. |
| **Concurrency (MEASURED, fixes the above)** | A live microbenchmark found near-linear scaling up to ~10 concurrent connections and continued sub-linear gains to 30-40, with **zero errors or rate-limit responses at any tested concurrency**. Chose 30 workers. **Validated at full production scale** (100,000-record pilot, not a small benchmark): **~885 rec/s sustained**, cutting Stage 1's full-run estimate to **~47 minutes** (0.79hr) — all Stage time estimates below are revised accordingly, superseding the original sequential-based figures. Data integrity checked on the 100K-record pilot output: 0 malformed rows, 88 exact-duplicate rows (0.088%, a minor concurrent-pagination-boundary artifact, harmless — Stage 7's aggregation absorbs exact duplicates naturally). |

**Known limitation, disclosed not hidden**: this plan pulls via ChEMBL's live REST API over an estimated (revised, see above) 2-4 hours total across all stages, not the originally estimated 6-10+. ChEMBL's underlying database does not change mid-release, so this is not expected to cause real drift — but it is not a true atomic snapshot the way a single bulk database dump download would be. Pin the release number at the start (above) and treat that as the version of record; do not re-verify it obsessively mid-pull.

---

## Stages

### Stage 0 — Pin release + verify allowlist (~1 min, MEASURED — actually run)
- `GET status.json` → record `chembl_db_version`, `chembl_release_date`. **Done: `ChEMBL_37`, released 2026-05-01.**
- Rebuild `single_protein_target_ids.json` fresh (not preserved from v1's build). **Done: 5,869 ids, ~51s — matches v1's historical count exactly.**
- Output: `phase2/data/release_manifest.json`, `phase2/data/single_protein_target_ids.json`. **Scripts written and run: `stage0_pin_release.py`.**

### Stage 1 — Base positive evidence, `=` relation (~47 min, MEASURED at production scale — 2.6M records, 30 concurrent workers)
- Same filter v1 already used (`target_organism=Homo+sapiens`, `standard_type__in=IC50,Ki,Kd,EC50`, `confidence_score__gte=8`, `pchembl_value__isnull=false`), bulk-paginated (limit=1000/page), **fetched with 30 concurrent connections** (see Pre-flight table — sequential pagination measured only ~54 rec/s; concurrency validated at ~885 rec/s at 100K-record production scale, 0 malformed rows, 0.088% harmless duplicate rows from concurrent-pagination boundaries).
- Reuses `pchembl_value` directly (safe here — this is exactly the case where ChEMBL populates it).
- Tag every row `is_censored=False`.
- **Script written and pilot-validated**: `stage1_fetch_base_evidence.py` — per-page-offset checkpointing (resumable, safe under concurrency), `--pilot N` flag for capped test runs. Full run (`python3 stage1_fetch_base_evidence.py`, no `--pilot`) not yet launched.

### Stage 2 — Censored evidence, `<`/`<=` (~3-5 min, estimated at Stage 1's measured rate — same query pattern, ~174K records)
- Same base filter, but `standard_relation__in=<,<=` and `standard_value__isnull=false` (NOT `pchembl_value__isnull=false` — confirmed no-op above). **Same 30-worker concurrent pattern as Stage 1**, not yet built as its own script (trivial adaptation of `stage1_fetch_base_evidence.py`'s query string).
- Compute pChEMBL locally from `standard_value`+`standard_units` (convert to molar, -log10) per §4's L4 correction.
- Tag every row `is_censored=True`. These are **lower bounds**, not exact values — carry that distinction through to Stage 7's aggregation (item 3's decision: don't silently treat as exact).

### Stage 3 — `Potency`-type records (~55-60 min raw pull, estimated at Stage 1's measured rate — ~2.98M records before filtering)
- `standard_type=Potency`, same organism/confidence filter, `standard_value__isnull=false`. **Same 30-worker concurrent pattern.**
- **Local filtering after pull** (item 3's decision, not done server-side): keep only `assay_type in {B,F}` and `standard_units='nM'` — expect meaningful shrinkage from the raw 2.98M once the non-quantitative/wrong-unit/wrong-assay-type records are dropped; the real usable count won't be known until this filter actually runs, worth reporting once it does.
- Tag `source_type='Potency'` for provenance (distinct from the four traditional bioactivity types, per the corrected framing that this isn't narrowly PubChem-sourced).

### Stage 4 — Measured-inactive evidence (~9-10 min, estimated at Stage 1's measured rate — ~490K records)
- `standard_relation__in=>,>=`, `standard_value__gte=10000` (nM, i.e. ≥10µM — symmetric with the ground-truth active boundary), `standard_units=nM`, same organism/confidence filter. **Same 30-worker concurrent pattern.**
- Tag `is_measured_inactive=True`.
- Secondary `activity_comment` pass (item 6's noted casing issue) is **optional/deferred** — the numeric-relation pull above is the primary source; don't block the main rebuild on resolving the casing inconsistency.

### Stage 5 — Orthologue tier (time TBD, design choice below)
- **Design decision for execution, not yet made in items 1-6**: pull **per major species in bulk** (one paginated pull per organism: *Mus musculus, Rattus norvegicus, Bos taurus, Sus scrofa, Oryctolagus cuniculus, Canis familiaris, Cavia porcellus, Macaca mulatta, Gallus gallus* — P0b-19's species list), each as a single big paginated query (`target_organism=<species>&standard_type__in=...&confidence_score__gte=8&standard_value__isnull=false`), **using the same 30-worker concurrent pattern as Stages 1-4** — **not** ~2,600 separate per-target REST calls, which the item-4 pilot's per-target `compounds_for_target()` approach would make impractically slow at this scale (confirmed slow in the item-4 pilot: each target needed its own paginated call).
- After each species pull, locally intersect against a pre-built gene-symbol → human-target map (built once via the confirmed `target_synonym__iexact` method from item 4/the adjudication pilot) to keep only records whose target maps onto one of our 4,658 human single-protein targets' orthologues.
- Tag `species_provenance=<organism>`, kept as a **separate stratum**, never merged into the native-human rows (item 4's decision).
- **Estimate genuinely uncertain** — depends on per-species activity volume, unmeasured at this scale. Recommend running this stage last and independently timed, not blocking the rest of the rebuild.

### Stage 6 — Standardization (CPU-bound, not network-bound; est. 30-60 min for ~6M rows)
- Swap `SaltRemover.StripMol()` → `rdMolStandardize.FragmentParent()` (item 5's confirmed fix) in the standardization step, reusing v1's fingerprinting/Murcko-scaffold code otherwise unchanged.
- Deduplicate by standardized SMILES as before; retain the **list** of contributing original `molecule_chembl_id`s per collapsed structure (item 5's provenance decision), not just one winner.

### Stage 7 — Aggregate to (compound, target) pairs (CPU-bound; est. 15-30 min)
- Reuse v1's Pass-3/Pass-4 aggregation pattern with these changes:
  - **Native-human tier** and **orthologue tier** aggregated as separate, never-merged output strata (item 4).
  - **Max pChEMBL per pair, not mean** (item 3's correction).
  - Per-pair fields retained: `n_documents`, `assay_types_seen`, `best_relation`, `confidence_score` (item 2), `pchembl_value` (items 1/3), `is_censored` (item 3), `is_measured_inactive` (item 6), `contributing_molecule_chembl_ids` (item 5), `species_provenance` (item 4, native-human tier omits this or sets `"Homo sapiens"`).

### Stage 8 — Freeze (item 7's protocol; ~10 min)
- `CHECKSUMS.sha256` over every output artifact.
- `manifest.json` with the Stage 0 release pin + all fetch filters + counts per stage.
- `V2_FREEZE.md` narrative doc, following `V1_FREEZE.md`'s exact structure (architecture, dataset table, frozen decisions, "what re-opens this freeze," regenerate/verify commands).
- Confirm `BUILD_PLAN.md` §2's findings table shows every row Closed or explicitly-still-open-with-reason.

---

## Total estimate and sequencing

| Stage | Est. time | Status | Bound by |
|---|---:|---|---|
| 0 | ~1 min | **DONE** (measured) | network |
| 1 | ~47 min | **Script ready, pilot-validated at production scale (885 rec/s)**, full run not launched | network, 30 concurrent workers |
| 2 | ~3-5 min | estimated (Stage 1's measured rate) | network, 30 concurrent workers |
| 3 | ~55-60 min | estimated (Stage 1's measured rate) | network, 30 concurrent workers |
| 4 | ~9-10 min | estimated (Stage 1's measured rate) | network, 30 concurrent workers |
| 5 | unmeasured, run independently | not started | network, 30 concurrent workers |
| 6 | ~30-60 min | not started | CPU |
| 7 | ~15-30 min | not started | CPU |
| 8 | ~10 min | not started | disk I/O |

**Revised sequential total: roughly 2.5-3.5 hours excluding Stage 5** — down from the original plan's 7-8+ hour estimate, entirely because of the concurrency fix found and validated in this pass (Stage 1's own sequential-vs-concurrent measurement was a 16x difference: ~13.3hr projected sequentially vs. ~47min measured concurrently at the same 2.6M-record scale). Stages 2-4 use estimated, not yet independently measured, throughput at Stage 1's validated rate (same query pattern and page size, so a reasonable extrapolation — should still be confirmed with their own pilots before the full run, not assumed).

**Parallelization — REVISED, measured, superseded the original plan**: the original text here recommended testing *between-stage* concurrency (Stage 1 + Stage 2 running as separate concurrent processes). That turned out to matter far less than *within-stage* concurrency: a single stage's own sequential pagination was the real bottleneck (~54 rec/s), and running 30 concurrent connections *within* Stage 1 alone gave a measured ~16x speedup (885 vs. 54 rec/s) with zero errors or rate-limit responses observed through 40 concurrent connections in a live microbenchmark. **Decision: give every stage its own 30 concurrent workers; run the stages themselves sequentially** (simpler to reason about, avoids untested compounding risk from running multiple 30-worker pools against ChEMBL at once, and the sequential total is now only ~2.5-3.5hr anyway — not worth the added complexity/risk of also parallelizing across stages for a further, unmeasured, likely-marginal gain).

**Checkpointing, required from the start** (the lesson from Phase 1's mechanism-support fetch, which crashed near the end and lost ~140 rows of unsaved progress before checkpointing was added) — **implemented and validated in Stage 1's script**: per-page-offset checkpointing (not a single cursor, since pages complete out of order under concurrency), safe under 30 concurrent writers via a single-writer lock, confirmed resumable. Every future stage's script follows the same pattern.

---

## What this plan does not do / current status

- **Stage 0: done.** **Stage 1: script written, pilot-validated at full production scale (100K records, 885 rec/s, 0 data-integrity issues) — the full ~2.6M-record run has not been launched yet.**
- Stages 2-4: not yet built as their own scripts (trivial adaptations of Stage 1's query string + tags), not yet run even as pilots.
- Does not resolve Stage 5's real runtime — flagged as genuinely unknown pending a first measurement.
- Does not decide whether to run the full Stage 1 (and later stages) now or wait — that decision belongs to whoever approves committing the ~47-minute Stage 1 run (and beyond).
