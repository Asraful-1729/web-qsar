# Target Prediction v2 — Data Rebuild Freeze

**Frozen:** 2026-09-22
**Status:** Locked. Do not change without a new validation cycle (see "What re-opens this freeze" below).
**ChEMBL release:** `ChEMBL_37`, released 2026-05-01 — the first time this program has pinned an exact release (a disclosed gap in every prior phase, `PHASE1_BENCHMARK_CARD.md` §1/§10).

This is the frozen output of `BUILD_PLAN.md` Phase 2's data rebuild — every decision in items 1-7 physically executed and verified, per `PHASE2_REBUILD_EXECUTION_PLAN.md`'s staged plan. `phase2/data/CHECKSUMS.sha256` lets anyone verify later whether the artifact files have drifted from what's described here.

## Dataset (frozen)

| | |
|---|---:|
| Distinct standardized compounds | 1,590,210 |
| Native-human (compound, target) pairs | 3,973,653 |
| Orthologue-tier (compound, human_target, species) pairs | 34,080 |
| Human single-protein targets | 5,869 |
| Human targets with ≥1 confirmed orthologue | 385 |
| Distinct non-human orthologue target ids | 513, across 9 species |

Full provenance: `phase2/data/manifest.json`.

## What changed from v1, and why (the actual decisions, physically verified)

1. **Potency as a graded index weight** (item 1): `pchembl_value` retained per pair (max, not mean — see item 3 below); the graded weight *function* is explicitly deferred to Phase 3A (`BUILD_PLAN.md` §5's own scope), not fixed here.
2. **Confidence as an index weight** (item 2): confidence 4-7 is not a continuum for single-protein targets (7/6 are protein complexes, 5/4 target-ambiguous — confirmed live against ChEMBL's `confidence_description` values). No fetch-filter relaxation; `confidence_score` retained per pair for a future two-level (8 vs. 9) weight. **Found during execution, not anticipated in planning**: `activity.json` does not serialize `confidence_score` in its own output at all (confirmed empirically — a full unrestricted record has no such key). Backfilled via a separate `assay.json` lookup, batched by `assay_chembl_id` (270,225 distinct ids resolved).
3. **Censored relations, `Potency` records, functional assays, most-potent aggregation** (item 3): ~174K censored (`<`/`<=`) records included with an `is_censored` flag (lower bound, not exact); ~2.96M `Potency`-type records included after quality filtering (`assay_type∈{B,F}`, `standard_units=nM`) — corrected from the plan's "PubChem" framing, predominantly Scientific Literature-sourced; functional assays already confirmed present in v1; **max pChEMBL per pair, not mean** — v1's own code was found to use mean, contradicting its own recipe text.
4. **Orthologue tier** (item 4): built via gene-symbol matching (`target_synonym__iexact`), not the weaker name-matched proxy Phase 0b used. **385 human targets have a confirmed orthologue — far fewer than Phase 0b's 2,600-candidate estimate**, which used a looser method. Kept as a fully separate stratum (`stage7_orthologue_pairs.jsonl`), never merged into native-human rows, granular per species (not collapsed across species).
5. **Salt/parent collapse** (item 5): `rdMolStandardize.FragmentParent()` replaces `SaltRemover.StripMol()`. **Found and fixed a live bug in v1's own shipped index**: `SaltRemover` cannot collapse n:1 stoichiometry salts (confirmed on the exact ladostigil tartrate case the Phase 1 adjudication study flagged as a "salt-form gap") — 2,167 of v1's 1,312,849 rows (~0.17%) carried this exact defect. `contributing_molecule_chembl_ids` retained per standardized structure (5,203 structures have >1 contributor, confirming the fix works at scale).
6. **Measured-inactives term** (item 6): ~490,821 genuinely tested-and-inactive records (`>`/`>=`, ≥10µM) retained with an `is_measured_inactive` flag; the weighting function itself is Phase 3A's job.
7. **Freeze protocol** (item 7): this document + `CHECKSUMS.sha256` + `manifest.json`, following `V1_FREEZE.md`'s exact structure.

## A second bug found only during execution (not anticipated in planning)

`confidence_score__gte=8` at the ChEMBL API's activity-filter level does **not** reliably restrict results to genuine single-protein targets — the same failure mode `build_target_fishing_index.py`'s own docstring already warned about for v1. Confirmed directly: `CHEMBL372` (`target_type=ORGANISM`, `pref_name="Homo sapiens"`) passed every Stage 1-4 fetch's filter despite its actual assay confidence being 1 ("Target assigned is non-molecular"). **18.99% of the pre-fix native-human pairs (931,457 of 4,905,110) referenced an invalid target.** Fixed by applying the same `single_protein_target_ids.json` cross-reference v1 already used — Stage 0 built this file but Stage 7's first pass hadn't applied it. Re-run after the fix: 3,973,653 clean pairs, 0 remaining invalid targets (verified).

## Known limitations (frozen, disclosed, not hidden)

- **Not a true atomic snapshot.** Fetched via ChEMBL's live REST API over several hours, not a single bulk database dump. ChEMBL's release doesn't change mid-cycle, so no material drift is expected, but this is weaker than a true point-in-time snapshot.
- **API instability during fetch.** ChEMBL returned intermittent HTTP 500s throughout this rebuild (confirmed broad, not specific to any one query — even bare `status.json` failed at points). All fetch stages are checkpointed and resumable; self-healing wrapper scripts (`run_stage_until_done.sh`) retried through these episodes. No data loss — every stage's final count matches its independently-probed `total_count`.
- **Residual null `confidence_score`** on 271,821 of 3,973,653 native-human pairs (6.8%) — the winning (max-pchembl) record for those pairs lacks an `assay_chembl_id` to look up; not a bug in the backfill (0 of the 270,225 resolved assay ids themselves have a null value).
- **Document/assay-campaign and temporal split types remain blocked** — this rebuild retains `n_documents`/`assay_types_seen` per pair but not per-pair dates, so the temporal holdout Phase 1's G1 gate needs is still not buildable from this freeze. A future rebuild pass would need to add activity-level date retention.
- **The density-adaptive rule (density≥12, `PHASE1_DENSITY_ADAPTIVE_RULE.md`) must be re-fit against this index before use** — this rebuild materially changes reference density (the `Potency` and censored-relation additions roughly triple native-human pair volume vs. v1's 1,312,849), so the old threshold no longer describes the same distribution.

## Files (in `target_prediction_v2/phase2/data/`)

```
manifest.json                      (build provenance, this freeze's authoritative numbers)
release_manifest.json              (ChEMBL release pin)
CHECKSUMS.sha256                   (sha256 of every file below)
stage7_native_human_pairs.jsonl    (3,973,653 pairs — the primary deliverable)
stage7_orthologue_pairs.jsonl      (34,080 pairs — separate stratum, never merged)
stage6_fingerprints.npz            (1,590,210 packed Morgan fingerprints, radius=2/2048 bits)
stage6_fp_meta.json                (standardized SMILES -> scaffold, fp index, contributing molecule_chembl_ids)
single_protein_target_ids.json     (5,869 verified human single-protein target ids)
orthologue_target_map.json         (385 human targets -> confirmed non-human orthologues)
```

Raw intermediate pulls (`stage1_activities.jsonl` through `stage5b_orthologue_activities.jsonl`, `assay_confidence_map.json`, `stage6_std_cache.json`) remain in the same directory for auditability but are not part of the checksummed freeze — they're reproducible from the scripts in `target_prediction_v2/phase2/` given the pinned release.

## What re-opens this freeze

Only a new validation/rebuild pass should change:
- Any of the filter/method decisions in items 1-7 above.
- The single-protein-target allowlist or orthologue gene-symbol map (both tied to the pinned ChEMBL release).
- The salt-collapse or max-pchembl aggregation logic.

Adding per-pair document dates (unblocking the temporal/document split types) or re-fitting the density-adaptive rule against this index are legitimate, planned next steps (Phase 1's remaining item and Phase 3A respectively) — they don't need to re-open this freeze, but they DO need their own validation pass before being trusted.

## Regenerating / verifying

```bash
# Verify no drift from this freeze:
cd target_prediction_v2/phase2/data && sha256sum -c CHECKSUMS.sha256

# Rebuild from scratch (requires re-fetching from ChEMBL):
cd target_prediction_v2/phase2
python3 stage0_pin_release.py
python3 stage1_fetch_base_evidence.py   # or ./run_stage_until_done.sh for resilience
python3 stage2_fetch_censored.py
python3 stage3_fetch_potency.py
python3 stage4_fetch_inactives.py
python3 stage5a_build_orthologue_map.py
python3 stage5b_fetch_orthologue_activities.py
/path/to/rdkit-env/python3 stage6_standardize.py   # requires RDKit
python3 stage_confidence_backfill.py
python3 stage7_aggregate.py
```
