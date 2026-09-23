# Target Fishing — v1 Freeze

**Frozen:** 2026-09-20
**Status:** Locked. Do not change without a new validation cycle (see "What re-opens this freeze" below).

This directory is a complete, standalone snapshot of everything that makes up Target Fishing v1 — code, the built index it runs against, and the validation report that justified its current ranking method. `CHECKSUMS.sha256` lets anyone verify later whether the live files in the repo have drifted from what was actually validated.

## Architecture (frozen)

```
Query molecule
   -> ChEMBL similarity (Morgan/ECFP4, packed-bit Tanimoto)
   -> Minimum Tanimoto threshold filter
   -> per-target aggregation (n_similar_actives, scaffolds, pChEMBL, evidence_score)
   -> ranked by BEST SIMILARITY (not evidence_score — see "Why best_similarity" below)
   -> reference-depth badge + full evidence metrics shown alongside
```

## Dataset (frozen)

| | |
|---|---:|
| Distinct compounds | 857,232 |
| Compound–target pairs | 1,312,849 |
| Human single-protein targets | 4,658 |
| Source | ChEMBL REST API — human targets, IC50/Ki/Kd/EC50, assay_confidence_score≥8, pChEMBL present, cross-verified against ChEMBL's own single-protein target list |

Full provenance: `target_fishing_index/manifest.json` in this snapshot.

## Why `best_similarity`, not `evidence_score` (frozen ranking decision)

A scaffold-split validation benchmark (`TARGET_FISHING_BENCHMARK.md`, `target_fishing_benchmark_report.json`, 7,295 novel-scaffold evaluation instances) found `best_similarity` beats `evidence_score` on **every** metric:

| Metric | evidence_score | best_similarity |
|---|---:|---:|
| Top-1 | 55.27% | **56.41%** |
| Top-10 | 91.12% | **92.15%** |
| MRR | 0.6927 | **0.7048** |

`evidence_score` is still computed and shown (as "Evidence score", explicitly labeled non-predictive), because it's a real, differently-purposed corroboration signal — just not the ranking key.

## Known limitation (frozen, documented, not silently hidden)

Recovery is heavily dependent on how much ChEMBL evidence exists for a target:

| Reference evidence | Top-10 recovery |
|---|---:|
| 1–2 similar actives ("Limited") | 0% |
| 3–10 ("Low") | 21% |
| 11–50 ("Moderate") | 72% |
| 51+ ("High") | 92% |

Surfaced to users directly as the "Reference evidence" badge on every result (`TargetFishingTab.tsx`), with a tooltip citing these exact numbers — not a vague confidence label.

## Files in this snapshot

```
target_fishing_v1_freeze/
├── V1_FREEZE.md                              (this file)
├── CHECKSUMS.sha256                          (sha256 of every file below, at freeze time)
├── TARGET_FISHING_BENCHMARK.md               (full validation methodology + results)
├── target_fishing_benchmark_report.json      (machine-readable benchmark results)
├── target_fishing_benchmark_report_instances.csv  (all 7,688 raw per-instance results)
├── target_fishing_index/
│   ├── fingerprints.npz                      (857,232-compound packed Morgan fingerprints)
│   ├── compounds.csv.gz                      (target_chembl, target_pref_name, target_id, smiles, pchembl_value, murcko_scaffold)
│   └── manifest.json                         (build provenance)
└── code/
    ├── target_fishing.py                     (-> backend/target_fishing.py)
    ├── build_target_fishing_index.py         (-> backend/scripts/build_target_fishing_index.py)
    ├── target_fishing_benchmark.py           (-> backend/scripts/target_fishing_benchmark.py)
    ├── TargetFishingTab.tsx                  (-> frontend/src/tabs/TargetFishingTab.tsx)
    └── types.ts                              (-> frontend/src/lib/types.ts)
```

## What re-opens this freeze

Only a new validation result should change any of the following:
- The ranking method (currently `best_similarity`).
- The `evidence_score` formula (`0.6·similarity + 0.2·scaffold_bonus + 0.2·potency_bonus`).
- The reference-evidence bucket boundaries (1–2 / 3–10 / 11–50 / 51+) or their labeled recovery rates.

Adding target-family stratification or 3D/pharmacophore similarity are legitimate future extensions (deliberately deferred, not rejected) — they don't need to re-open this freeze, but they DO need their own validation pass before being trusted the way this v1 report is.

## Regenerating / verifying

```bash
# Verify no drift from this freeze:
sha256sum -c CHECKSUMS.sha256   # run from inside this directory, paths relative to it

# Rebuild the index from scratch (requires re-fetching ChEMBL data — see
# scripts/fetch_chembl_activities.py's documented filter):
cd backend && python -m scripts.build_target_fishing_index \
    --input activities.jsonl --single-protein-ids single_protein_target_ids.json

# Re-run the validation benchmark:
cd backend && python -m scripts.target_fishing_benchmark --n-heldout 5000 --seed 42
```
