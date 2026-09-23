# Target Fishing — Scaffold-Split Validation Benchmark

**Date:** 2026-09-20
**Script:** `backend/scripts/target_fishing_benchmark.py`
**Full report (machine-readable):** `target_fishing_benchmark_report.json`
**Raw per-instance results:** `target_fishing_benchmark_report_instances.csv` (7,688 rows)

## The question this answers

> When the reference pool has never seen the query compound, or anything sharing its Murcko scaffold, can target-fishing still recover a target independently known to interact with that compound?

This is the only benchmark design that tests scaffold-hopping generalization rather than database lookup. A random compound split would let near-identical analogs leak between "reference" and "test," inflating apparent performance.

---

## 1. Methodology

### Dataset
- **857,232** distinct compounds
- **1,312,849** compound–target pairs
- **4,658** human single-protein targets
- Source: broad ChEMBL pull (IC50/Ki/Kd/EC50, assay confidence ≥8, pChEMBL present, cross-verified against ChEMBL's own single-protein target list) — see `TARGET_FISHING.md`-adjacent history / `target_fishing_index/manifest.json` for full provenance.

### Scaffold-level split (leakage control)
1. Every distinct compound grouped by its precomputed Murcko scaffold. A compound with no ring system (~0.36% of the pool) is treated as its own singleton scaffold group.
2. **Eligible holdout scaffolds: groups with 1–50 members only.** This cutoff is not arbitrary — the real scaffold-size distribution is catastrophically skewed: median group size is 1, but the single largest scaffold is *benzene* (`c1ccccc1`) with **14,381 members** (1.6% of the entire pool), and several other generic fragments have 1,000+ members. Holding out a mega-scaffold would gut the reference pool for unrelated targets and produce an incoherent test set. Verified empirically before writing the benchmark, not assumed.
3. Scaffolds randomly sampled (seed=42) until reaching **5,001 held-out compounds** across **1,961 scaffold groups**. A scaffold group is always held out whole, never split.
4. **Every row for a held-out compound — across ALL targets, not just the one being evaluated — is removed from the reference index.** Ground truth still comes from the original, pre-removal data.

### Ground truth & evaluation instances
- For each held-out compound, every distinct target it was originally associated with is a **separate evaluation instance** (a 5-target compound = 5 independent recovery opportunities, not one pass/fail).
- Capped at **20 targets per compound** (randomly sampled if more) so one extremely promiscuous compound (the real max in this pool: **1,810 distinct targets for a single compound**) can't dominate the aggregate statistics.
- Result: **5,001 held-out compounds → 7,688 evaluation instances.**
- Removing the held-out compounds caused **7 targets to lose all their evidence entirely** (4,658 → 4,651 targets remaining in the reduced reference) — an expected, informative consequence for very sparse targets.

### Search
- Each held-out compound queried at `threshold=0.0` (so a true target's rank is always defined, never silently excluded) against the reduced reference, via `target_fishing.py`'s own `_search_against()` — the exact production aggregation/ranking code, not a reimplementation.
- **Two rankings scored for comparison:** by `evidence_score` (what real users see) and by `best_similarity` alone (naive-similarity baseline) — this directly tests whether the evidence-score formula improves ranking quality over raw Tanimoto.

### Primary / secondary split (near-exact control)
If a held-out compound's own top hit — even after removing its entire scaffold group — still comes back at `best_similarity ≥ 0.99`, that compound's instances are routed to a **secondary** bucket (near-exact retrieval) instead of **primary** (novel-scaffold generalization). A folded 2048-bit Morgan fingerprint can collide for distinct structures, and a stereoisomer/tautomer/salt-form variant can get a technically-different Murcko scaffold despite being chemically near-identical to something still in the reference. Lumping those into the primary number would inflate it with what's still effectively lookup.

- **Primary (novel scaffold): 7,295 instances (94.9%)**
- **Secondary (near-exact): 393 instances (5.1%)** — reassuringly small, meaning the primary bucket isn't inflated by leakage.

### Metrics
Top-1/5/10/20 recovery, Mean Reciprocal Rank (MRR), random-baseline expected recovery (K / n targets in reduced reference = 4,651), and enrichment (observed / random) — computed for both rankings. Stratified by (a) the target's original reference depth and (b) the achieved similarity band.

---

## 2. Primary headline results (novel scaffold)

| Metric | by `evidence_score` | by `best_similarity` alone |
|---|---:|---:|
| Top-1 | 55.27% | **56.41%** |
| Top-5 | 88.06% | **89.51%** |
| Top-10 | 91.12% | **92.15%** |
| Top-20 | 92.78% | **93.57%** |
| MRR | 0.6927 | **0.7048** |

Enrichment over random chance ranges from ~216× (Top-20) to ~2,624× (Top-1) — the system is overwhelmingly concentrating true targets near the top of a 4,651-target ranking, far beyond chance.

### ⚠️ Key finding: `evidence_score` underperforms raw similarity

**On every single metric, in both the primary and secondary buckets, ranking by raw Tanimoto similarity alone beats ranking by `evidence_score`.** The gap is modest (~1–1.5 percentage points) but consistent and reproducible across all four Top-K cuts and MRR — not noise.

This makes mechanistic sense: `evidence_score`'s scaffold-diversity and potency bonuses can promote a target with several mediocre analog matches above a target with one truly excellent match. That's a reasonable thing to want when asking "how much should a human trust this evidence," but it actively works against pure target-identification accuracy.

**Conclusion:** `evidence_score` is **not shown to improve target-recovery accuracy** over naive similarity ranking — if anything, it costs a little. It should not be presented as a better *predictor*, only (at most) as a differently-purposed trust/robustness signal.

## 3. Secondary results (near-exact retrieval, n=393)

| Metric | by `evidence_score` | by `best_similarity` alone |
|---|---:|---:|
| Top-1 | 59.29% | 59.80% |
| Top-10 | 90.59% | 90.59% |
| MRR | 0.7158 | 0.7184 |

Only marginally better than primary — confirms the primary ("novel scaffold") number is a fair, non-inflated headline rather than being artificially starved by pathologically hard cases.

## 4. Stratification by original reference depth (primary only)

**The most actionable finding in this benchmark.**

| Supporting compounds for target (pre-removal) | n instances | Top-1 | Top-10 | Top-20 |
|---|---:|---:|---:|---:|
| 1–2 | 12 | 0% | 0% | 0% |
| 3–10 | 19 | 10.5% | 21.1% | 31.6% |
| 11–50 | 110 | 37.3% | 71.8% | 79.1% |
| 51+ | 7,154 (98.1%) | 55.8% | 91.8% | 93.3% |

**Target-fishing essentially cannot generalize for sparse targets.** With only 1–2 known actives for a target, holding out even one compound often removes all its evidence, and recovery is a flat zero. 98.1% of primary evaluation instances fall in the well-studied (51+) bucket — meaning the strong ~91% headline Top-10 number is really a statement about **well-characterized targets**, which is exactly where target-fishing is least needed. For a genuinely novel or under-studied target — arguably where this tool would be most valuable for a natural product — **do not trust it.**

## 5. Stratification by achieved similarity band (primary only)

| Similarity band | n instances | Top-1 | Top-10 |
|---|---:|---:|---:|
| 0.0–0.5 | 798 | 8.3% | 32.5% |
| 0.5–0.6 | 624 | 37.5% | 88.6% |
| 0.6–0.7 | 1,415 | 54.7% | 98.2% |
| 0.7–0.8 | 2,295 | 65.8% | 99.7% |
| 0.8–0.9 | 1,647 | 65.5% | 99.7% |
| 0.9–0.99 | 516 | 71.5% | 100% |

Clean, monotonic relationship between achieved similarity and recovery — validates that exposing "Minimum Tanimoto similarity" as a user-facing control is meaningful, not cosmetic.

**Quotable, defensible number: at achieved similarity ≥ 0.6, Top-10 recovery exceeds 98%.**

---

## 6. Recommendations

1. **Don't market `evidence_score` as more predictive than similarity** — the data shows it's slightly worse at pure target recovery. Either revert to similarity-primary sorting as the default, or keep `evidence_score` explicitly framed as a corroboration/trust signal, not an accuracy improvement.
2. **Surface reference depth to users.** A hit backed by 1–2 known actives should be flagged very differently from one backed by 50+ — the benchmark shows these have wildly different reliability (0% vs. 92% Top-10 recovery).
3. **The Minimum Tanimoto control is validated** — encourage users toward ≥0.6 for higher-confidence results; below 0.5, recovery is poor (8% Top-1, 33% Top-10).
4. **Not yet done:** target-family stratification (needs an additional ChEMBL protein-classification pull) and 3D/pharmacophore similarity — both reasonable follow-ups, not blockers.

---

## Appendix: reproducing this benchmark

```bash
cd backend
python -m scripts.target_fishing_benchmark --n-heldout 5000 --seed 42 \
    --out ../target_fishing_benchmark_report.json
```

Runtime: ~110 minutes for 5,001 held-out compounds (~1.3s/query against the 1.3M-row reduced reference; see `target_fishing.py`'s `_search_against()` performance notes — an earlier, non-vectorized version of this aggregation took ~20s/query at `threshold=0.0`, which would have made this benchmark impractical).
