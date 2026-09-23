# Phase 3A — Popularity Correction: Tested, Rejected

**Status: DONE — genuine negative result, do not ship.** `score.py::popularity_correct()` already existed (observed/expected vote-mass correction, `score / (popularity + eps)`, rev5 §3's "same logic as SEA's background model") but was flagged in `BUILD_PLAN.md` as "untested lever." Tested here against the n=1,600 pooled v2 captures (same data as the density re-fit), at k=10 unweighted pooling (matching the density-adaptive rule's own comparator).

**Real target popularity computed from the rebuilt index** (`phase3/results/target_popularity_v2.json`): 4,882 targets, distinct-compound-count range 1–195,809 (mean 814, median 8) — a severe, real skew, confirming the correction targets a real problem in principle.

## Result: significantly negative at every tested eps (0.1, 1.0, 5.0, 20.0, 100.0)

| Tier | Scope | Mean paired RR diff (corrected − raw) | 95% CI (scaffold-clustered) |
|---|---|---:|---|
| Any-annotated | all n=1600 | −0.040 to −0.058 | excludes zero, every eps |
| Any-annotated | dense only (density≥8) | −0.059 to −0.076 | excludes zero, every eps |
| Primary-target | all n=518 | **−0.105 to −0.125** | excludes zero, every eps |
| Primary-target | dense only n=381 | **−0.124 to −0.142** | excludes zero, every eps |

Direction is consistent and monotone in eps (smaller eps → stronger penalty → larger harm) but **never crosses into positive territory at any tested value** — this isn't a tuning problem, the correction is harmful in this form regardless of strength.

## Why (plausible mechanism, not just an empirical shrug)

The any-annotated ground-truth tier is deliberately broad — and popular targets are popular *because* they're genuinely promiscuous/commonly active, not only because they're over-represented in the index. Dividing by popularity penalizes exactly the targets most likely to be **true** annotations under this tier's own definition. The effect is worse on the primary-target tier, consistent with popular targets often *being* the correct primary mechanism for well-studied compound classes (e.g., kinases), not just index noise.

## Disposition

- **Do not ship `popularity_correct()` as part of the retrieval pipeline.** `BUILD_PLAN.md` §3's architecture diagram position (between density-adaptive retrieval and potency/confidence weighting) is not adopted.
- **D2's "de-flooded by 3A's popularity correction" framing (`BUILD_PLAN.md` §Phase 5) needs revisiting** — polypharmacology grouping cannot rely on this mechanism to de-flood popular targets; a different de-flooding approach (or none) is needed before D2 ships. Flagged, not solved here — out of this test's scope.
- The function stays in `score.py` (harmless, unused) with this finding referenced in its docstring for anyone tempted to re-enable it without re-checking.
