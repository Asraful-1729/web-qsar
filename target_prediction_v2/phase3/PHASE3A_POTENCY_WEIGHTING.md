# Phase 3A — Potency Weighting: Tested, Adopted

**Status: DONE — genuine positive result, adopt.** Distinct lever from H4b's already-closed `s^α` similarity weighting (operates on pchembl, not Tanimoto) — legitimately open, now tested against n=1,600 pooled v2 captures, k=10 unweighted pooling as the base comparator, scaffold-clustered BCa.

## Result: positive at every tested configuration, no exceptions

11 (threshold, floor, steepness) configurations tested across two rounds. **Every single one showed a statistically significant positive effect** on all three metrics (any-annotated all, any-annotated dense-only, primary-target) — unlike popularity correction, this is not parameter-sensitive in direction, only in magnitude.

| Config | any-tier (all) | any-tier (dense≥8) | primary-tier |
|---|---:|---:|---:|
| Function default (thr=6, floor=0.2, steep=2.0) | +0.0385 | +0.0280 | +0.0454 |
| thr=6, floor=0.0, steep=2.0 | +0.0408 | +0.0283 | +0.0515 |
| thr=6, floor=0.2, steep=5.0 | +0.0412 | +0.0296 | +0.0500 |
| **thr=6, floor=0.0, steep=5.0 (best)** | **+0.0431** | **+0.0293** | **+0.0548** |
| thr=6, floor=0.0, steep=10.0 | +0.0423 | +0.0285 | +0.0533 |

All values are 95% CIs excluding zero (`ci_excludes_zero: True` in every case, `phase3/results/potency_weighting_test.json`).

## Recommended parameterization

**threshold=6.0 (pChEMBL), floor=0.0 (true floor, not the function's softer 0.2 default), steepness=5.0 (steeper ramp than the 2.0 default)** — best or tied-best on all three metrics. A true floor (weak-potency evidence contributes ~0, not a guaranteed minimum) combined with a steep transition around pChEMBL=6 (~1µM) outperforms the original softer default.

## Why this differs from H4b's null result (not a contradiction)

H4b tested `similarity^α` — reweighting neighbours by how structurally close they are. Closed as a genuine null: no significant effect at any tested width. This tests something orthogonal — reweighting each neighbour's **vote** by how potently that specific compound was measured against that specific target, independent of how similar the neighbour is to the query. The two levers don't overlap, and only one is a null.

## Next step

Fold `potency_weight(threshold=6.0, floor=0.0, steepness=5.0)` into the production density-adaptive pooling formula for the dense (density≥8) regime — combinable with the density rule directly (both operate on the same k=10 pooled vote, one gates *whether* to pool, the other weights *how much each vote counts*).
