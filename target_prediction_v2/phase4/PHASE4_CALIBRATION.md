# Phase 4 — Calibration: Tested, L-score Wins, Ship L-score

**Status: DONE, with an honest negative result for the more complex approach.** `BUILD_PLAN.md` Phase 4 item 5's own gate: "Must beat the L-score baseline on Brier and ECE, or ship the simple L-score instead — no complexity credit without earning it." Tested for real, not assumed.

## Scope actually covered (disclosed, not a full Mondrian cross-conformal spec)

- **Stratification: density only** — the one dimension this program has rigorously validated (just re-fit, `phase3/PHASE3A_DENSITY_REFIT.md`). Target-size and drug/non-drug strata (item 2's full spec) are **not** layered in here — real scope reduction, not hidden.
- **Split types: scaffold-strict captures only.** Document/temporal splits remain genuinely blocked (no per-pair dates in the rebuilt index — disclosed since Phase 1). Random/scaffold holdouts (`phase1/data/holdouts.json`) were not used for a formal three-way train/calibrate/test split in this pass — a simplification, stated plainly.
- **Method: isotonic regression per density band**, out-of-fold (5-fold scaffold-grouped CV, same folds as the 3E stacker test). Platt/beta scaling for small strata (item 4's spec) was **not implemented** — bands under 30 training points fall back to a global fit instead.

## Result

| | Brier score | ECE |
|---|---:|---:|
| **L-score baseline** (l/10, PHASE0B_ADDENDUM.md §7's convention) | **0.1532** | 0.0384 |
| Density-stratified isotonic-calibrated potency-weighted score | 0.1604 | **0.0360** |

n=1,600, base rate (top-1 correct, any-tier) = 0.590. Density band sizes: 0-8 (n=585), 8-25 (n=439), 25-100 (n=427), 100+ (n=149).

**Does not beat L-score on both metrics** — worse on Brier, marginally better on ECE. Per the plan's own explicit rule, this fails the gate.

## Disposition

**Ship L-score (l/10) as the shipped confidence/reliability measure, not the density-stratified calibration built here.** This is exactly the outcome the plan's gate exists to catch — added statistical machinery that doesn't earn its complexity. Plausible reasons this particular calibration attempt underperformed: only one stratification dimension used (target-size/drug-non-drug might matter more than density alone for *this* specific task — predicting whether the single top-ranked target is correct, a coarser question than the ranking-quality metrics density stratification was validated against); Platt/beta not implemented for small strata; no formal held-out calibration set separate from the scaffold-strict captures.

**Not a dead end — a scoped, honest first attempt that didn't clear the bar.** A follow-up with the full stratification (target-size × drug/non-drug), Platt/beta for small strata, and a proper three-way split could plausibly close the gap — not attempted here given the scope already covered. Left as a next step, not silently abandoned.
