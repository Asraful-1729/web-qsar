# Phase 3E — Logistic Stacker: Tested, Adopted (small, real win)

**Status: DONE.** Pointwise logistic regression over three base scores (native k=10 pooling, adopted potency-weighted pooling, orthologue additive score) as per-(query,candidate-target) features, trained out-of-fold (5-fold scaffold-grouped CV, `analyze_h4.group_kfold`), compared against potency-weighted-alone (the current best single-lever) on the same n=1,600 captures.

36,137 (query, candidate) rows, 3,639 positive.

| Tier | Mean paired RR diff (stacked − potency-alone) | 95% CI | Significant? |
|---|---:|---|:---:|
| Any-annotated | +0.0064 | [0.0024, 0.0115] | Yes |
| Primary-target | +0.0163 | [0.0088, 0.0266] | Yes |

Small but genuine, on both tiers, no exceptions.

## Disposition, per `BUILD_PLAN.md`'s own gate

"Logistic regression on out-of-fold base scores first; GBM only if it wins and calibration survives." Logistic regression **wins** here — adopted. GBM escalation is explicitly gated on calibration surviving too, which is Phase 4's job (not yet run) — **not escalating to GBM until Phase 4's calibration pass confirms the current stack holds up**, per the plan's own sequencing, not a shortcut.
