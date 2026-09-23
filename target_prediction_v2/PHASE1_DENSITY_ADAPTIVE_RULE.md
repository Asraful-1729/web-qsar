# Phase 1 Entry Requirement: The Fitted Density-Adaptive Retrieval Rule

**Answers `BUILD_PLAN.md` Phase 1 item 11 — "a fitted, continuous density-adaptive retrieval rule... not a hand-set quartile threshold."** Built directly on `PHASE1_DENSITY_CONFIRMATORY.md`'s confirmed finding that density (not which eval set a query came from) is the driver, and that both sets show the effect — so this pools Set A + Set B (1,600 queries total) rather than fitting two separate per-set rules.
**Code:** `phase0b/fit_density_adaptive_rule.py`. **Raw results:** `phase0b/results/density_adaptive_rule.json`.

---

## 1. The rule

**Use unweighted k-NN pooling when a query has ≥12 neighbouring reference compounds at Tanimoto ≥0.5; use best-similarity (singleton retrieval) below that.**

This is not a hand-picked number — it is the smallest density in a fine grid (0, 1, 2, 3, 5, 8, 12, 18, 25, 35, 50, 70, 100, 150, 200, 300) at which a scaffold-clustered 95% bootstrap band on an isotonic (monotonic, nonparametric) fit of pooling's benefit first excludes zero **and stays excluding zero at every higher density tested** — i.e., the threshold doesn't just clear significance once and wobble back, it holds all the way to the top of the observed range.

## 2. How it was fit

Two complementary fits, both on the same pooled 1,600-query dataset (per-query outcome: reciprocal-rank of unweighted-10-NN minus reciprocal-rank of best-similarity, any-annotated tier — the tier with the full sample):

1. **Isotonic regression** (`sklearn.isotonic.IsotonicRegression`, monotonic non-decreasing) — the rigorous fit. It assumes only that expected benefit does not *decrease* as density increases, a weaker and more defensible assumption than any specific functional form, and one now consistent with the confirmed findings (the one claim that would have violated it — the Q2 "valley" — was retracted in the confirmatory pass). A **scaffold-clustered bootstrap** (2,000 resamples, refitting isotonic regression on each resampled set of scaffold groups, reusing the same cluster-resampling logic as every other confidence interval in this program) gives a genuine confidence band on the *curve itself*, not just on isolated point estimates.
2. **A smooth logistic mixing-weight formula**, least-squares fit to the same raw data, in the exact form Rev 7 itself suggested ("the pooling/best-similarity mixing weight as a smooth function of local density"): `w(density) = floor + (ceiling − floor) / (1 + exp(−(density − d0) / s))`.

## 3. The fitted curve

| Density | Fitted benefit | 95% CI | Significant? |
|---:|---:|---:|:---:|
| 0–3 | −0.010 | [−0.040, 0.019] | No |
| 5 | +0.007 | [−0.026, 0.035] | No |
| 8 | +0.018 | [−0.010, 0.038] | No |
| **12** | **+0.029** | **[0.007, 0.042]** | **Yes — first crossing** |
| 18–35 | +0.029 | [0.01–0.02, 0.04–0.05] | Yes |
| 50–300 | +0.031 | [0.02, 0.05–0.12] | Yes (plateaus) |

The curve is flat and slightly negative below density ~3, rises through densities 5–12, then **plateaus around +0.03** from density ~50 onward — pooling's benefit doesn't keep growing indefinitely with more neighbours, it saturates. This matches the earlier decile-level observation that the very top decile (D10) wasn't always the strongest individual result — there are diminishing returns past a moderate density, not an ever-increasing one.

## 4. An honest, informative surprise: the "smooth" fit found a near-hard step, not a gradual ramp

The logistic fit's own optimized scale parameter came out at **s = 0.04** — for a density axis running from 0 to 300+, this is effectively an instantaneous transition, not a gradual sigmoid. The least-squares optimizer, given full freedom to find a gentle ramp if the data supported one, instead converged on something indistinguishable from a hard step at density ≈5:

| Density | Logistic fit |
|---:|---:|
| 0–3 | −0.010 |
| 5 | +0.007 |
| **8+** | **+0.028 (flat)** |

**This is a real finding, not a fitting artifact to explain away**: the relationship between neighbourhood density and pooling's benefit is much closer to an on/off switch than a continuum requiring a nuanced mixing weight. A simple threshold rule (§1) is not a crude simplification of some richer underlying continuous relationship — it is, empirically, close to the actual shape of the relationship.

## 5. What this changes in the plan

- **Phase 3A's retrieval core** (`BUILD_PLAN.md` §3) should implement this as a genuine density check at query time — count neighbours at Tanimoto ≥0.5 post-leakage-removal (already computed for every query in this program's captures; cheap to compute in production) — and branch: **≥12 → pool (unweighted or lightly-weighted k-NN, per the still-closed H4b null); <12 → best-similarity.**
- **No mixing weight is needed.** Rev 7 anticipated possibly needing a smooth interpolated weight between the two strategies; §4's finding says a hard branch is empirically adequate — simpler to implement and to reason about than a continuous blend would have been, and the data doesn't support the extra complexity.
- **This closes `BUILD_PLAN.md` Phase 1 item 11.** Combined with `PHASE1_DENSITY_CONFIRMATORY.md` (items 8 and 10), three of Rev 7's four Phase 1 entry requirements are now done. Item 9 (similarity-spread diagnostic) and item 12 (explicit unification test against H1's bucket / H5's abstention sweep) remain open.
- **The G2 comparator** (`BUILD_PLAN.md` §6) can now cite an actual number instead of a placeholder: the density-adaptive baseline is expected to show roughly +0.03 MRR-scale benefit over best-similarity in the ≥12-density regime, ~0 (not negative) below it — this is the real number the power calculation (Phase 1 item 13, still open) should be built from.

## 6. Caveats, stated directly

- The plateau value (+0.03) and the threshold (12) are both fit on the **any-annotated** tier specifically — the primary-target tier (smaller n, per `PHASE1_DENSITY_CONFIRMATORY.md`) was not separately re-fit here; its effect sizes ran larger throughout this program (e.g. dense-half primary-target deltas of +0.09 to +0.11), so a primary-target-specific threshold might legitimately differ (likely lower, since that tier showed benefit even in some lower-density strata). Worth a follow-up fit before finalizing the production rule if primary-target performance is a separate optimization target.
- This rule was fit on Sets A and B only (Set D, the sparse-target set, was not included in the pooled fit — it wasn't part of the n=800 scaffold-strict recapture this pass used). Given Set D is specifically the sparse-target population, it is disproportionately likely to sit below the density=12 threshold; worth confirming this threshold still makes sense once Set D is captured at the same scale, rather than assuming it transfers unchanged.
- The threshold is specific to the exact density definition used throughout this program (count of neighbours at Tanimoto ≥0.5, post scaffold-strict leakage removal) — a different similarity cutoff or leakage protocol would need its own fit, not a reuse of "12."
