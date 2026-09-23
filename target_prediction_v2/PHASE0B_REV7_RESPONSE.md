# Response to Rev 7

**Answers Rev 7's three flags and states what's out of scope for now.** Code: `phase0b/test_rev7_flags.py` (flags 2/3), inline queries (flag-1). Raw results: `phase0b/results/rev7_flags_report.json`.

**G0b closure acknowledged**: agreed with Rev 7 §0 that the grid-widening null and the density finding together close H4b and give H4a a real mechanism. No dispute with §0 or §1.

---

## Flag-1 (blocking) — answered

Rev 7 correctly caught that the addendum reported the *pair-level* primary rate (2–6%) but not the *compound-level* coverage Rev 5 Decision 3 actually needs for scoping primary-target metrics. Computed directly, on the **full matched population** (not the 200-query sample):

| Set | n (full matched) | Compounds with ≥1 primary/intended target | Coverage |
|---|---:|---:|---:|
| A (approved) | 1,965 | 679 | **34.6%** |
| B (clinical) | 2,750 | 790 | **28.7%** |
| D (sparse-target) | 334 | 141 | **42.2%** |

(679 + 790 = 1,469 of 1,965 + 2,750 = 4,715 — matches the original Phase 0 combined 31.2% figure exactly, internally consistent.)

**All three sets sit well below the ~60% adequacy bar Rev 5 implies for primary-target metrics to be meaningfully powered.** This is a real constraint, not a rounding issue — primary/intended-tier Top-k/MRR/median-rank should be reported as a supplementary, lower-power check in Phase 1, not relied on as a primary metric.

The broader **mechanistically-supported** tier (Rev 5's 2-of-3 rule) does much better, computed on the 200/200/334 *sample* (mechanism-support pull was scoped to the sample, not the full population — noted, not hidden):

| Set | Compounds with ≥1 mechanistically-supported target |
|---|---:|
| A | 50.5% |
| B | **63.5%** — clears the bar |
| D | **74.3%** — clears the bar |

**Recommendation for Phase 1**: use mechanistically-supported (not strict primary) as the powered ground-truth tier for the harder ranking metrics on Sets B and D; treat Set A's primary/mechanistic metrics as underpowered on current coverage and report them with that caveat rather than silently proceeding as if 200+ queries meant adequate power.

---

## Flags 2 and 3 — preliminary only, and one honest methodological catch

Both came back with a genuine limitation in *how* I checked them, which matters more than the numbers themselves.

**What went wrong with the first attempt:** I tested both with linear correlation / OLS on the raw per-query outcome (reciprocal-rank difference, unweighted-10NN vs. best-similarity) against density and spread. Both came back essentially null — R² < 1%, no coefficient significant. **This does not mean density or spread don't matter — it means linear correlation is the wrong tool for a relationship already shown to be non-monotonic** (the Q2 valley: negative, then positive, then flattening). A linear fit across a valley-then-rise shape averages the negative and positive halves toward zero by construction. This would have been a misleading "null result" if reported without the caveat, so it's reported *with* the caveat instead.

**Redone properly, matching the method that actually found the density effect** (group comparison, not linear correlation) — within the Q2 band itself (density 1–7 neighbours ≥0.5, n=74), split by similarity spread (gap between the 1st- and 10th-ranked neighbour's similarity):

| | n | mean spread gap | mean Δ(reciprocal rank) | mean Δ(Top-1 hit) |
|---|---:|---:|---:|---:|
| Narrow-spread half of Q2 | 37 | 0.140 | −0.023 | 0.000 |
| Wide-spread half of Q2 | 37 | 0.298 | **−0.066** | **−0.081** |

**Directionally supports flag-2's dilution hypothesis**: within the already-identified danger zone, wider spread among the top-10 neighbours makes pooling's penalty worse, not better — consistent with "a couple of close neighbours diluted by several mediocre, spread-out ones" being the actual mechanism, not density alone. At n=37/half this is preliminary and underpowered, exactly as Rev 7 anticipated it would be at this sample size — reported as a directional signal worth confirming, not a settled finding.

**Flag-3 (joint density + depth model)**: the same linear-OLS limitation applies, so the near-zero, non-significant coefficients for both density and depth in this preliminary fit **should not be read as "depth adds nothing"** — the model itself likely can't detect either variable's real (non-monotonic) effect. This flag needs the same fix as flag-2: Phase 1 should test the joint model with a specification that can represent a non-monotonic density effect (e.g., density as a categorical/spline term, or the model fit within density bands) rather than a flat linear term. Redoing this properly is Phase 1 scope, not a quick rerun — flagged as still genuinely open, not preliminarily answered either way.

---

## Decision #18 (confirmatory re-test on the full benchmark) — correctly out of scope now

Agreed with Rev 7: this needs the full ~4,700-drug matched population (not the 200/200/334 exploratory sample), which is a substantially larger capture-and-analysis job than anything run in Phase 0b. Not attempted here. Explicit sign-off, as requested: **the density finding (and its flag-2/3 refinements) is a pre-registered hypothesis for Phase 1, not an established fact carried over from this exploratory round.**

## What's now settled vs. still open, going into Phase 1

| | Status |
|---|---|
| Compound-level primary coverage (flag-1) | **Answered.** 28.7–42.2% across sets, below the adequacy bar; mechanistically-supported tier clears it for B/D not A. |
| Density is the real variable behind H4a's set-dependence | Established in Phase 0b, **pending Phase 1 confirmatory re-test** (decision #18) |
| Similarity spread contributes beyond density (flag-2) | **Directionally supported**, preliminary, underpowered (n=37/half) — needs Phase 1's proper test |
| Target-depth adds nothing once density is modeled (flag-3) | **Not actually tested yet** — the linear check that was run couldn't have detected it either way; still open |
| Full-benchmark confirmatory rerun (decision #18) | **Out of scope for Phase 0b**, correctly assigned to Phase 1 |
