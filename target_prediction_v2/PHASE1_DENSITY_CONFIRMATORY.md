# Phase 1 Entry Requirement: Density Finding, Confirmed at Scale (n=800/set)

**Answers Rev 7 §2.2 items 1-3 / `BUILD_PLAN.md` Phase 1 items 8-10.** This is the pre-registered confirmatory re-test `PHASE0B_DENSITY_STRATIFICATION.md` and the Rev 7 response explicitly called for before treating the density finding as established — run now on 4x the original sample (800 queries/set, up from 200), under the same scaffold-strict leakage control.
**Code:** `phase0b/phase1_density_confirmatory.py`. **Raw results:** `phase0b/results/phase1_density_confirmatory_{A,B}.json`.

**Headline: the core finding holds. Two specific claims from the exploratory pass do not, and are retracted here.**

---

## 1. The core finding replicates — at a more modest, now-reliable effect size

Median split (400/half), same test as the original exploratory pass but 4x the data:

| | Set A sparse (n=400) | Set A dense (n=400) | Set B sparse (n=400) | Set B dense (n=400) |
|---|---:|---:|---:|---:|
| mean density (≥0.5) | 1.86 | 55.95 | 3.42 | 64.73 |
| 10-NN vs best-sim, Top-1 | +0.01, CI[−0.04, 0.06] — **not sig.** | **+0.05, CI[0.02, 0.08] — sig.** | +0.008, CI[−0.04, 0.05] — **not sig.** | **+0.063, CI[0.03, 0.10] — sig.** |
| 10-NN vs best-sim, MRR | −0.004, CI[−0.03, 0.03] — not sig. | **+0.019, CI[0.00, 0.04] — sig.** | −0.001, CI[−0.03, 0.03] — not sig. | **+0.037, CI[0.02, 0.06] — sig.** |
| primary-target Top-1 | +0.037, CI[−0.07, 0.12] — not sig. | **+0.088, CI[0.02, 0.15] — sig.** | +0.022, CI[−0.09, 0.12] — not sig. | **+0.114, CI[0.03, 0.19] — sig.** |

**Confirmed**: pooling's benefit is real and statistically robust in dense neighbourhoods, and genuinely absent (not just "smaller") in sparse ones, in **both** sets, at 4x the original sample. This is the load-bearing result and it survived the confirmatory re-test.

## 2. What changed from the exploratory (n=200) numbers — stated plainly

**Effect sizes were overstated in the exploratory pass, as expected for an underpowered exploratory split** (textbook regression-to-the-mean on replication, not a methodology error — this is exactly why Rev 7 required this re-test before treating the finding as established):

| | Exploratory (n=100/half) | Confirmatory (n=400/half) |
|---|---:|---:|
| Set A dense half, Top-1 Δ | +0.11 | **+0.05** (less than half) |
| Set A sparse half, Top-1 Δ | −0.06 (leaning negative) | **+0.01** (now flat, not negative) |

The dense-half effect is real but roughly half the size first reported. The sparse-half effect, which leaned negative and fed directly into the "pooling can actively hurt" framing, is now indistinguishable from zero — genuinely flat, not harmful.

## 3. Retracted: the Q2 "significant valley" does not replicate

The exploratory quartile split found pooling *significantly hurting* in the low-moderate density band (Q2 MRR CI entirely below zero). At decile resolution on 800 queries, the corresponding band (Set A, D2: mean density 0.26) still points negative (Top-1 Δ = −0.0625) but its CI is [−0.191, 0.027] — **crosses zero, not significant**. Rev 7 itself flagged this exact risk when reviewing the exploratory result: *"Q2's dip could still be quartile-boundary noise rather than a true local minimum."* It was. **This specific claim — that pooling actively hurts in a mid-density band — is retracted.** The honest current state is: sparse neighbourhoods show no benefit (§1), not a penalty.

## 4. New finding: decile resolution needs more than n=80/bin at this effect size

Fitted at 10-decile resolution (n=80/decile), only the top 1–2 deciles reach individual significance in either set (Set A: D9 significant, D8 borderline; Set B: only D6, inconsistently). Every other decile's point estimate is broadly consistent with the monotonic dense-helps-more pattern but individually too noisy to distinguish from zero at n=80. **Practical implication for Phase 1's benchmark card and Phase 4's calibration strata**: use coarser bins (halves, terciles, or a fitted continuous curve with pooled degrees of freedom) rather than fine deciles unless the full-population sample (thousands per set) is used — a decile split on a few hundred queries will systematically under-detect a real, monotonic effect of this size.

## 5. Set A's overall effect: real, but much smaller than Set B's

Full-population aggregate (not density-split), n=800 each:

| | Set A | Set B |
|---|---:|---:|
| Top-1 Δ | **+0.030, CI[0.001, 0.057] — sig.** | **+0.035, CI[0.008, 0.062] — sig.** |
| MRR Δ | +0.008, CI[−0.010, 0.026] — not sig. | **+0.018, CI[0.0002, 0.037] — sig.** |
| primary-target Top-1 Δ | **+0.071, CI[0.019, 0.126] — sig.**| **+0.076, CI[0.013, 0.138] — sig.** |

**Revision to the Rev 6/Rev 7 framing**: at n=200, Set A's H4a effect was read as "not established" (no significant any-annotated result in either leakage setting). At 4x the sample, it *is* now detectable — small, but real (Top-1 and primary-target both clear significance). **H4a is not "Set B only."** It is present in both sets; it is simply smaller in Set A and concentrated almost entirely in the denser half (§1) — exactly consistent with, not contradicting, the density mechanism. The earlier "set-dependent" framing should be read as "density-dependent, and Set A happens to have less density-favourable mass," not "absent in Set A."

## 6. Flag-3, done properly this time: depth remains predictive even at high density, within every decile

The earlier linear joint model (Rev 7 response) was the wrong tool and settled nothing. Redone as a within-decile stratified check (does target reference depth still separate recovered/not-recovered queries *inside* each density decile):

| Decile (Set A) | mean density | depth gap (recovered − not-recovered) |
|---|---:|---:|
| D1 (sparsest) | 0.0 | 1,204 |
| D5 | 5.1 | 1,022 |
| D9 | 53.1 | **1,737** |
| D10 (densest) | 179.4 | 719 |

**Depth stays predictive across the entire density range, including the two highest deciles** — it does not wash out once density is high, at least not within Set A. This refines (not contradicts) the earlier finding that depth was flat within Set B but predictive within Set A: at finer resolution, depth carries real, independent information throughout Set A's density spectrum. **Conclusion for the joint model question (flag-3): density and depth are not interchangeable — both should be retained as separate stratification variables**, not collapsed to density alone.

## 7. Spread diagnostic: still genuinely open

Overall (n=800, narrow vs. wide top-10 similarity spread): mean outcome +0.0093 vs. +0.0064 — negligible difference, no clear signal at this resolution. Within the empirically-identified (not assumed) valley decile D2 specifically (n=80): narrow half −0.026 vs. wide half −0.055 — still directionally consistent with the original dilution hypothesis, but D2's own aggregate effect is no longer significant (§3), so this sub-analysis is even lower-powered than before. **Unresolved, carried forward exactly as flagged last round** — not strengthened, not weakened by this pass.

---

## 8. What this means for the build plan

1. **§1's core mechanism is confirmed and should be built on** — the density-adaptive retrieval architecture in `BUILD_PLAN.md` §3 stands.
2. **§2/§3's corrections must propagate into the benchmark card and G2 comparator language**: do not describe the dense-half effect as "+0.11 Top-1" (the real, replicated figure is roughly half that), and do not describe sparse/mid-density regimes as "actively harmful" — describe them as "no established benefit," which is a different and more defensible claim.
3. **§5 changes the retrieval-rule framing**: rather than "pool in dense regimes, use best-similarity in sparse ones, with Set B benefiting and Set A not," the correct framing is "pool where density is high, regardless of which set a query came from" — Set A does benefit, just less often, because less of its mass sits in the dense regime. The density-adaptive rule (§ Phase 1 item 11 in `BUILD_PLAN.md`) should condition on density and (per §6) depth, never on which evaluation set a query happens to belong to.
4. **§4's methodological finding is a real, standalone deliverable**: Phase 1's benchmark card and Phase 4's calibration strata should use coarse (half/tercile) or continuous-fitted density bins, not fine deciles, unless working with the full population.
5. **§6 settles part of decision #17** (`BUILD_PLAN.md` §9 was silent on this, `rev7` flag-3): the joint density+depth model should ship as two retained stratification variables, not one.
6. **§7 remains an open item**, correctly still flagged rather than forced to a conclusion either way — carry into Phase 3A's first-week work as originally planned.

**This confirmatory pass is now the authoritative version of the density finding.** `PHASE0B_DENSITY_STRATIFICATION.md` (the original 200-query exploratory pass) should be read as superseded on the specific numeric claims corrected in §2–3 above; its methodology and the core discovery it made remain valid and are what led here.
