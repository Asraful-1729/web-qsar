# Density-Stratification Test: Confirmed on Set A Alone

**Tests the hypothesis from `PHASE0B_DIAGNOSIS_A_VS_B.md`**: that neighbourhood density, not "approved vs. clinical" as such, is the real variable behind H4a (pooling beats best-similarity) being significant for Set B and not for Set A.
**Code:** `phase0b/test_density_stratification.py`. **Raw results:** `phase0b/results/density_stratification_report.json`. Runs entirely on **Set A's own 200 queries**, under the scaffold-strict leakage control (the same control the H4a/H4b conclusions rest on) — no Set B data involved at all.

**Confirmed.** Splitting Set A purely by neighbourhood density (n neighbours at Tanimoto ≥0.5) reproduces the Set A/B pattern almost exactly, entirely within Set A.

## Median split: clean, significant confirmation

| | Sparse half (n=100) | Dense half (n=100) |
|---|---|---|
| mean neighbours ≥0.5 | 1.79 | 50.9 |
| mean max-similarity-to-index | 0.551 | 0.716 |
| **10-NN vs best-sim, Top-1** | −0.06, CI [−0.15, 0.01] — **not significant, leans negative** | **+0.11, CI [0.04, 0.19] — significant** |
| **10-NN vs best-sim, MRR** | −0.047, CI [−0.11, 0.01] — **not significant, leans negative** | **+0.052, CI [0.007, 0.11] — significant** |
| 10-NN vs best-sim, primary Top-1 | −0.048, CI [−0.35, 0.10] — not significant (n=21) | **+0.209, CI [0.07, 0.36] — significant (n=43)** |

The dense half of Set A behaves like Set B did in the Rev 6 report — a real, bootstrap-confirmed pooling benefit on every metric. The sparse half behaves like Set A's aggregate did — no significant benefit, point estimates actually negative. **This split was made entirely within one set, using only a neighbourhood-density count**, with no reference to drug-development phase at all.

## Quartile dose-response: mostly monotonic, with one genuine wrinkle

| Quartile | mean density (≥0.5) | Top-1 Δ | MRR Δ | Significant? |
|---|---:|---:|---:|:---:|
| Q1 (sparsest) | 0.38 | −0.06 | −0.023 | No |
| Q2 | 3.2 | −0.06 | **−0.070** | **Yes — significantly negative** |
| Q3 | 12.7 | **+0.18** | **+0.089** | **Yes — significantly positive** |
| Q4 (densest) | 89.1 | +0.04 | +0.014 | No (any-annotated) — but primary Top-1 **+0.269, CI [0.12, 0.48], significant, n=26** |

Not a clean monotonic ramp. Two things worth stating plainly rather than smoothing over:

- **Q2 is a real, statistically significant negative result** (MRR CI entirely below zero) — at low-to-moderate density, pooling doesn't just fail to help, it measurably hurts. A plausible mechanism: with only a handful of genuinely close neighbours plus several mediocre ones filling out the k=10 pool, the vote gets diluted by noise that best-similarity's single-best-pick correctly ignores. Only once density is high enough (Q3+) does the k-pool consist mostly of *good* neighbours, flipping the noise/signal ratio in pooling's favour.
- **Q4's any-annotated effect is not significant, but its primary-target effect is the strongest and tightest of any stratum.** At very high density, best-similarity itself gets good enough that the any-annotated gap narrows — but pooling still measurably sharpens the specific, harder primary-target ranking.

## What this settles, and what it doesn't

**Settled:** the H4a effect is a density effect, not a drug-phase effect. Any Phase 1/3A baseline recommendation should be conditioned on neighbourhood density (or an equivalent evidence-richness signal), not on which evaluation set a query happens to belong to. This directly supports D5 (explicit abstention) and the reliability/L-score work (Phase 4) being aimed at the *right* variable — "how much evidence exists near this query" — rather than something set-specific and non-generalizable.

**Not settled:** why the relationship dips at Q2 before rising at Q3 rather than climbing smoothly. Worth a finer-grained sweep (deciles rather than quartiles) in Phase 1 once the full benchmark is built and sample size supports it — at n=50/quartile here, Q2's dip could still be quartile-boundary noise rather than a true local minimum, though its CI is clean and entirely negative, which argues against dismissing it outright. Flagged as open, not resolved.

**Practical recommendation:** the benchmark card and any G2 comparator should stratify by neighbourhood density (a continuously available, query-time-computable signal, unlike "is this an approved drug") as a primary reporting axis — this is likely the same variable underlying the original H1 "reference-evidence bucket" finding and the H5 threshold-sweep coverage tradeoff, and unifying them under one density axis is more useful than three separate observations.
