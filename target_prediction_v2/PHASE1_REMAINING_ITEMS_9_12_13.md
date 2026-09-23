# Phase 1 Entry Requirements: Items 9, 12, 13 — Closed

**Code:** `phase0b/phase1_remaining_items.py`. **Raw results:** `phase0b/results/phase1_remaining_items.json`. All computed from the existing n=800/set scaffold-strict captures (pooled, n=1,600) — no new data pull needed.

---

## Item 9 — similarity-spread diagnostic, properly powered: **the dilution hypothesis is not supported**

The earlier preliminary check (n=37/half) found wide top-10 similarity spread associated with a worse outcome specifically in the low-density "valley" band, and flagged it as directional-but-underpowered. Re-run banded by the actual fitted threshold (density <12 vs. ≥12, from `PHASE1_DENSITY_ADAPTIVE_RULE.md`) with real power (n=435–436 per half in the low band, not 37):

| Band | Narrow-spread outcome | Wide-spread outcome |
|---|---:|---:|
| Low density (<12), n=871 | +0.008, CI[−0.019, 0.035] — not sig. | −0.006, CI[−0.036, 0.024] — not sig. |
| High density (≥12), n=729 | +0.023, CI[0.002, 0.050] — sig. | +0.031, CI[0.012, 0.056] — sig. |

**In the low-density band, narrow and wide spread are statistically indistinguishable** — no detectable dilution effect. In the high-density band, both narrow and wide spread show the expected positive pooling benefit, with wide spread if anything numerically *higher* (opposite the dilution direction, though the two CIs overlap so this isn't a real reversal, just not the hypothesized pattern either).

**Conclusion: closed as a genuine negative result.** With adequate power, similarity spread does not explain or modulate pooling's benefit in the direction the original small-sample check suggested. Density remains the load-bearing variable; spread adds nothing detectable on top of it. Do not spend further Phase 3A time on a spread-aware weighting scheme.

## Item 12 — unification: density vs. max-similarity-to-index vs. target depth

| Comparison | Result |
|---|---|
| density vs. max-similarity-to-index | Pearson r = 0.491, **Spearman ρ = 0.827** |
| density vs. target reference depth | Already shown independently predictive within every density decile (`PHASE1_DENSITY_CONFIRMATORY.md` §6) |

**Density strongly rank-correlates with max-similarity-to-index (ρ=0.83)** — the two are measuring closely related phenomena (both "how close is the nearest real evidence"), so rev5's original instinct to use max-similarity as a stratification variable was pointing at the right general effect. Density is the better-justified variable of the two going forward: it counts *how much* evidence clears a real threshold, not just the single closest point, and it's what the actual retrieval rule (`PHASE1_DENSITY_ADAPTIVE_RULE.md`) is built on — but the strong correlation means most of what max-similarity captured is not lost by switching.

**Density does not subsume target depth.** Depth remains predictive even after conditioning on density (shown directly, not assumed).

**Conclusion: partial unification, not full.** Density + depth should both ship as retained stratification variables (as already recommended); max-similarity-to-index can be retired in favor of density without losing much, given the strong rank correlation between them.

## Item 13 — power calculation for G2's margin

Using rev5 §8.7's formula with real numbers from this program's own data (not assumed values):

- σ_diff = **0.2481** (paired reciprocal-rank difference, unweighted-10NN vs. best-similarity, any-annotated tier, pooled n=1,600).
- 1,263 distinct scaffold clusters, mean cluster size 1.27, **estimated ICC (ρ) ≈ 0.279** — real, non-trivial clustering, not negligible.
- **Design effect = 1.074** (a modest, ~7% sample-size inflation from scaffold clustering — smaller than it might have been, because mean cluster size is close to 1).

| Margin (Δ) | Required n per arm | Feasible given ~4,715-drug matched population? |
|---:|---:|:---:|
| 0.01 | 5,190 | **No** |
| 0.02 | 1,297 | Yes |
| 0.03 | 577 | Yes |
| 0.05 | 208 | Yes |
| 0.08 | 81 | Yes |

**Recommendation for decision #6 (`BUILD_PLAN.md` §9): set G2's margin at Δ ≥ 0.02.** A margin finer than 0.01 cannot be reliably detected even with the entire matched population this program has assembled — per rev5's own rule, that means widening the margin (done here) rather than running an underpowered gate and reading the point estimate.

**Explicit caveat, not glossed over**: σ_diff here is a **proxy** — it's the variance of the H4a (pooling-vs-best-similarity) paired difference, the best real data available, not the variance of the actual future G2 comparison (a v2 candidate vs. the density-adaptive baseline). Re-run this exact calculation once real v2-candidate paired differences exist, before finalizing the margin for the real gate — this is a placeholder built from real data, not a placeholder built from a guess.

---

**Status: all three items closed.** Combined with `PHASE1_DENSITY_CONFIRMATORY.md` and `PHASE1_DENSITY_ADAPTIVE_RULE.md`, all six of Rev 7's Phase 1 entry requirements (`BUILD_PLAN.md` items 8–13) are now done.
