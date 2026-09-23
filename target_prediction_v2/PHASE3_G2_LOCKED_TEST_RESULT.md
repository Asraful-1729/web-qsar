# G2 — The Locked Scaffold Test, Opened

**Status: OPENED, per explicit user authorization, after the fresh-holdout sanity check showed real positive signal.** This is a one-time, irreversible event. 795 compounds, sealed since Phase 1 (`phase1/build_holdouts.py`).

## A critical finding, caught before scoring anything

**476 of the 795 locked compounds (60%) had already been used somewhere in Phase 3 tuning** — `capture_scaffold_strict_v2_{A,B}.jsonl` and the fresh-holdout check both sampled directly from `eval_sets.json`, not from `holdouts.json`'s tuning-only subset, so nothing excluded the locked test's scaffold groups from earlier work. Checked explicitly before running anything on the locked set, not discovered after the fact.

**This report treats the 433 (of which 370 successfully resolved against the v2 index) never-before-touched compounds as the real, trustworthy G2 answer.** The full 846-query capture (clean + contaminated) is reported alongside for context only — it is not a valid unbiased test and should not be quoted as one.

## The result — clean subset (n=370, the real answer)

Same assembled v2 pipeline as the fresh-holdout check: density-adaptive gate (pool if density≥8 else best-similarity), potency weighting (threshold=6.0/floor=0.0/steepness=5.0), additive orthologue term.

| Comparison | Effect | 95% CI | Significant? |
|---|---:|---|:---:|
| v2 vs. v1 `best_similarity`, any-tier | **+0.034** | [0.012, 0.059] | **Yes** |
| v2 vs. v1 `best_similarity`, primary-tier (n=103) | **+0.072** | [0.021, 0.128] | **Yes** |
| v2 vs. internal unweighted-10NN baseline, any-tier | **+0.031** | [0.004, 0.059] | **Yes** |
| v2 vs. internal unweighted-10NN baseline, primary-tier | +0.011 | [-0.039, 0.062] | No |

228/370 queries fell in the pooling regime (density≥8), 142 fell back to best-similarity.

**G2's core bar — beat frozen v1 with a paired CI excluding zero — is cleared, on genuinely held-out data, on both metrics.** v2 also clears the internal "always pool" baseline on the broad any-tier task here (stronger than the fresh-holdout check showed, though same direction). The primary-tier-vs-baseline comparison is not significant — most plausibly underpowered (n=103; `PHASE1_REMAINING_ITEMS_9_12_13.md`'s own power calculation flagged this scale as adequate only for margins ≥0.03-0.05, and this comparison's point estimate, +0.011, sits below that), not evidence of a real null.

## Full set, for context only (n=846, 60% contaminated) — do not treat as a G2 answer

| Comparison | Effect | 95% CI | Significant? |
|---|---:|---|:---:|
| v2 vs. v1 `best_similarity`, any-tier | +0.034 | [0.017, 0.050] | Yes |
| v2 vs. v1 `best_similarity`, primary-tier (n=245) | +0.054 | [0.016, 0.092] | Yes |
| v2 vs. internal baseline, any-tier | +0.008 | [-0.014, 0.028] | No |
| v2 vs. internal baseline, primary-tier | +0.002 | [-0.030, 0.035] | No |

Directionally consistent with the clean subset (reassuring — the contamination didn't flip any conclusion), but not cited as the real answer given the known overlap with tuning.

## What this means

- **v2 is real.** On the one test built specifically to give an unbiased answer, using only the genuinely never-touched portion of it, v2 beats v1's actual shipped production system by a real, statistically supported margin, on both the broad and the harder primary-target task.
- **The locked test is now spent.** Both the full 795 and the clean 370-compound subset have been looked at. Neither should be reused for further tuning — doing so would retroactively contaminate this exact result. Any future refinement needs a fresh, newly-sealed holdout.
- **G1 still isn't satisfied** (adjudication study still needs a human reviewer) — this result answers "does the method work," not "is the benchmark's ground truth itself fully validated." Both matter for a full release, but this closes the more consequential of the two open questions for anyone asking "should I trust this enough to build on it."

## Disposition

**Green light to proceed with Phase B — building the live serving engine.** The method has now cleared the one test that matters most. Move to `backend/target_prediction_v2.py` next.
