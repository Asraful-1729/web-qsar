# Phase 3A — Orthologue Ablation: Tested, Mixed Real Result

**Status: DONE.** Tests item 4's orthologue tier (34,080 pairs, 385 human targets with a confirmed gene-symbol-matched orthologue) as an additive term on top of native-human k=10 pooling, on the same n=1,600 captures. Simple union rule: `combined_score(t) = native_score(t) + orthologue_score(t)`, both on the unweighted-vote scale, orthologue neighbours found via direct Tanimoto search (Tanimoto≥0.4, top-300 cap) against a dedicated orthologue-tier index.

**Coverage**: 788/1,600 queries (49%) had ≥1 qualifying orthologue neighbour — real, substantial applicability, not a rare edge case.

## Result: no broad effect, but a real, specific win on the harder task

| Tier | Scope | Mean paired RR diff | 95% CI | Significant? |
|---|---|---:|---|:---:|
| Any-annotated | all n=1600 | +0.0057 | [−0.0037, 0.016] | No |
| Any-annotated | dense only (density≥8) n=1015 | −0.0074 | [−0.0194, 0.0064] | No (slightly negative point estimate) |
| **Primary-target** | **all n=518** | **+0.0194** | **[0.0012, 0.041]** | **Yes** |

## Interpretation

Adding orthologue evidence doesn't move the easy, broad "any known target" task — sensible, since dense native-human evidence (which most any-tier hits already have) leaves little room for a small additive term to matter, and may even introduce mild noise in the dense regime (the slightly negative, non-significant dense-only point estimate). But it **does** produce a real, if modest, improvement on the harder, more specific primary-target recovery task — exactly where rev5's original orthologue-tier rationale expected it to matter: filling gaps where native human evidence for the *true* mechanism is thin. The effect is real but fragile (CI lower bound 0.0012, barely above zero) — a genuine positive, not a strong one.

## Disposition

**Adopt the orthologue term for primary-target scoring specifically, not as a blanket addition to all ranking.** Simple additive union is sufficient for this first pass — no evidence here that a more elaborate down-weighted combination is needed, but this wasn't tested against alternatives (e.g., a orthologue-specific discount factor) and could be revisited if the primary-target gain needs to be made more robust.
