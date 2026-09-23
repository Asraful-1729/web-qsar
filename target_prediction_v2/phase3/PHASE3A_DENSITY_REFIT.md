# Phase 3A — Density-Adaptive Rule Re-Fit Against the Rebuilt Index

**Status: DONE.** This closes the mandatory prerequisite `BUILD_PLAN.md` Phase 3A's own text stated: the original threshold (density≥12, `PHASE1_DENSITY_ADAPTIVE_RULE.md`) was fit against v1's frozen 1,312,849-pair index; Phase 2's rebuild roughly tripled native-human pair volume (3,973,438 pairs) via `Potency`/censored/inactive record inclusion, which was flagged at the time as certain to shift the density distribution enough to invalidate that threshold. It has — re-fit below, same methodology, real numbers.

---

## What changed to make this possible

Reused the *exact* fitting methodology (`phase0b/fit_density_adaptive_rule.py`: isotonic regression + scaffold-clustered bootstrap, plus a logistic mixing-weight approximation) — nothing new decided here, only re-run against new data. Getting to that point required building a parallel infrastructure stack, since the density-adaptive fit depends on the full retrieval pipeline, not just the aggregated pairs:

1. **`phase2/build_v2_index_files.py`** — converts Phase 2's `stage7_native_human_pairs.jsonl` + `stage6_fingerprints.npz` into v1's exact index file format (`v2_index/compounds.csv.gz` + `fingerprints.npz`), so the already-validated `target_fishing.py`/`leakage.py` search pipeline could be reused unchanged rather than rebuilt from scratch. 3,973,438 pairs, 1,347,013 distinct compounds.
2. **`phase2/build_v2_leakage_fp2.py`** — a second, independent fingerprint (RDKit topological, 2048-bit) per distinct compound, mirroring `phase0/build_leakage_fp2.py` exactly, for near-duplicate leakage detection under two independent fingerprint families.
3. **`phase3/capture_scaffold_strict_v2.py`** — `phase0b/capture_scaffold_strict.py`'s exact logic (near-dup removal under both fingerprints, then whole-scaffold-group removal, then neighbour capture), pointed at the new index.
4. **`phase3/fit_density_adaptive_rule_v2.py`** — the re-fit itself.

## Two real bugs found and fixed during this build (not visible at planning time)

1. **Query-smiles resolution mismatch.** `eval_sets.json`'s query drug SMILES were matched against v1's (buggy `SaltRemover`-based) standardization; Phase 2 uses the fixed `rdMolStandardize.FragmentParent()` (item 5's fix). Confirmed empirically before writing the capture script: 88.9%/93.4% of Set A/B query drugs are present in Phase 2's data at all (the rest have no activity data in the Phase 2 pull's scope — a real, disclosed coverage change, not a bug); of those present, 98%+ have identical resolved SMILES. Fixed by resolving each query's current standardized SMILES via `stage6_std_cache.json` rather than trusting `eval_sets.json`'s stale field.
2. **Index-membership mismatch, found only when Set B crashed mid-run.** Checking `stage6_std_cache.json` alone isn't sufficient — a compound can have a valid standardized SMILES there but zero surviving native-human pairs after Stage 7's single-protein-target filtering, making it absent from `v2_index/compounds.csv.gz` specifically. Fixed by cross-checking against the *actual* index's distinct-SMILES set before including a query drug (coverage dropped slightly further, e.g. Set B 93.4%→93.0%, once this stricter, correct check was applied).

## The re-fit result

**n=1,600 pooled queries (Set A + Set B, n=800 each), 1,295 distinct scaffold groups.**

| Density | Fitted benefit | 95% CI (scaffold-clustered bootstrap, 2000 resamples) | Significant? |
|---:|---:|---|:---:|
| 0 | −0.044 | [−0.083, −0.010] | **Yes — significantly negative** |
| 1 | −0.044 | [−0.083, −0.009] | **Yes — significantly negative** |
| 2 | −0.009 | [−0.061, 0.024] | No |
| 3 | −0.004 | [−0.041, 0.044] | No |
| 5 | −0.004 | [−0.038, 0.046] | No |
| **8** | **0.051** | **[0.032, 0.065]** | **Yes — first significant crossing** |
| 12 | 0.051 | [0.032, 0.065] | Yes |
| 18–70 | 0.051 | [0.032–0.034, 0.065–0.069] | Yes |
| 100–300 | 0.054 | [0.037–0.040, 0.094–0.114] | Yes |

**New rule: pool (unweighted k-NN) when density ≥ 8; use best-similarity below that** — down from the previous density≥12 threshold, and holds significant at every higher grid point through density=300 (`threshold_holds_for_all_higher_grid_points: true`).

## What's different from the original v1-index fit, and why it matters

1. **Threshold moved down (12→8).** A denser reference pool (3x more native-human pairs) means a given raw neighbour count corresponds to a different effective evidence level than before — expected direction, now confirmed with real numbers rather than assumed.
2. **Effect size at plateau is larger** (~0.051–0.054 vs. the original fit's ~0.031) — pooling's benefit in the dense regime is *stronger* against the richer index, not weaker. Plausible mechanism: more reference evidence per dense query means the pooled vote aggregates over a larger, more reliable neighbour set.
3. **A genuinely new finding, not clearly present before**: at very low density (0–1 neighbours), pooling now shows a **significant negative** effect (CI excludes zero on the negative side), not just "no significant benefit." The original v1-index fit's low-density band was non-significant in both directions; here there's real signal that pooling actively hurts when evidence is nearly absent. Worth carrying forward as an explicit abstention/floor consideration for Phase 3A's negative-evidence and D5 abstention-policy work, not just noting a null zone.

## Logistic mixing-weight approximation

Collapsed to a near-hard step again (as the original fit did): `s_scale≈0.0`, `d0≈5.25`, floor=−0.026, ceiling=0.052 — consistent with the isotonic curve's own sharp jump between density 5 and 8. Same conclusion as before: the practical rule is closer to a hard threshold than a smooth sigmoid, and the logistic form's `s_scale≈0` is itself informative (not a fitting failure).

## What this does NOT change

- This re-fit used the **same tier (any-annotated), same k=10 unweighted pooling vs. best-similarity comparison** as the original — no new hypothesis introduced, purely a re-measurement of the already-established comparison against new data.
- Popularity correction, potency/confidence weighting, the orthologue arm, and the negative-evidence term (all still-open Phase 3A scope items per `BUILD_PLAN.md` §5) are unaffected by this re-fit and remain to be built.
- The **density≥8 threshold from this re-fit is now the number to use for any Phase 3A work going forward** — supersedes `PHASE1_DENSITY_ADAPTIVE_RULE.md`'s density≥12 for anything touching the rebuilt index.
