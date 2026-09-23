# Target Prediction v2 — Baselines Consolidated

**`BUILD_PLAN.md` Phase 1 item 5.** Every baseline listed in rev5 §5's own list ("frozen v1; popularity; internally-reproduced 10-NN; [1]'s L-score; external servers where a common universe can be built"), pulled together from where each was actually measured (Phase 0 and Phase 0b), with nothing re-run here — this is a consolidation pass, not a new experiment. See `PHASE1_BENCHMARK_CARD.md` for the evaluation-set/ground-truth context these numbers sit inside.

---

## 1. Frozen v1 (best_similarity, production default) — the floor every candidate must clear

Leakage-controlled (Tanimoto ≥0.95, near-duplicate removal), one search per query, threshold=0.0. Source: `PHASE0_AUTOPSY.md` §4.

| | Set A (n=200) | Set B (n=200) | Set D (n=334, full) | Set C (reused v1 benchmark) |
|---|---:|---:|---:|---:|
| Any-annotated Top-1 | 46.0% [39–53%] | 70.0% [63–76%] | 66.5% [61–71%] | 56.4% |
| Any-annotated Top-10 | 75.0% [69–80%] | 90.5% [86–94%] | 86.8% [83–90%] | 92.2% |
| Any-annotated MRR | 0.559 | 0.768 | 0.737 | 0.705 |
| Primary-target Top-1 | 48.4% (n=64) | 64.9% (n=57) | 58.2% (n=141) | — |
| Primary-target Top-10 | 79.7% | 94.7% | 89.4% | — |
| Primary-target MRR | 0.600 | 0.758 | 0.695 | — |

**Note, carried from Phase 0 and unresolved since**: Set B beats Set A on every metric — the opposite of the "drugs are harder than screening compounds" expectation. Plausible cause (dense clinical-program analogue series vs. sparser older approved-drug chemotypes), not re-tested since Phase 0's first observation.

## 2. Internally-reproduced weighted k-NN (H4) — the actual baseline to beat, not `best_similarity`

`score(t) = Σ over k nearest neighbours i whose annotation set contains t of tanimoto_i^α`. Source: `PHASE0_AUTOPSY.md` §7.

| | best_similarity (frozen) | k=10, α=1.0 | k=25, α=1.0 | k=50, α=1.0 |
|---|---:|---:|---:|---:|
| Set A primary-target Top-1 (n=64) | 48.4% | 60.9% | 62.5% | **65.6%** |
| Set A primary-target MRR | 0.600 | 0.707 | 0.709 | **0.719** |
| Set A % never found | 0.0% | 7.8% | 7.8% | 6.2% |
| Set B primary-target Top-1 (n=57) | 64.9% | **80.7%** (k=10, α=0.5/1.0) | 79.0% | 79.0% |
| Set B primary-target MRR | 0.758 | 0.852 | 0.846 | 0.846 |

**This is the single strongest, cleanest result across Phase 0/0b — the +12 to +17pt Top-1 gain is concentrated on the primary-target task, costs a real 6-8% coverage gap (best_similarity always returns a rank; weighted k-NN can return "not found" when no top-k neighbour carries the target), and any-annotated metrics barely move.** Confirmed not to be a Set A/B artifact but a density effect (`PHASE0B_DENSITY_STRATIFICATION.md`) and refined into the density-adaptive rule at n=800 scale (`PHASE1_DENSITY_ADAPTIVE_RULE.md`): use pooling above density 12, best-similarity below. **This is what "beat the baseline" means going forward — best_similarity alone is not a serious comparator once weighted k-NN exists.**

## 3. Popularity floor (sanity check, not a real comparator)

Rank targets by raw global popularity, ignoring the query molecule. Two independent measurements, consistent:

| Source | n | Any-annotated Top-1 | Any-annotated Top-10 | Any-annotated MRR |
|---|---:|---:|---:|---:|
| `PHASE0_AUTOPSY.md` §9 | 734 (A+B+D) | 6.5% | 30.5% | 0.163 |
| `PHASE0B_ADDENDUM.md` §9 | 734 (A+B+D) | 6.5% | 30.1% | 0.160 |

Primary-target tier (n=262, `PHASE0_AUTOPSY.md` §9): Top-1 3.4%, Top-10 14.9%, MRR 0.086. Every v1 configuration tested clears this floor by a wide margin on every metric — confirms the method does real structure-based work, not popular-target recovery. **Not informative beyond this sanity check — do not use as a G-gate comparator.**

## 4. L-score reliability baseline ([1]'s own comparator — Phase 4's bar to beat, not Phase 3's)

`l` = count of a query's 10 nearest neighbours annotated to the predicted target; `L = l/10`, used as a confidence/reliability score. Source: `PHASE0B_ADDENDUM.md` §7.

| L | n predictions (Set A) | Mean precision (ours) | [1]'s published figure |
|---:|---:|---:|---:|
| 0.1 | 1,343 | 0.176 | 0.191 |
| 0.4 | 69 | 0.522 | 0.505 |
| 0.7 | 32 | 0.594 | 0.714 |
| 1.0 | 44 | 0.886 | 0.929 |

Same monotone shape as [1]'s table; our numbers track slightly below at the top end on Set A (closer on B/D, see raw JSON). **This is the specific number `BUILD_PLAN.md` §9 gate G3 binds to: Phase 4's calibration layer must beat L-score on Brier score and ECE, or ship L-score instead — no complexity credit without earning it against this exact table.**

## 5. Potency-floor variants (H3) — a tunable hyperparameter, not a separate baseline

Set A, n=200, reduced reference restricted to `pChEMBL ≥ floor`. Source: `PHASE0_AUTOPSY.md` §6.

| | No floor (index as-is) | ≥10 µM floor | ≥1 µM floor |
|---|---:|---:|---:|
| Any-annotated Top-10 | 75.0% | 71.5% | 68.0% |
| Primary-target Top-1 (n=64) | 48.4% | 54.7% | 51.6% |
| Primary-target MRR | 0.600 | 0.642 | 0.627 |

Two-directional: helps primary-target ranking (removes noisy weak evidence), hurts any-annotated breadth-recall (filters out real, weaker-potency annotated off-targets along with noise). **Ship as a tunable hyperparameter (Phase 2 already plans this), not a fixed cutoff.**

## 6. Leakage-control sensitivity — confirms baselines above aren't leakage artifacts

| Level | Result | Source |
|---|---|---|
| Near-duplicate cutoff (0.90 vs. 0.95 Tanimoto) | No material difference — every gap within sampling noise at n=31-100 | `PHASE0_AUTOPSY.md` §8 |
| Scaffold-group removal (whole analogue series stripped, not just near-dupes) | **Material**: Top-1 −4pts, Top-10 −7pts vs. near-dup-only control (n=100, Set A) | `PHASE0B_ADDENDUM.md` §8 |

**All baseline numbers in §1-§2 above use near-duplicate-only leakage control (Phase 0's original protocol), not the scaffold-level control.** The scaffold-level drop is real and material — any Phase 3 candidate comparison against these baseline numbers must use the *same* leakage protocol on both sides (scaffold-safe, per `PHASE1_HOLDOUTS_AND_SPLITS.md`'s partition), not compare a scaffold-safe candidate against a near-dup-only baseline number, which would understate the candidate's true improvement.

## 7. External servers (SwissTargetPrediction, SEA) — not run, disclosed as a real gap

**Never executed in this program.** Both lack a public bulk API; a comparison run would require a small manual/scripted per-compound query on a common query set. Per `PHASE0_AUTOPSY.md` §11 and rev5's own phase boundary, this is explicitly **Gate G5's job** (a later phase), not a Phase 0/1 requirement — flagged here so its absence from this consolidated table is a documented decision, not an oversight.

---

## Summary — what "beating the baseline" means at each future gate

| Gate | Baseline to clear | Table |
|---|---|---|
| G2 (retrieval/ranking candidates) | Weighted k-NN / density-adaptive rule (§2), **not** raw best_similarity | §2, `PHASE1_DENSITY_ADAPTIVE_RULE.md` |
| G3 (calibration) | L-score (§4) on Brier + ECE | §4 |
| G5 (external validity) | SwissTargetPrediction/SEA — **not yet measured**, still open | §7 |
| Sanity-only, any gate | Popularity floor (§3) — should never be close | §3 |
