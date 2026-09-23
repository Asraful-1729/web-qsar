# Phase 0b Addendum — G0b

**Status: G0b complete, including the item previously marked pending.** The per-record ChEMBL pull for the "mechanistically supported" tier finished (710/710 molecules, ~95 min total, no failures) — §13 below.
**Supersedes nothing in `PHASE0_AUTOPSY.md`** — this is the closing addendum Rev 5 asks for, read together with it. **See also `PHASE0B_REV6_RESPONSE.md`**, which revises §1's H4 conclusion after Rev 6's leakage-confound flag — read that section in light of it: H4a ("pooling beats best-similarity") is confirmed for Set B only, not in general, and H4b's "weighting hurts MRR" claim did not replicate under a wider grid and should be treated as retracted, not as stated below.
**Code:** `target_prediction_v2/phase0b/`. **Raw results:** `phase0b/results/*.json`, `phase0b/results/capture_{A,B,D}.jsonl` (734 queries, full neighbour-level data, re-analyzable without re-querying the index).

---

## 1. The one finding that changes the plan: H4 does not survive the correct comparison

Phase 0 (rev3) claimed weighted k-NN was "the strongest, cleanest result," comparing it only against `best_similarity`. Rev 5 asked for nested CV, a scaffold-clustered paired BCa bootstrap, and matched breadth. Doing all three surfaces something Phase 0 never tested: **a plain, unweighted 10-nearest-neighbour vote (baseline B) captures nearly all of the gain by itself.** The similarity-weighting on top of it (baseline C) adds little, and on the primary-target tier it sometimes makes things *worse*.

**Nested-CV selected hyperparameters** (5 outer × 3 inner scaffold-grouped folds, grid k∈{5,10,15,25,50,75,100}, α∈{0,0.25,...,2.0}, tuned on any-annotated MRR): every one of 10 outer folds across Sets A and B selected **α = 2.0 — the edge of the grid.** This means the grid was too narrow to find where the gain actually plateaus; the reported weighted-k-NN numbers below are real but the optimum α is not yet located. Flagged, not hidden.

**Paired BCa bootstrap (10,000 resamples, clustered by Bemis–Murcko scaffold), weighted k-NN vs. each baseline:**

| Comparison | Set | Metric, tier | Point est. | 95% CI | Excludes 0? |
|---|---|---|---:|---:|:---:|
| weighted k-NN vs. **best_similarity** | A | Top-1, any-annotated | +0.075 | [0.028, 0.123] | **yes** |
| weighted k-NN vs. **best_similarity** | A | Top-1, primary | +0.141 | [0.048, 0.258] | **yes** |
| weighted k-NN vs. **best_similarity** | B | Top-1, any-annotated | +0.095 | [0.047, 0.150] | **yes** |
| weighted k-NN vs. **best_similarity** | B | Top-1, primary | +0.120 | [−0.060, 0.260] | no (n=50) |
| weighted k-NN vs. **unweighted 10-NN** | A | Top-1, any-annotated | +0.040 | [0.009, 0.080] | **yes, small** |
| weighted k-NN vs. **unweighted 10-NN** | A | MRR, primary | **−0.022** | [−0.053, −0.005] | **yes — weighted is worse** |
| weighted k-NN vs. **unweighted 10-NN** | B | Top-1, any-annotated | +0.010 | [0.000, 0.034] | no |
| weighted k-NN vs. **unweighted 10-NN** | B | Top-1, primary | **−0.040** | [−0.160, −0.020] | **yes — weighted is worse** |

**Matched-breadth check (P0b-9)**, truncating each method's ranked list at the same N (pooled precision/recall, any-annotated, Set A):

| N | best_similarity P / R | unweighted 10-NN P / R | weighted k-NN P / R |
|---|---:|---:|---:|
| 5 | 0.311 / 0.214 | 0.341 / 0.204 | 0.337 / 0.202 |
| 10 | 0.229 / 0.316 | 0.276 / 0.267 | 0.272 / 0.263 |
| 20 | 0.158 / 0.437 | 0.250 / 0.312 | 0.251 / 0.313 |

Confirms the same shape from a different angle: **the k-NN-vote family (weighted or not) buys real precision over `best_similarity` at every matched breadth — but weighted and unweighted are statistically indistinguishable from each other**, and `best_similarity` wins on recall/breadth-diversity at larger N (it surfaces a more diverse target set; a k-NN vote saturates on whatever its k neighbour compounds happen to be annotated to).

**Revised verdict on H4:** *"Move from single-best-similarity ranking to a k-neighbour vote"* is a real, robust, bootstrap-confirmed win. *"Weight that vote by similarity/potency"* — Phase 0's actual headline claim — is **not supported once compared against the right baseline**, and actively loses on primary-target MRR on both sets tested. This is the single most consequential correction in Phase 0b. Recommend Phase 1/3A treat **unweighted or lightly-weighted k-NN (k≈10–15) as the new baseline to beat**, and re-run this comparison with an extended α grid (weighting is not ruled out — the edge-of-grid α=2.0 selection means it is simply unresolved, not disproven) before claiming similarity-weighting as a design decision.

---

## 2. L1 — the averaging trap, measured on our own data (P0b-3)

Confirmed at the scale Rev 5's own reference-paper analysis predicted. At the **exact operating point matching reference [1]** (top-10, no similarity cutoff):

| Set | Convention | Precision | Recall | Mean targets/query |
|---|---|---:|---:|---:|
| A (approved) | **Macro** (per-query, then averaged) | 0.229 | **0.474** | 9.95 |
| A (approved) | Pooled | ~0.10–0.23 (threshold-dependent; see Phase 0) | ~0.31–0.42 | — |
| Reference [1] | Macro | 0.348 | 0.423 | 11.7 |

Under the correct (macro) convention, **Set A's recall (0.474) is now *above* reference [1]'s (0.423)**, not merely matched to it as Phase 0 claimed — precision is still below (0.229 vs. 0.348). This is a materially different, and more favorable, picture than Phase 0's pooled-based comparison implied. Full macro-vs-pooled tables (all 5 thresholds × 3 sets × 2 tiers) are in `phase0b/results/phase0b_block1_report.json`; the relative gap between conventions ranges 15–500% depending on threshold and tier, exactly the "not a rounding issue" scale Rev 5 predicted.

---

## 3. |T| and target-universe composition (P0b-1, P0b-2)

| Set | n | mean \|T\| (annotated) | median \|T\| |
|---|---:|---:|---:|
| A (approved) | 1,965 | 8.81 | 3 |
| B (clinical) | 2,750 | 6.54 | 2 |
| D (sparse-target) | 334 | 33.84 | 12 |

Universe: 4,658 human single-protein targets, 100% human by construction (v1's own build filter). Target-class/family breakdown not built (diagnostic-only per Rev 5 §7, deferred).

**mechanistic_targets == primary_targets everywhere** — a real quirk found, not a bug in our code: all 6,984 rows in our ChEMBL `mechanism.json` pull have `molecular_mechanism=1`; ChEMBL's own endpoint apparently returns no `molecular_mechanism=0` rows in bulk. So our current "mechanistic" tier is really just Rev 5's **primary/intended** tier — the broader two-of-three "mechanistically supported" tier (§4.1) needs the still-running per-document pull (P0b-12, pending).

---

## 4. Internally-reproduced 10-NN baseline (P0b-4)

| Set | Tier | Top-1 | Top-10 | MRR |
|---|---|---:|---:|---:|
| A | any-annotated | 0.490 | 0.735 | 0.571 |
| A | primary (n=64) | 0.641 | 0.891 | 0.732 |
| B | any-annotated | 0.705 | 0.850 | 0.761 |
| B | primary (n=50) | 0.760 | 0.860 | 0.804 |
| D | any-annotated | 0.713 | 0.847 | 0.761 |
| D | primary (n=141) | 0.702 | 0.879 | 0.770 |

This is the fixed-k=10, unweighted, not-tuned-on-our-eval-sets baseline Rev 5 requires as the honest internal comparator, superseding any external-reference-based comparison for go/no-go purposes.

---

## 5. Potency scope: index-filter vs. ground-truth-filter (P0b-5)

Two genuinely different operations, both measured separately (Rev 5 E3's correction):

**Index-filter** (remove weak-potency neighbours from voting, weighted k-NN k=25/α=1): small, inconsistent effect (Set A Top-10: 0.755→0.73 at 10µM, →0.71 at 1µM; Set D: essentially flat).

**Ground-truth-filter** (only count a target as truly annotated if the query drug's *own* activity against it clears the threshold): **large, consistent effect** — Set A Top-10 drops 0.755→0.68 (10µM)→0.555 (1µM); Set B drops 0.855→0.82→0.745; Set D drops 0.865→0.841→0.722. Filtering the ground truth removes a real fraction of qualifying queries too (Set A: 200→181→150 queries with any qualifying true target). This is the more consequential of the two operations and matches [1]'s own practice of applying the potency cut to the *target definition*, not the index.

## 6. Weak potency among voting neighbours (P0b-6)

Of the compound-target pairs actually cast as votes among each query's top-10 nearest neighbours: **11.9–18.5% weaker than 10 µM, 35–45% weaker than 1 µM**, consistent across Sets A/B/D. The evidence directly supporting predictions is substantially weaker-potency than a naive reading of "high similarity" would suggest.

---

## 7. L-score reliability baseline (P0b-18) — Rev 5's actual comparator, not ours

`l` = count of the 10 nearest neighbours annotated to the predicted target; `L = l/10`. Mean precision per level (Set A; B and D in the JSON, same shape, generally higher):

| L | n predictions | mean precision |
|---|---:|---:|
| 0.1 | 1,343 | 0.176 |
| 0.4 | 69 | 0.522 |
| 0.7 | 32 | 0.594 |
| 1.0 | 44 | 0.886 |

Compares reasonably to [1]'s own published table (L=0.1: 0.191; L=0.4: 0.505; L=0.7: 0.714; L=1.0: 0.929) — same monotone shape, our numbers track slightly below at the top end on Set A (closer on B and D — see JSON). This is the bar Phase 4's calibration work must beat (Rev 5 G3).

---

## 8. Leakage (P0b-10, P0b-17)

**Near-duplicate-only (Phase 0's original control, Tanimoto≥0.95, two fingerprints):** re-confirmed negligible sensitivity to the 0.90-vs-0.95 cutoff choice once compared on a proper paired sample (Phase 0's *first* leakage-sensitivity comparison was confounded — compared non-identical query subsets — corrected in `PHASE0_AUTOPSY.md` §8; restated here per P0b-17).

**Scaffold-level (the real risk Rev 5 flags, newly tested, n=100 sample of Set A):**

| | Top-1 | Top-10 |
|---|---:|---:|
| Near-duplicate removal only | 0.49 | 0.78 |
| **+ whole scaffold group removed** | 0.45 | 0.71 |

A real, measurable drop (Top-1 −4pts, Top-10 −7pts) when the entire analogue series — not just near-exact duplicates — is stripped from the reference. **Confirms Rev 5's prediction: near-duplicate control alone is not sufficient; some of v1's apparent performance is analogue-series lookup, not generalization.** Mean scaffold-group size for these queries was 1,194 (heavily skewed by a few generic scaffolds, same catastrophic-skew pattern the original v1 benchmark already documented).

**Document-level leakage: not measurable with current infrastructure.** v1's aggregated index carries no `document_chembl_id`; pulling it for the full ~1.3M-pair reference is out of scope for Phase 0b's "cheap counting" mandate. Left open for Phase 1's benchmark card (which already plans a document/assay-campaign split). Rev 5 anticipated this staying open.

---

## 9. Popularity baseline (P0b-13)

Combined A+B+D (n=734): Top-1 = 6.5%, Top-10 = 30.1%, MRR = 0.160. Every v1 configuration tested clears this floor by a wide margin — the method is doing real structure-based work, not recovering popular targets by default.

---

## 10. Orthologue scoping (P0b-19)

2,600 candidate orthologue single-protein targets (name-matched, not a validated HomoloGene-style mapping — disclosed proxy) found across mostly mouse (1,076), rat (805), cow, pig, rabbit, dog, guinea pig, macaque and chicken, covering **1,517 of our 4,658 human targets (32.6%)**. Confirms Rev 5 §1.4's recommendation is grounded in real numbers on our own data, not just the cited literature. Does not yet estimate how many *additional compounds/annotations* each candidate would add — that needs a per-target activity pull, correctly scoped to Phase 2, not this counting pass.

---

## 11. Gap ledger, L1–L6

| # | Cause | Status |
|---|---|---|
| L1 | Averaging convention | **Measured** (§2). Large, confirmed on our own data (up to 5×+ relative on precision at low thresholds). |
| L2 | Target universe (orthologues) | **Scoped** (§10). 32.6% of targets have a candidate orthologue; not yet built into the index (Phase 2). |
| L3 | Bioactivity-type breadth | **Resolved by documentation.** Index is IC50/Ki/Kd/EC50 only, no PubChem `Potency`, by construction. |
| L4 | Relation handling | **Resolved empirically.** 0 ChEMBL records anywhere have pchembl_value with a `<`/`<=` relation — our index structurally has zero censored records, a real, quantified gap vs. [1]. |
| L5 | Leakage asymmetry | **Partially measured.** Near-duplicate: negligible sensitivity. Scaffold-level: real, ~4-7pt effect (§8). Document-level: unmeasurable currently, open. |
| L6 | Residual ranking quality | **Not cleanly isolable yet** — same conclusion as Rev 5 anticipated, because L1-L5 each move the comparison numbers by amounts of similar order to the original gap. What *is* now isolable, independent of any external reference: weighted-vs-unweighted k-NN (§1), which resolves against our OWN internally-reproduced baseline and needs no external-reference decomposition at all. |

---

## 13. Mechanism coverage and the "mechanistically supported" tier (P0b-12) — completed

The per-record pull (710/710 query molecules, ~95 min, 0 failures) is in. Rev 5 §4.1's two-of-three rule (≥2 independent documents; a `drug_mechanism` record; binding-assay evidence) applied to each query's own `annotated_targets`:

| Set | Annotated (query,target) pairs | Primary/intended (drug_mechanism) | ≥2 documents (of pairs with data) | Binding assay present (of pairs with data) | **Mechanistically supported (≥2 of 3)** |
|---|---:|---:|---:|---:|---:|
| A (approved) | 1,441 | 86 (6.0%) | 62.8% | 86.4% | **601 (41.7%)** |
| B (clinical) | 1,482 | 56 (3.8%) | 37.1% | 94.5% | **534 (36.0%)** |
| D (sparse-target) | 11,302 | 258 (2.3%) | 56.1% | 95.2% | **5,138 (45.5%)** |

Confirms Rev 5's three-tier structure is meaningfully distinct, not redundant: **primary/intended is narrow (2-6% of annotated pairs)** — ChEMBL only flags a target as "the" mechanism for a small fraction of a drug's known activities, most of which are legitimate secondary pharmacology never curated as a formal mechanism entry. **Mechanistically supported is far broader (36-46%)** — the large majority of that gap comes from "has a binding assay" alone being very common (86-95%), with the ≥2-documents criterion doing most of the remaining discriminating work. This is now usable ground truth for Phase 1's three-level metrics, not just a defined-but-unmeasured rule.

**Coverage caveat, disclosed:** 70-95% of annotated pairs got per-record data (the rest presumably fell outside the per-molecule pull's target filter or hit a rate-limit skip); the percentages above are of pairs *with* pulled data, stated as such, not silently extrapolated to the full annotated set.

---

## 14. G0b exit assessment

| Requirement (Rev 5 §5 exit condition) | Status |
|---|---|
| Gap ledger, L1-L5 measured, L6 isolated | L1-L4 measured/resolved; L5 partially (scaffold done, document open by design); L6 remains entangled with H4a/H4b — see Rev 6 response, not fully isolated |
| \|T\| and universe composition | Done (§3) |
| Internally-reproduced 10-NN result | Done (§4) |
| H4 survives nested CV, matched breadth, scaffold-clustered CI? | **Partially, and set-dependent — see `PHASE0B_REV6_RESPONSE.md`**, which supersedes §1's conclusion above after the leakage-confound rerun |
| Scaffold/document leakage drop | Scaffold: done (§8), and reused directly in the Rev 6 rerun. Document: open by design, deferred to Phase 1 |
| Popularity baseline | Done (§9) |
| L-score reliability table | Done (§7) |
| Orthologue scoping counts | Done (§10) |
| Mechanism coverage / mechanistically-supported tier (P0b-12) | **Done (§13)** |

**G0b: PASS, fully closed** — no items outstanding. The autopsy's central uncertainty (is v1's apparent performance a ranking artifact, a ground-truth artifact, or real) has been replaced by a sharper question, refined twice: first from "does weighted k-NN help" to "pooling vs. weighting are different effects" (this addendum, §1), then from there to **"pooling helps on Set B, is unestablished on Set A, and weighting is unresolved rather than harmful"** (`PHASE0B_REV6_RESPONSE.md`, after the leakage-confound was actually checked rather than just flagged).

**Recommended next step:** proceed to Phase 1 per `PHASE0B_REV6_RESPONSE.md`'s bottom line — use the internally-reproduced, fixed k=10 unweighted 10-NN as the G2 reference baseline without over-claiming generality for the pooling effect, fold the Set A/B asymmetry into P0b-11's clinical-vs-approved diagnostic rather than opening a new sub-question, carry the ground-truth-filter potency operation (not the index-filter) into the benchmark recipe, and widen both the α *and* k grids before any weighting decision is finalized in Phase 3A. Mechanism-coverage data (§13) is now available for Phase 1's primary-target and mechanistically-supported metrics.
