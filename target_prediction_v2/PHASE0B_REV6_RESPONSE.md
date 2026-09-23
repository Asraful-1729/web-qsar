# Response to Rev 6: H4 Decomposition Re-run Under Corrected Leakage Control and Widened Grid

**Answers Rev 6's three flags** (§1.1-1.2, decisions #11-13). **Code:** `phase0b/capture_scaffold_strict.py`, `phase0b/analyze_h4_rev6.py`. **Raw results:** `phase0b/results/phase0b_h4_rev6_report.json`.

**Headline: Rev 6's own executive summary needs a correction too.** Rev 6 stated H4a ("pooling beats best-similarity") was "large, robust" and separable from H4b. Once flag-3 is actually resolved (not just flagged), **that isn't uniformly true — H4a holds for Set B and does not hold for Set A under the stricter, more correct leakage control.** And H4b's earlier "significant negative" finding does not survive the wider grid: it is now genuinely unresolved, not leaning negative. This is not a failure of Rev 6's reasoning — flagging the confound before trusting the numbers was exactly right — it's what that check was for.

---

## Addendum: k and α grids widened further, rerun — all 20 fold selections unchanged

Following a direct request to widen both grids, the new question was whether k=100 (selected in 3 of 5 outer folds for Set A under scaffold-strict leakage, in the run above) was a real optimum or another boundary artifact like α=2.0 turned out to be. This required recapturing with a larger neighbour cap (150 → 300 stored neighbours per query — `capture.py`/`capture_scaffold_strict.py`, since k cannot be meaningfully tested past however many neighbours were actually stored), then rerunning nested CV over k ∈ {5,...,300} and α ∈ {0,...,24} (`analyze_h4_widegrid.py`, raw results `phase0b/results/phase0b_h4_widegrid_report.json`).

**Result: all 20 fold selections (5 outer folds × 2 sets × 2 leakage variants) are identical to the narrower-grid run**, k for k, α for α, to the decimal. k=100 is a genuine optimum, not a ceiling effect — nothing between 100 and 300 ever scored better on inner-fold MRR. α's selections (3–8 across folds) stayed just as far from their new ceiling (24) as before. Every downstream metric and BCa interval in this document is therefore unchanged and now **confirmed grid-robust** rather than provisional on an under-tested search range. No further grid widening is warranted from this evidence — the H4a/H4b conclusions above stand as the stable answer, not as an artifact of where the search happened to stop looking.

---

## Flag-1 (selection-metric consistency) — resolved, no new compute needed

Nested CV was tuned on **any-annotated MRR**. At that exact tier, weighted beats unweighted cleanly (Set A: 0.598 vs 0.571; Set B: 0.769 vs 0.761). The negative result reported in the original addendum was on the **primary-target** tier, never part of the selection loop. Not a contradiction — a tier-generalization gap: optimizing α for the broad/easy task does not transfer to the narrow/hard one, and (as flag-2's results below show) may not even be the right conclusion once the grid is right-sized.

## Flag-2 (grid too narrow) — resolved, and it changes the reading

Widened grid: α ∈ {0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12}. **The new edge-of-grid unanimity from Rev 5/original G0b (α=2.0 in all 10 folds) does NOT repeat.** Selected α now lands in the **interior** (mostly 3–8) across both sets and both leakage settings — Set A: α∈{6,8}; Set B: α∈{3,4,6}. None hit the new edge (12). This confirms Rev 6's diagnosis was correct: the old grid really was too narrow, and there is a genuine, non-monotonic optimum somewhere in the 3–8 range, not "weighting wants to go to infinity."

**New issue surfaced, not previously visible:** with the wider α grid, k now sometimes hits **its own grid edge (100)** — 3 of 5 outer folds for Set A under scaffold-strict leakage selected k=100. The k-grid (max 100) needs widening too before the pooling-size question is fully resolved. Flagging this now rather than letting it hide behind the α question that just got fixed.

## Flag-3 (leakage control, blocking) — resolved; this is the important part

Re-ran the full H4a/H4b decomposition under `capture_scaffold_strict.py` (near-duplicate removal **plus** the query's entire Bemis–Murcko scaffold group stripped from the reference), matching the stricter control that separately showed a −4/−7 point effect. Comparing side by side with the same wide grid under the original near-duplicate-only control:

### H4a: does pooling (unweighted 10-NN) beat best_similarity?

| Set | Leakage control | Tier | Δ Top-1 | 95% CI | Significant? |
|---|---|---|---:|---:|:---:|
| A | near-dup only | any-annotated | +0.035 | [−0.019, 0.086] | **No** |
| A | **scaffold-strict** | any-annotated | +0.025 | [−0.034, 0.080] | **No** |
| A | near-dup only | primary | +0.156 | [0.063, 0.275] | Yes |
| A | **scaffold-strict** | primary | +0.125 | [−0.015, 0.246] | **No — loses significance** |
| B | near-dup only | any-annotated | +0.085 | [0.036, 0.141] | Yes |
| B | **scaffold-strict** | any-annotated | +0.075 | [0.020, 0.133] | **Yes — holds** |
| B | near-dup only | primary | +0.160 | [0.000, 0.280] | borderline |
| B | **scaffold-strict** | primary | +0.200 | [0.020, 0.340] | **Yes — strengthens** |

**H4a does not survive uniformly.** It is robust for Set B (clinical) under both leakage settings — if anything it gets *stronger* under the stricter control. For Set A (approved), it was never significant on the any-annotated tier at all, and its one significant result (primary-target, near-dup-only) **disappears under the correct, stricter leakage control.** Rev 6's framing of H4a as "large, robust" was itself built on a comparison (weighted-vs-best-similarity, not unweighted-vs-best-similarity) that conflates the two effects it was trying to separate — once genuinely isolated, "pooling helps" is a Set-B finding, not a general one.

### H4b: does similarity-weighting beat unweighted pooling?

| Set | Leakage control | Tier | Δ MRR (weighted − unweighted) | 95% CI | Significant? |
|---|---|---|---:|---:|:---:|
| A | near-dup only | primary | −0.021 | [−0.063, 0.011] | **No** (was significant in the narrow-grid run: −0.022, [−0.053,−0.005]) |
| A | **scaffold-strict** | primary | **+0.040** | [−0.020, 0.106] | No |
| B | near-dup only | primary | −0.020 | [−0.080, 0.000] | No (borderline in narrow-grid run) |
| B | **scaffold-strict** | primary | −0.013 | [−0.072, 0.020] | No |

**The original "weighting is significantly worse on primary MRR" finding does not replicate.** Under the wider grid the point estimate for Set A even flips sign (−0.021 → +0.040), and none of the four cells reach significance. The honest reading is unchanged from Rev 6's own instinct — **H4b is genuinely unresolved** — but the specific negative claim from the first addendum should be **retracted**, not softened. It was an artifact of forcing α to the edge of a too-narrow grid, exactly as flag-2 suspected.

---

## Revised bottom line

| Question | Status after this rerun |
|---|---|
| H4a (pooling > best-similarity) | **Set-dependent, not general.** Robust for Set B under both leakage settings. Not statistically supported for Set A under the correct (scaffold-strict) leakage control. |
| H4b (weighting > unweighted) | **Unresolved**, and the earlier negative point estimate does not replicate under a wider grid — retract the "weighting hurts MRR" claim from the previous addendum. |
| Which leakage control should Phase 1 trust? | **Scaffold-strict**, per Rev 6 §2.5 — already the plan; this rerun is the confirmation that switching control changes conclusions enough to matter, not just magnitudes. |
| k-grid | **Also too narrow** — new finding, not previously visible. Needs widening (past 100) before H4a/H4b can be called settled even on Set B. |

**Why Set A and Set B disagree is still not diagnosed.** This is the same open question as Rev 5/G0b's "clinical > approved" puzzle (P0b-11), and this rerun makes it sharper: it's not just that Set B scores higher, it's that the *entire pooling-vs-singleton effect* is only detectable on Set B. Recommend folding this directly into Phase 1's P0b-11 diagnostic (max-similarity distribution, |T|, document holdout behavior, A vs B) rather than opening a fourth H4 sub-question — the mechanism is likely the same one already flagged (Set B's ChEMBL coverage is denser/more series-like), and one diagnostic should explain both puzzles or neither.

**Practical recommendation for Phase 1 (revising Rev 6 §2.1-2.3 slightly):** do not yet commit to "unweighted k-neighbour pooling" as a universal baseline-to-beat claim. Use the **internally-reproduced, fixed k=10 unweighted 10-NN** as the reference baseline (well-defined, simple, already computed on both sets) for G2's comparator, exactly as Rev 5 originally specified — but drop any claim that pooling's advantage over best-similarity is itself an established, general result. That claim is real for Set B and unestablished for Set A, and the benchmark card (Rev 6 §2.4's own standard) should say so explicitly rather than generalize from the set where it happens to hold.
