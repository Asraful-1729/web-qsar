# Diagnosis: Why Set B (clinical) Outperforms Set A (approved)

**Answers P0b-11** (Rev 5 §5, Block 3), open since Phase 0 first noticed the effect and sharpened by the Rev 6 rerun (the entire H4a pooling-vs-singleton effect turned out to hold for Set B and not Set A under correct leakage control — this is why the diagnosis, not just the observation).
**Code:** `phase0b/diagnose_A_vs_B.py`. **Raw results:** `phase0b/results/diagnose_A_vs_B_report.json`. Uses the already-captured neighbour data (`capture_{A,B}.jsonl`, widened to 300 neighbours/query) plus one direct lookup against the full v1 index and the completed mechanism-support pull — no new ChEMBL queries.

---

## The answer: Set B has denser structural neighbourhoods, not richer documentation

| | Set A (approved) | Set B (clinical) | Direction |
|---|---:|---:|---|
| Top-10 recovery (any-annotated) | 74.5% | **86.5%** | B higher |
| Top-1 recovery | 45.5% | **62.0%** | B higher |
| max similarity to index (mean / median) | 0.677 / 0.700 | **0.712 / 0.750** | B higher |
| neighbours at Tanimoto ≥0.4 (mean) | 73.0 | **91.6** (+25%) | B denser |
| neighbours at Tanimoto ≥0.5 (mean) | 31.7 | **48.6** (+53%) | B denser |
| neighbours at Tanimoto ≥0.7 (mean) | 3.28 | **9.00** (+174%) | B much denser |
| \|T\| (annotated targets, mean / median) | 7.2 / 3.0 | 7.4 / **2.0** | B has fewer typical targets to find |
| documents per annotated target (mean) | **4.59** | 2.63 | **A higher** — ruled out below |
| own scaffold-group size (mean, median) | 905.6 / 8.0 | 334.5 / 7.0 | A has a heavier tail of huge generic scaffolds |

**The wrong hypothesis, tested and ruled out:** if Set B won simply because it has *more total bioactivity data or literature confirmation*, documents-per-target should be higher for B. It's the opposite — Set A's annotated targets have **more** independent documents on average (4.59 vs 2.63). Whatever is driving the gap, it isn't "B is just better-studied."

**The right explanation:** Set B's queries sit in **structurally denser neighbourhoods** — far more close analogs already in the reference pool, especially at high similarity (nearly 3× as many neighbours ≥0.7 Tanimoto). This is a similarity-search method, so what matters isn't how much has been published about a target in general, it's how many near-identical compounds already sit next to the query in chemical space.

## Why this is causal, not coincidental: within-set validation

The same features were checked as predictors of success *within* each set (recovered-top-10 vs. not, same set) — if a feature predicts success inside both A and B with the same sign, a between-set difference in that feature is a legitimate causal candidate, not a spurious correlation:

| Feature | Set A: recovered vs. not | Set B: recovered vs. not | Consistent? |
|---|---:|---:|:---:|
| max similarity to index | 0.708 vs 0.586 | 0.736 vs 0.562 | **Yes, both large** |
| neighbours ≥0.4 | 83.7 vs 41.9 | 98.3 vs 48.2 | **Yes, ~2× in both** |
| neighbours ≥0.5 | 38.8 vs 11.0 | 53.2 vs 18.9 | **Yes, ~3× in both** |
| mean target reference depth | 2364 vs 1490 | 2648 vs 2654 | No — matters in A, not in B |
| own scaffold-group size | 865 vs 1023 (smaller wins) | 286 vs 643 (smaller wins) | Yes, same direction, smaller/more-specific scaffolds do better |

Similarity/neighbour-density features predict recovery cleanly in **both** sets, with the same sign and comparable magnitude — exactly the signature of a real mechanism, not noise. Target reference depth is revealing in a different way: it strongly separates recovered/not-recovered **within Set A** but is essentially flat **within Set B** (2648 vs 2654) — i.e., for Set B, target-side evidence is already abundant enough everywhere that it stops being the bottleneck; for Set A it's often the limiting factor. That is consistent with, not contradictory to, the main finding: Set A's problem is chemical-neighbourhood sparsity for *some* queries specifically, not a uniform data shortage.

## What this means for the plan

1. **The H4a set-dependence (Rev 6) now has a mechanism, not just an observation.** Pooling votes across k neighbours only helps when there are enough genuinely close neighbours to pool — in a sparse neighbourhood (typical of Set A), pooling averages in noise from weaker matches instead of reinforcing signal from strong ones. In a dense neighbourhood (typical of Set B), pooling has real signal to average over. **Confirmed directly — see `PHASE0B_DENSITY_STRATIFICATION.md`: splitting Set A's own 200 queries by neighbourhood density alone reproduces the Set A/B split almost exactly, entirely within Set A.** "Set A vs Set B" was a proxy for density, not the real variable.
2. **This is likely the same mechanism as H1's original "benchmark too easy" question and H5's abstention design.** A density-based stratum (already informally present as the "reference-evidence bucket" in `V1_FREEZE.md`) may be the right unifying variable for several previously-separate observations (Set A/B gap, H1's bucket skew, H5's threshold-sweep coverage tradeoff) rather than three different phenomena.
3. **Practical, near-term recommendation:** report Set A and Set B results *and* a neighbourhood-density stratification together in the Phase 1 benchmark card, rather than as two unrelated per-set numbers. If density (not drug-development phase per se) turns out to be the real driver, that's a more useful, generalizable finding than "clinical beats approved," and it directly motivates D5 (explicit abstention) as a targeted fix for the sparse-neighbourhood regime rather than a blanket floor.
4. **Not tested here:** true document-level leakage/holdout behavior (the other half of P0b-11's original ask) — the documents-per-target comparison above is a related but different measurement (evidence richness, not leakage). Document-level leakage remains open per the existing gap-ledger entry (L5), unchanged by this diagnosis.
