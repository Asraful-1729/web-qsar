# Phase 3B — SEA (Keiser et al. 2007) Comparison

**Status: DONE.** A genuine SEA reimplementation, a real bug found and fixed in its null-distribution model, and a final, sanity-checked comparison against v2 and v1.

## What was built

`build_sea_null_distributions_v2.py` + `test_sea_vs_v2_fixed.py`. Raw score: `sum(Tanimoto(query, ligand))` over every distinct ligand indexed against a target (uncapped, no threshold — same "active set" definition used everywhere else in this program). Raw scores are only meaningful once converted to a significance value against a target-size-specific null distribution — the same logic as a BLAST E-value — fit by Monte Carlo simulation (random query vs random same-size ligand-set draws from the real reference pool), per Keiser et al.'s original design.

## A real bug found, diagnosed, and fixed

**First attempt** (`build_sea_null_distributions.py`, discrete size bins, top bin `1280-999999`) gave an implausible result: SEA scored *catastrophically* worse than v1's crude `best_similarity` (any-tier effect −0.538, CI excludes zero). That contradicts everything known about SEA's real-world performance, so it was treated as a suspected implementation bug and not written into any result document until diagnosed.

**Diagnosis**: inspecting individual true-positive (query, annotated-target) pairs found z-scores up to **636–641** — targets with **106,614** and **114,949** ligands were being scored against a null distribution whose top bin used a fixed "representative size" of only **~1,640** (computed as the bin's midpoint, capped at 2000). A direct check of the real index (`phase2/v2_index/compounds.csv.gz`) confirmed the scale of the problem: ligand-set sizes range up to **195,809**; **63 targets exceed 10,000 ligands, 15 exceed 50,000**. Scoring a 100K+-ligand target against a null fit for ~1,640 ligands produces an astronomically over-significant p-value regardless of true relevance — these few mega-promiscuous targets would then dominate SEA's ranked output for nearly every query, explaining the catastrophic result.

**Fix** (`build_sea_null_distributions_v2.py`): replaced the discrete, capped bins with a **continuous power-law null model**. Null mean and std were sampled at 16 log-spaced representative sizes spanning the full real range (7 to 196,608 ligands, 80 null queries × 15 subsets each), then fit as power laws of size:

- `mean(N) = 0.11812 * N^1.0001` (log-log R² = 1.00000)
- `std(N)  = 0.01645 * N^0.9813` (log-log R² = 0.99968)

Both fits are essentially exact — `mean` scales linearly with target size as expected for a sum over N similarity terms; `std` also scales nearly linearly rather than as `sqrt(N)`, consistent with positive correlation among reference-compound similarities (chemical-space clustering), not independence. Gumbel `loc`/`scale` are derived from `mean`/`std` via the standard moment relations (`scale = std / (π/√6)`, `loc = mean − 0.5772×scale`) for any real target size, in-range.

**Validation of the fix**: re-inspecting the same true-positive pairs post-fix, z-scores dropped to single digits (mostly ±3, one legitimate high-confidence hit at z≈19, one degenerate single-atom-fragment query at z≈−8) — no more scale-driven outliers.

## Final result (n=150, 75/set, fresh-holdout population, seed=999 — same query set as the original buggy run, so before/after is a clean comparison of the null-model fix alone)

| Comparison | Effect (paired reciprocal-rank, bootstrap) | 95% CI | Significant? |
|---|---:|---|:---:|
| v2 (full pipeline) vs. SEA, any-target | **+0.399** | [0.333, 0.467] | **Yes** |
| v2 (full pipeline) vs. SEA, primary-target (n=54) | **+0.569** | [0.464, 0.670] | **Yes** |
| SEA vs. v1 `best_similarity`, any-target | **−0.371** | [−0.442, −0.300] | **Yes** |

Fixing the null-model bug moved the "SEA vs v1" effect from −0.538 to −0.371 — a real, substantial correction — but the direction did not flip: **SEA still significantly underperforms v1's simple best-similarity search on this benchmark**, even with a correctly-calibrated null model.

## Honest interpretation — not a v2-favoring artifact

This is very plausibly a genuine property of the evaluation setup, not a remaining bug:
- The benchmark task is "does the model recover the single documented target for a query with known ChEMBL activity," and — because the reference pool is large and shares substantial chemical space with any given query — a close structural analog of the true target commonly exists in the pool. `best_similarity` (and v2, which pools weighted nearest-neighbour votes) is built to exploit exactly that: a single very close match wins immediately.
- SEA's ensemble raw-score-sum, even correctly null-normalized, dilutes a single strong analog's signal across the full (often large) ligand set of a target — a design intended for robustness to weak, diffuse evidence across many imperfect analogs, not for maximizing precision when one near-exact match already exists.
- In other words: this benchmark's structure (near-neighbour-rich ChEMBL data) favors nearest-neighbour-style methods over ensemble-normalization methods — a known trade-off in the target-prediction literature, not a defect in either implementation.

## Disclosed limitations of this comparison

- This SEA reimplementation is a good-faith one, not a byte-for-byte reproduction of Keiser et al.'s original tuning (bin/size choices, exact `MIN_LIGANDS=5` cutoff, and Monte Carlo sample counts are this program's own disclosed choices).
- Evaluated on the fresh-holdout population (n=150 subset), not the locked scaffold test (already fully spent on the v1-vs-v2 G2 result — cannot be reused).
- The original (buggy, discrete-bin) result and null file are retained as `sea_null_distributions.json` / `sea_vs_v2_result.json` for transparency but are **superseded** by `sea_null_distributions_v2.json` / `sea_vs_v2_result_fixed.json` — only the fixed numbers should be cited.

*Sources: `build_sea_null_distributions_v2.py`, `test_sea_vs_v2_fixed.py`, `results/sea_null_distributions_v2.json`, `results/sea_vs_v2_result_fixed.json`.*
