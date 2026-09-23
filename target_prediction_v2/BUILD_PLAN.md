# Target Prediction v2 — Final Build Plan

**Synthesizes:** procedure revisions 3, 5, 6, 7 (rev4 superseded, folded into rev5's errata) and all Phase 0 / Phase 0b empirical work (`PHASE0_AUTOPSY.md`, `PHASE0B_ADDENDUM.md`, `PHASE0B_REV6_RESPONSE.md`, `PHASE0B_DIAGNOSIS_A_VS_B.md`, `PHASE0B_DENSITY_STRATIFICATION.md`, `PHASE0B_REV7_RESPONSE.md`).
**Status:** G0b closed. This document is the authoritative plan going into Phase 1 — it does not re-derive findings already established; it cites them and says what to build.
**v1 stays frozen** (`target_fishing_v1_freeze/`) throughout everything below.

---

## 1. Objective (unchanged since rev3)

General target-identification performance — precision, recall, ranking, coverage, calibration, behaviour on unfamiliar chemistry. Not a single accuracy number. Ceiling expectation, grounded in [1]: a well-tuned ligand-centric method reaches roughly 0.35 macro precision / 0.42 macro recall on approved drugs at ~11.7 predicted targets/query over 4,167 targets. Approach or modestly exceed it; make reliability and incompleteness visible; don't oversell.

---

## 2. What's actually known now (do not re-litigate these)

| # | Finding | Evidence | Status |
|---|---|---|---|
| H1 | Original v1 benchmark was moderately (not drastically) easier than drug-realistic queries | 72–84% vs. 98% in the "51+ evidence" bucket [Phase 0] | Closed |
| H3 | No potency floor in v1's index; ground-truth-filtering it is large, index-filtering it is small/mixed | 8.4%/29.2% of pairs weaker than 10µM/1µM; ground-truth filter drops Top-10 recovery 8–14pts, index filter ≈flat [G0b] | Closed — use ground-truth filter (10µM), not index filter, for comparability with [1] |
| L1 | Macro vs. pooled averaging differs by up to ~500% relative; at [1]'s own operating point, our macro recall (0.474) *exceeds* [1]'s (0.423) | [G0b §2] | Measured. **Not fully closed** — Rev 6 flagged this reversal may be a ground-truth-*scope* artifact (broadening the definition mechanically converts unlabelled predictions to true positives) rather than a modelling win. Must be stated as such in the benchmark card, not marketed as "beats the reference." |
| L2 | Orthologue tier is a real, sizeable opportunity | 2,600 candidate targets (mouse/rat/cow/pig/rabbit/dog/guinea pig/macaque/chicken), covering 32.6% of our 4,658 human targets [G0b] | Scoped, green light, **not yet built** |
| L3 | Index is IC50/Ki/Kd/EC50 only, no PubChem `Potency` | Documentation | Closed |
| L4 | Index structurally has **zero** censored (`<`/`<=`) records | ChEMBL never populates `pchembl_value` for non-`=` relations — verified empirically, 0 counter-examples across the full API | Closed. **Technical consequence for Phase 2**: including `<` relations requires computing potency from raw `standard_value`/`standard_relation`, not from ChEMBL's precomputed `pchembl_value` field. |
| L5 | Near-duplicate leakage control (Tanimoto≥0.95) ≈ no sensitivity to 0.90 vs 0.95. Scaffold-level leakage is real: −4pts Top-1, −7pts Top-10 when the whole analogue series is removed | Paired comparisons [Phase 0, G0b] | Near-dup: closed. Scaffold: closed, **now mandatory control**. Document-level: still unmeasurable, open |
| H4a | "Pooling beats best-similarity" is **not general** — holds for Set B under both leakage settings (gets *stronger* under scaffold-strict); does not hold for Set A under scaffold-strict control | Paired BCa, both leakage settings [Rev 6 response] | Closed as set-dependent |
| — | **The real variable is neighbourhood density, not drug phase — CONFIRMED AT SCALE (n=800/set).** Dense half beats best-similarity in both sets (Set A Top-1 +0.05 CI[0.02,0.08]; Set B +0.063 CI[0.03,0.10]); sparse half shows no significant benefit in either (not negative — genuinely flat) | Within-set stratification, n=200 exploratory then n=800 confirmatory [`PHASE1_DENSITY_CONFIRMATORY.md`] | **Closed — confirmed, with corrections.** Effect sizes in the original 200-query pass were ~2x too large (regression to the mean on replication); H4a is present in **both** sets, not "Set B only" — Set A's is just smaller because less of its mass sits in the dense regime |
| — | ~~Dose-response is non-monotonic: a density "valley" (Q2) where pooling *significantly hurts*~~ | Quartile BCa at n=200 | **Retracted at n=800.** The Q2 band still points negative but is no longer statistically distinguishable from zero — exactly the small-sample artifact Rev 7 itself flagged as a risk when reviewing the original result |
| — | Depth and density are NOT interchangeable — target reference depth stays predictive within every density decile, including the two highest | Within-decile stratified check, n=800 [`PHASE1_DENSITY_CONFIRMATORY.md` §6] | Closed — retain both as separate stratification variables, don't collapse to density alone |
| — | Similarity-spread as a driver of the (now-retracted) valley effect | Preliminary, n=37/half then n=40/half at n=800 | Still genuinely open — carry into Phase 3A as planned |
| H4b | Similarity-weighting (`s^α`) on top of pooling: **no significant effect, either direction, at any tested width** | Grid widened twice (α to 24, k to 300); all 20 fold selections identical to the narrower grid — k=100 and α∈[3,8] are real interior optima, not boundary artifacts [Rev 6 response, Rev 7 response] | **Closed as a genuine null.** Do not spend further engineering time searching the `s^α` family or the Itskowitz–Tropsha variational alternative — the effect being searched for has been shown not to exist in the tested range. |
| — | Compound-level primary-annotation coverage: 28.7% (A) / 34.6%(A, full-pop, corrected)/ 28.7%(B) / 42.2%(D) — **below the ~60% adequacy bar** for primary-target metrics on all sets. Mechanistically-supported tier, **now full-population for A/B (Phase 1 item 7, superseding the sample estimate below)**: A 50.4%, B 57.2%, D 74.3% (already full-pop) — clears the bar clearly for D only | Full-population query for both tiers as of Phase 1 [`PHASE1_BENCHMARK_CARD.md` §5] | Closed, revised. Use mechanistically-supported as the powered ground-truth tier for **D only**; report A's and B's primary/mechanistic metrics with an under-coverage caveat (B's sample estimate of 63.5% did not hold at full-population scale) |
| — | L-score reliability baseline (Peón et al.'s vote-fraction score) computed on all 3 sets, tracks [1]'s published shape | [G0b §7] | Closed — this is Phase 4's comparator to beat |
| — | Popularity baseline cleared by wide margin in every configuration tested | [G0b §9] | Closed — v1 (any variant) does real structure-based work |
| — | Internally-reproduced, fixed k=10 unweighted 10-NN computed on all 3 sets | [G0b §4] | Closed — this is the honest internal comparator pending the density-adaptive rule |

---

## 3. Architecture to build

The retrieval core is **not** "best-similarity" and **not** "a single globally-tuned weighted k-NN." It is a **density-conditioned retrieval rule**, the direct consequence of §2's central finding:

```
Query SMILES
   |
Standardisation / QC, salt-parent collapse
   |
Domain check: per-query density vector
   (neighbour counts at Tanimoto >=0.4/0.5/0.6/0.7, computed post-leakage-removal)
   |
Density-adaptive retrieval  <-- Phase 1 fits this; Phase 0b only located the regimes
   - sparse / Q2 "valley" regime -> best-similarity (or near-k=1), NOT pooling
   - dense (Q3+) regime          -> unweighted (or near-uniform) k-neighbour pooling
   - functional form: fitted continuous curve (decile-resolution or spline),
     not the hand-set quartile thresholds Phase 0b used to find the shape
   |
Popularity correction (observed/expected vote mass) -- untested lever, still planned
   |
Potency weight (smooth, index-level) + confidence weight (index-level,
   only meaningful once Phase 2 relaxes the confidence>=8 filter -- see below)
   |
SEA branch (targets with >=5 ligands only) -- combine only if it beats
   tuned k-NN at matched breadth
   |
Orthologue-derived candidates -- separate reported stratum, never merged silently
   |
Calibrated scorer -> prediction set | explicit abstention
   (abstention driven by the density signal, not a fixed 0.4 floor)
   |
Reliability + evidence panel + incompleteness range
   |
Grouped output: individual targets + families ("high-evidence predicted
   targets", never "primary" unless from drug_mechanism)
   |
Final report
```

**Explicitly NOT building:** a similarity-weighting search or a data-derived weighting function (H4b closed as null — §2). A hard potency or confidence filter inside the retrieval index (H3 — weight, don't filter). A single fixed k or α applied everywhere (superseded by the density-adaptive rule).

---

## 4. Ground-truth recipe (final, for Phase 1)

Adopts rev5 §4.2, with the technical corrections §2 forced:

- ChEMBL release stated by version.
- Activity ≤10µM (ground truth), relations `=` and `<`/`<=` **computed from raw `standard_value`+`standard_relation`, not from ChEMBL's `pchembl_value`** (L4's consequence), taking the most potent value per pair.
- Bioactivity types: IC50, Ki, EC50, Kd, **+ ChEMBL `Potency`-type records** (closes L3; **corrected** per Phase 2 item 3 — this is not narrowly PubChem-sourced as originally described, predominantly `src_id=1` Scientific Literature; ~2.98M records at confidence≥8/human, roughly doubling index volume — see Phase 2 item 3 for the quality filters and the resulting density-adaptive-rule re-fit dependency in Phase 3A).
- Assay types: binding **and** functional (confirm explicitly against the current pull, not assumed).
- **Corrected** (Phase 2 item 2's decision, §5): confidence ≥7 was rev5's original spec, but ChEMBL's confidence_score is not a continuum below 8 — 7/6 assign activity to a **protein complex** (a different target entity than any single-protein target in our universe), 5/4 are explicitly target-ambiguous. **Threshold stays ≥8** (matches v1's current, already-correct-for-a-different-reason practice), with `confidence_score` itself retained per pair and used as a **two-level weight** (full weight @9 direct, downweighted @8 homology-assumed) — not a widened filter.
- Exclude non-specific/multi-protein targets by name and type.
- Deduplicate compound-target pairs; **collapse salt/parent forms to one compound** (rev5 §4.1's caution: ChEMBL mechanism/indication annotations already propagate across salt/parent, so failing to collapse is a leakage vector, not just a counting nuisance).
- **Retain, per compound-target pair: `n_distinct_documents`, `assay_types_seen`, `best_relation`.** This is a concrete schema change from v1's aggregated index — its absence is exactly what forced the slow, query-scoped side-pull for mechanism-support data in Phase 0b. Build it in from the start this time.
- Run twice: human-only and human+orthologue tier, reported as **separate strata**, never silently merged.
- Three ground-truth levels stored at build time, not computed ad hoc:
  - **Annotated**: any qualifying activity row.
  - **Mechanistically supported**: ≥2 of {≥2 independent documents, `drug_mechanism` record, binding-assay evidence} — rev5's 2-of-3 rule, now with real coverage numbers (§2) to scope expectations.
  - **Primary/intended**: from `drug_mechanism` directly. Use for Sets B/D headline primary-target metrics; report Set A's as under-powered (§2).

---

## 5. Phases

### Phase 1 — Benchmark rebuild, lockbox, entry requirements (weeks 2–4, revised up from 2–3 given the added density work)

**Standard Phase 1 scope (rev5 §6.1, unchanged):**
1. **DONE.** Locked scaffold test built (795 compounds, 598 scaffold groups, seed=42, sealed) — `phase1/build_holdouts.py`, `phase1/data/holdouts.json`. Renewable temporal holdout: **BLOCKED**, no per-pair date retained in v1's index — deferred to Phase 2's rebuild (§ `PHASE1_BENCHMARK_CARD.md` §7).
2. **DONE.** Set C verified drug-*like* but not literal real drugs: 0.42% direct overlap with Sets A/B (21/5,001 compounds), 84.7% Lipinski pass rate on a 2,000-compound sample. Acceptable for its stated purpose, described precisely going forward (`PHASE1_BENCHMARK_CARD.md` §8).
3. **Partially done.** Random and scaffold split types built for the full 3,975-drug tuning/validation population (`phase1/data/holdouts.json`). Document/assay-campaign and temporal split types: **BLOCKED**, same data gap as item 1 — deferred to Phase 2.
4. **DONE.** `PHASE1_BENCHMARK_CARD.md` — full recipe, `|T|` distributions, universe composition, averaging convention, all four leakage levels' status, evaluation-set sizes, Set C finding, holdout/split numbers.
5. **DONE.** `PHASE1_BASELINES.md` — frozen v1, popularity floor, internally-reproduced weighted k-NN (H4), L-score reliability baseline, potency-floor variants, leakage sensitivity, all consolidated. External servers (SwissTargetPrediction/SEA): confirmed **not run**, explicitly deferred to G5 per rev5's own phase boundary, not a Phase 1 gap.
6. **Prepared, not run** — genuinely cannot be completed autonomously. Candidate worksheet built and ready (153 top-ranked false positives, 54 drugs, `phase1/data/adjudication_candidates.csv`) via `phase1/build_adjudication_candidates.py`; the actual adjudication (literature/database lookup, human judgment) needs a real domain reviewer (rev5's own estimate: 2-3 days). See `PHASE1_ADJUDICATION_STUDY_STATUS.md` — this is the one Phase 1 item this program cannot close by itself.
7. **DONE.** Full-population mechanistically-supported coverage recomputed: A 50.4%, B 57.2% (down from a 63.5% sample estimate — B no longer clearly clears the ~60% bar at scale), D 74.3% (already full-population). Primary/intended-tier coverage: 31.2% overall, full-population from the start. See `PHASE1_BENCHMARK_CARD.md` §5 — **this revises the Phase 0b recommendation**: mechanistically-supported is now the adequately-powered tier for D only, not B.

**New Phase 1 entry requirements, from Rev 7, gating Phase 3A specifically (not the rest of Phase 1) — all six closed:**
8. **DONE.** `PHASE1_DENSITY_CONFIRMATORY.md` — n=800/set confirmatory re-test, decile resolution attempted (found underpowered), Q2 valley claim retracted with corrected effect sizes.
9. **DONE.** `PHASE1_REMAINING_ITEMS_9_12_13.md` §Item 9 — properly powered (n=871/729 vs. the original n=37), dilution hypothesis explicitly not supported.
10. **DONE.** Non-linear (within-decile) joint density+depth check in `PHASE1_DENSITY_CONFIRMATORY.md` §6 — depth remains predictive even within top deciles (D9 gap=1737, D10 gap=719).
11. **DONE.** `PHASE1_DENSITY_ADAPTIVE_RULE.md` — isotonic + scaffold-clustered bootstrap fit; rule: pool when density ≥12, else best-similarity.
12. **DONE.** `PHASE1_REMAINING_ITEMS_9_12_13.md` §Item 12 — density vs. max-similarity Spearman ρ=0.827 (partial unification); density does not subsume depth.
13. **DONE.** `PHASE1_REMAINING_ITEMS_9_12_13.md` §Item 13 — real σ_diff=0.2481, ICC≈0.279, design effect=1.074; recommendation: G2 margin Δ≥0.02.

### Phase 2 — Data rebuild (weeks 4–5)

**STATUS: EXECUTED AND FROZEN.** All seven items below were not just decided but physically run against ChEMBL_37 (pinned 2026-09-22) — full staged execution log in `PHASE2_REBUILD_EXECUTION_PLAN.md`, frozen deliverables in `V2_FREEZE.md`. Final numbers: 1,590,210 distinct standardized compounds, 3,973,653 native-human (compound,target) pairs, 34,080 orthologue-tier pairs, 385 human targets with a confirmed orthologue. Two additional bugs were found and fixed only during execution, not anticipated in planning — both disclosed in full in `V2_FREEZE.md`: (a) `activity.json` never serializes `confidence_score` in its own output (backfilled via a separate `assay.json` lookup, 270,225 assay ids resolved); (b) `confidence_score__gte=8` alone does not reliably exclude non-single-protein targets — 18.99% of pairs (931,457/4,905,110) referenced an invalid target (e.g. `CHEMBL372`, `target_type=ORGANISM`) before the `single_protein_target_ids.json` cross-reference was actually applied at aggregation time; fixed, verified 0 remaining. **The density-adaptive rule (density≥12) must be re-fit against this rebuilt index before Phase 3A uses it** — this was already flagged as a consequence of item 3 below, now actually due.

1. **DECIDED (scope).** Potency as a graded index weight. Per §5's own Phase 3A scope line ("potency/confidence weighting (index-level)" is explicitly listed there, not here), **Phase 2's job is the data commitment, not the weight function** — fixing a formula now without evidence would repeat the mistake this program has refused everywhere else (density-adaptive rule was fit, not guessed). Settled now: (a) **10µM ground-truth membership threshold is unchanged** (§4, still a hard boundary — a pair must clear it to count as "annotated" at all); (b) **retain the actual measured potency (pChEMBL, computed from raw `standard_value`+`standard_relation` per §4's L4 correction, not ChEMBL's own `pchembl_value` field) per pair**, not just the binary ≤10µM gate — this is the concrete schema addition, alongside item 5's other per-pair fields. Deferred to Phase 3A, explicitly: fitting the graded weight function itself, via the identical isotonic-regression + scaffold-clustered-bootstrap methodology already built for the density-adaptive rule (`phase0b/fit_density_adaptive_rule.py`), against real paired-difference data. Grounded in why this lever matters at all: H3 (`PHASE0_AUTOPSY.md` §6) showed a hard potency cutoff is a two-directional tradeoff — helps primary-target Top-1 (+6.3pts) but hurts any-annotated recall (−3.5pts) by discarding real weak-potency annotations outright; a graded weight is structurally better than either extreme because it lets weak evidence still count while counting for less. Real scale of what's being weighted, already measured on our own data (`PHASE0B_ADDENDUM.md` §6): 11.9–18.5% of pairs actually cast as votes among top-10 neighbours are weaker than 10µM, 35–45% weaker than 1µM — not a marginal edge case.
2. **DECIDED.** Confidence as an index weight — **not** by relaxing the fetch filter to admit confidence 4-7 data. Confirmed live against ChEMBL's own confidence_description values (not assumed): confidence 4-9 is not a continuum of "same evidence, decreasing certainty" for a single-protein-target benchmark. Only 9 ("Direct single protein target assigned," 375,124 human assays) and 8 ("Homologous single protein target assigned," 1,764) are single-protein-target assignments at all. 7/6 ("Direct/Homologous protein complex subunits assigned," 13,649/31) assign activity to a **protein complex**, a different target entity than any single-protein `target_chembl_id` in our 4,658-target universe. 5/4 ("Multiple direct/homologous protein targets **may be** assigned," 23,346/84) are explicitly target-ambiguous — ChEMBL itself won't commit to one target. Relaxing to 4-7 would be a no-op (those records mostly don't map onto our existing single-protein target IDs) or would require an unprincipled complex→subunit or multi-target→single-target inference, contaminating ground truth the way H3/scaffold-leakage/mechanism-support work has consistently refused to. **Revised spec**: a two-level weight within data already in the index (v1 already hard-filters to ≥8, so both 9 and 8 are already included, just currently treated identically) — full weight for confidence=9, downweighted for confidence=8. No fetch-filter change. Requires `confidence_score` itself retained per pair (add to item 5's field list below).
3. **DECIDED, with real scale/quality findings.** Four sub-items, checked against live ChEMBL data and v1's actual current code (not assumed):
   - **Censored `<`/`<=` relations**: real volume, ~174K records (165,645 `<` + 8,589 `<=`) vs. 2.6M `=`-relation baseline (+6.7%). Safe for ground-truth membership (a `≤10µM` record is at least that potent). Retain an explicit `is_censored` flag per pair — a censored value is a lower bound, not exact, and item 1's graded potency weight must not silently treat it as an exact measurement (would understate weight for the assay's most potent hits, which are disproportionately the ones that hit a detection limit).
   - **`Potency` records (closes L3)**: much larger than the plan's original phrasing assumed — ~2.98M records at confidence≥8/human, roughly **doubling** current index volume, not a small gap-fill. Real quality issues found: `assay_type` includes non-binding/functional codes (e.g. `T`) with null values (filter to `assay_type in {B,F}` + non-null `standard_value`); source is predominantly `src_id=1` Scientific Literature, not narrowly PubChem as the plan's rationale text claimed (correct that framing). Must reuse v1's existing `single_protein_target_ids.json` cross-reference (`build_target_fishing_index.py`'s own docstring: confidence≥8 alone does not guarantee single-protein target assignment) rather than re-deriving it. **Consequence that must not be silently skipped**: this volume roughly doubles the index, which will materially shift the density distribution Phase 1's density-adaptive rule (density≥12 threshold, `PHASE1_DENSITY_ADAPTIVE_RULE.md`) was fit against — **that threshold must be re-fit in Phase 3A against the rebuilt index**, not carried forward unchanged.
   - **Functional assays**: **already confirmed true, no new work.** v1's actual fetch/build code (`build_target_fishing_index.py`) has no `assay_type` filter at all — functional (`assay_type='F'`) records are already included today whenever ChEMBL populated a `pchembl_value`. Settles the plan's own "confirm explicitly, not assumed" instruction.
   - **Most-potent value per pair**: real correction needed. v1's current code aggregates multiple measurements per (compound, target) pair by **mean** pChEMBL (`build_target_fishing_index.py` Pass 4: `mean_pchembl = round(sum(pchembl_vals)/len(pchembl_vals), 2)`), not max — contradicting both §4's recipe and this item's own stated intent. Every pChEMBL value used anywhere in this program to date (H3, density work, etc.) is mean-based. Phase 2 switches to max/most-potent aggregation — disclosed as a real behavior change; not expected to overturn density/ranking conclusions (those depend on neighbour-count structure, not exact pChEMBL magnitude) but changes every reported potency value going forward.
4. **DECIDED, with a real correction to the plan's own cited prior.** Orthologue tier, built as a first-class experiment: species-provenance column per row, orthologue-derived hits as a separate reported stratum (never merged into the human-only headline).
   - **Method upgrade over P0b-19**: use ChEMBL `target_synonym__iexact` gene-symbol matching (confirmed precise on 6 real targets below — exact gene-symbol hits for F2/F10/MMP9/CTSK/AR/NR3C1), not P0b-19's "name-matched, disclosed proxy" (`PHASE0B_ADDENDUM.md` §10). This is the same method validated during the Phase 1 adjudication pilot's orthologue-gap check (`phase1/adjudication_pilot.py`).
   - **Real incremental-coverage numbers, measured live** (not assumed) on 6 targets — new distinct compounds gained per target from non-human single-protein orthologues, over and above compounds already annotated to the human target: Prothrombin +18.8% (1,047 new/5,583 human), Coagulation factor X +2.6%, MMP-9 +1.6%, Cathepsin K +3.1%, Androgen receptor +25.2% (744 new/2,953 human), Glucocorticoid receptor +3.1%. **Confirms the tier is worth building** (single-digit-to-double-digit percent new compound coverage per target, not negligible) but with **high per-target variance that swamps any class-level pattern in this small sample**.
   - **Correction to [23]'s cited prior**: this 6-target sample shows the **opposite** direction from "expect proteases > nuclear receptors" — mean gain 6.5% (proteases) vs. 14.1% (nuclear receptors), driven by androgen receptor's outlier +25.2% (plausibly heavy rodent-model use in AR pharmacology, a real confound distinct from target class). **n=6 is nowhere near enough to overturn [23]** — but the plan should not state "expect proteases > nuclear receptors" as a given. Revised: **test the target-class-dependent-benefit hypothesis empirically at full scale (Phase 3A ablation), do not assume it** — and stratify by research-tradition/species-availability confounds, not target class alone, since this small sample suggests those may matter more.
5. **DECIDED, with a confirmed live bug found in v1's current (already-shipped) index, not just a Phase 2 design question.**
   - **v1's reference index already does salt-stripping** (`SaltRemover.StripMol()`, keyed by standardized SMILES, confirmed by reading `build_target_fishing_index.py`) — but this is **not sufficient**. Confirmed via a direct test on the exact case the Phase 1 adjudication pilot flagged as a "salt-form gap" (ladostigil tartrate, `PHASE1_ADJUDICATION_STUDY_STATUS.md`): `SaltRemover` correctly strips the tartrate counter-ion but cannot collapse **n:1 stoichiometry salts** (the parent fragment appears twice in a 2:1 salt) down to one canonical structure — the stripped result keeps two disconnected copies joined by `.`, which does not match the true single-copy parent SMILES. `rdMolStandardize.FragmentParent()` was tested and correctly resolves this exact case. **This directly explains at least one of the 15 salt-form gaps the adjudication study already found** — not a coincidence, a confirmed causal link between two independent findings.
   - **Real, measured scale of the live bug**: scanned v1's actual production index (`compounds.csv.gz`, 1,312,849 rows) directly — **2,167 rows (~0.17%) are un-collapsed multi-fragment artifacts** (e.g., the identical structure written twice joined by `.`), sitting as orphaned entries disconnected from their true parent's existing index row. Small relative to total volume, but real and live, not hypothetical. Separately, in our own `eval_sets.json` query-drug population (4,715 drugs): 43 (0.9%) have multi-fragment SMILES, 15 (0.3%) show the exact duplicate-fragment stoichiometry pattern.
   - **Fix**: switch standardization from `SaltRemover.StripMol()` to `rdMolStandardize.FragmentParent()` (RDKit's MolVS-derived standardization, designed for exactly this) for the Phase 2 rebuild's one standardisation pipeline.
   - **Provenance per row**: since collapsing genuinely merges multiple distinct original `molecule_chembl_id`s into one canonical row, retain the list of contributing original `molecule_chembl_id`s (not just a single winner) per collapsed row — needed for auditability and for any future re-derivation of per-molecule-ID facts (e.g., the mechanism-support pulls, which are scoped per `molecule_chembl_id`).
   - **Sparse targets kept and flagged**: unchanged policy carried forward from Phase 0's Set D framing — no new work, confirmed as still the right default (never drop a target for having few compounds; flag it, per `PHASE0_AUTOPSY.md` §4's Set D caveat).
   - Per-pair retained fields (§4): `n_documents`/`assay_types_seen`/`best_relation`/**`confidence_score`** (item 2)/**`pchembl_value`** (item 1, max/most-potent not mean — item 3's correction)/**`is_censored`** (item 3)/**`contributing_molecule_chembl_ids`** (this item).
6. **DECIDED (scope), real volume confirmed.** Matches the item-1/item-3 pattern: **Phase 2's job is the data commitment; the functional form (how measured-inactive evidence gets used in scoring, e.g. an explicit negative term) is explicitly Phase 3's job to test**, per this item's own text — not re-litigated here. Settled now:
   - **Data source, confirmed real and substantial via live query**: `>`/`>=` relation records with `standard_value ≥10µM` (symmetric with our existing active/ground-truth boundary, human, confidence≥8, our 4 bioactivity types) — **490,821 records** (487,802 `>` + 3,019 `>=`). These are exactly the records item 3 correctly excluded from ground-truth *positive* evidence (uncertain how much less potent than the bound) — here they become the right thing: genuine measured-and-tested negatives, not "never tested" compounds treated as presumed-inactive (the known bias PIDGIN's convention exists to avoid).
   - **Secondary signal, noted with a real data-quality caveat**: categorical `activity_comment` inactive flags (e.g. "Inactive"/"inactive") also exist, but casing is inconsistent across ChEMBL depositors — "Inactive" and "inactive" return materially different counts (64,582 vs. 103,518) under exact string match. If used, normalize case and treat as a secondary/supplementary signal to the numeric-relation approach above, not the primary source.
   - **Schema**: retain an `is_measured_inactive` flag per pair (alongside item 5's other per-pair fields), scoped through the same `single_protein_target_ids.json` cross-reference already established (confidence≥8 alone never guarantees single-protein target assignment — applies to any activity pull, not just Potency records).
   - Citation stays exactly as rev5 framed it: an adaptation of PIDGIN's established practice [19], never claimed as novel — Phase 2's contribution here is confirming the data actually exists at real, usable scale, not inventing the method.
7. **DECIDED (protocol), execution deferred — a real scope boundary, stated plainly.** This item is the freeze *protocol*, not the freeze itself: the actual rebuild (items 1-6's decisions, physically executed — pulling item 3's ~3M `Potency` + ~490K measured-inactive + ~174K censored records, re-running standardization with `FragmentParent`, etc.) is a large, separate engineering task not undertaken in this planning session, likely requiring ChEMBL's bulk database dump rather than REST API calls given the volume (item 3's Potency-record pull alone, at the ~1.8 molecules/sec rate achieved for the much smaller mechanism-support pull, would take unreasonably long via REST). Settled now, so freezing is mechanical once the rebuild completes:
   - **Reuse v1's proven pattern exactly** (`V1_FREEZE.md`, `CHECKSUMS.sha256`, `manifest.json`) — don't invent a new scheme. `CHECKSUMS.sha256` = per-file SHA256 of every artifact in the snapshot (index files + the new per-pair-field files item 5 added); `manifest.json` = build provenance (counts, fetch filters, timestamps); `V2_FREEZE.md` = the narrative doc (architecture, dataset table, frozen decisions, "what re-opens this freeze," regenerate/verify commands) — matching `V1_FREEZE.md`'s exact structure.
   - **Closes a real, long-standing gap**: this program has never pinned an exact ChEMBL release (`PHASE1_BENCHMARK_CARD.md` §1/§10's disclosed gap; v1's own `manifest.json` doesn't have one either — confirmed by reading it). Fixed mechanism, confirmed live: `GET https://www.ebi.ac.uk/chembl/api/data/status.json` returns `chembl_db_version` (`"ChEMBL_37"` as of this check) and `chembl_release_date` (`"2026-05-01"`) — Phase 2's rebuild fetch script must capture this at fetch time and record it in `manifest.json`. First time this program pins a release.
   - **"Gap ledger (closed)"** (§9's V2_FREEZE spec) = `BUILD_PLAN.md` §2's own findings table. At freeze time, confirm every row shows Closed, or an explicitly-still-open item is disclosed with its reason (matching this program's standing practice) — not a new ledger built from scratch.
   - Every later Phase 3/4 result cites the resulting version hash + ChEMBL release, matching how `V1_FREEZE.md` is already cited as this program's fixed reference point throughout Phase 0/0b/1.

### Phase 3 — Model bake-off (weeks 5–8)

- **3A — density-adaptive retrieval core, the highest-value item in the whole bake-off. STATUS: DONE.**
  - **Density re-fit** (`phase3/PHASE3A_DENSITY_REFIT.md`): new rule, **pool when density ≥ 8** (down from 12), significant through density=300, n=1,600/1,295 scaffold groups. Plateau effect larger than the original fit (~0.051-0.054 vs. ~0.031); new significant-negative low-density band (0-1) — carry into D5's abstention policy.
  - **Popularity correction** (`phase3/PHASE3A_POPULARITY_CORRECTION.md`): **tested, rejected.** Significantly negative at every eps tested (0.1-100), any-tier and primary-tier alike. Do not ship. D2's "de-flooded by popularity correction" framing needs a different mechanism.
  - **Potency weighting** (`phase3/PHASE3A_POTENCY_WEIGHTING.md`): **tested, adopted.** Significantly positive at every configuration tested. Recommended: threshold=6.0, floor=0.0, steepness=5.0 (any_all +0.043, any_dense +0.029, primary +0.055, all CIs exclude zero).
  - **Orthologue ablation** (`phase3/PHASE3A_ORTHOLOGUE_ABLATION.md`): **tested, mixed real result.** No significant effect on any-tier (even slightly negative in the dense regime), but a real, if modest, significant positive effect on primary-target recovery (+0.019, CI [0.001, 0.041], 49% query coverage). Adopt for primary-target scoring specifically, not blanket.
  - **Negative-evidence term** (`phase3/PHASE3A_NEGATIVE_EVIDENCE_STATUS.md`): **blocked, disclosed, not faked.** Current capture schema can't distinguish "no evidence" from "measured inactive" per neighbour-target link. Needs a schema extension + a third capture run + a new scoring function (none of which exist yet) — scoped as its own follow-up, not attempted here.
  - Confidence weighting (item 2's two-level 8-vs-9 weight): not yet tested — `v2_index`'s schema doesn't carry `confidence_score` per pair (only in Stage 7's raw output); needs the same kind of schema extension as the negative-evidence term. Not attempted in this pass.
  - **Do not** re-open the `s^α` weighting search (H4b closed) — potency weighting above is a distinct, legitimately-open lever, not a reopening of that closed question.
- **3B — SEA and combination. DONE.** A genuine SEA implementation was built (target-size-stratified Gumbel/EVD null distributions fit by Monte Carlo simulation) — `phase3/PHASE3B_SEA_COMPARISON.md`. A real bug was found and fixed along the way (the first null model broke down for mega-promiscuous targets up to 195,809 ligands, producing an implausible result; fixed with a continuous power-law null, re-validated). Gate condition does **not** pass: v2's density-adaptive core significantly beats SEA (any-tier +0.399, primary +0.569, both CI excluding zero, n=150) — SEA is not adopted into the production score, kept only as the published-method comparator it was built to be.
- **3C — multi-fingerprint ladder. Correctly not pursued**, per the plan's own de-prioritization below the orthologue tier (now closed, §3A) — no work needed here, this is compliance with the plan's own instruction, not a gap.
- **3D — neural challenger. NOT ATTEMPTED — disclosed scope boundary.** Training a re-ranker requires its own labeled-candidate-set construction, architecture choice, training loop, and validation pass — a real ML engineering project, not a quick add-on. Left for dedicated follow-up.
- **3E — stacker. DONE** (`phase3/PHASE3A_STACKER.md`) — logistic regression on out-of-fold base scores, tested against the tuned single-lever pipeline (density-adaptive + potency weighting + orthologue-for-primary).

Selection on the validation split; the locked scaffold test stays shut throughout.

### Phase 4 — Reliability and conformal (weeks 8–9)

**STATUS: item 5 tested — DOES NOT beat L-score, ship L-score.** (`phase4/PHASE4_CALIBRATION.md`.) Real Brier/ECE computed, out-of-fold, density-stratified isotonic calibration: Brier 0.1604 vs. L-score's 0.1532 (worse), ECE 0.0360 vs. 0.0384 (marginally better) — fails the "beat on both" gate. **Disposition: ship L-score (l/10) as the confidence measure**, per the plan's own no-complexity-credit rule. Scope actually covered: density-only stratification (not the full density × target-size × drug/non-drug per item 2 below); scaffold-strict captures only (not a formal three-way split against `holdouts.json`); isotonic only (no Platt/beta for small strata, item 4). A follow-up with the full stratification could plausibly close the gap — not attempted in this pass, disclosed as a real next step.

1. Mondrian cross-conformal, not marginal-only [21]; no data rebalancing (Mondrian alone handles imbalance [22]). **Not implemented — the calibration test above used out-of-fold isotonic regression, not formal Mondrian cross-conformal p-values/prediction sets.** A real scope reduction for this pass.
2. Strata: **density** (replacing a coarser similarity/domain proxy, per §2's unification finding) × target-size × drug/non-drug. Coarse; merge until ≥100 calibration points per stratum. **Only density implemented and tested** (4 coarse bands, n=149-585 each). Target-size and drug/non-drug dimensions not layered in.
3. Coverage under all four split types, including temporal — under-coverage there is a deliverable, not a failure. **Document/temporal remain blocked** (no per-pair dates, disclosed since Phase 1/2); random/scaffold holdouts exist (`phase1/data/holdouts.json`) but weren't used for this pass's calibration test.
4. Calibration method by stratum size: isotonic ≥1000 points, Platt/beta below, test Venn-ABERS. **Only isotonic implemented** (used uniformly, with a global-fit fallback for small bands rather than Platt/beta). Venn-ABERS not tested.
5. **Must beat the L-score baseline** (already measured, §2) on Brier and ECE, or ship the simple L-score instead — no complexity credit without earning it. **Tested — fails. Ship L-score** (see status line above).
6. D5 abstention policy set from the fitted density-adaptive rule (§ Phase 1 item 11), not a fixed 0.4 floor. **Policy stated, not newly computed**: abstain (or flag low-confidence) below density=8, the re-fit threshold (`phase3/PHASE3A_DENSITY_REFIT.md`) — supersedes both the old density≥12 number and any fixed 0.4 similarity floor.

### Phase 5 — Differentiators, rescoped after rev5's E5 correction (weeks 9–10)

**STATUS: assessed, see `PHASE5_STATUS.md` for full detail.** D5 done. D2's stated de-flooding mechanism (3A's popularity correction) was tested and rejected — needs a different approach, recommendation given. D1 and the Novelty position remain blocked on human steps (adjudication reviewer pass; literature review) that cannot be fabricated. D3 remains an explicitly unresolved human decision (rev5 #9). D4/D7 have their underlying data ready; both need separate frontend UI work outside this pass's scope.

- **D1 — reframed.** The *concept* (reliability score) is [1]'s and MolTarPred's, not ours [20]. Our contribution: calibration, conformal sets, density-based stratification, and a stated incompleteness range from the adjudication study. Wording: "for queries like yours, about X% of predictions at this evidence level are already-supported targets; adjudication suggests a further Y% are plausible but untested."
- **D2 — polypharmacology grouping**, de-flooded by 3A's popularity correction. Top group named "high-evidence predicted targets," never "primary."
- **D3 — QSAR/docking consensus.** Gated experiment on the 64-target subset only, like-for-like (G4). In/out decision still needed (rev5 decision #9, unresolved).
- **D4 — evidence panel**: neighbour structures, potency, assay type, species, ChEMBL links.
- **D5 — abstention**, density-driven per Phase 4.
- **D7 — orthologue transparency.** Surface when a prediction rests on non-human data — [1]'s own error analysis shows this is exactly where apparent false positives concentrate.

**Novelty position (rev5 §9, still required before any external "first" claim):** formal novelty check, not a limited search — three of Rev 4's claims fell to one afternoon of reading the primary literature; assume the same base rate applies here.

### Phase 6 — Pathogen module (V2.1, post-release)

**STATUS: real-data scoping done, no production build** (`phase6/PHASE6_PATHOGEN_MODULE_STATUS.md`) — appropriate depth for an explicitly post-release, off-critical-path phase. Key numbers: 665 bacterial single-protein targets (6.0% of the 11,055-target cross-organism ChEMBL universe, found via NCBI taxonomy classification — one real bug fixed, NCBI's esummary has no `lineage` field, use `genbankdivision` instead); only 24,897 base activities across them (~37/target vs. human's ~427/target) — confirms bacterial queries will rarely reach the density≥8 pooling threshold, validating the plan's own mechanism-class-fallback requirement. Human-homology proxy (reusing the orthologue tier's gene-symbol method): 1.8% overlap (12/665), favorable for target selection. Essentiality: not attempted, needs genuinely external data (DEG/OGEE) with no integration path built. Mechanism-class fallback: ChEMBL's own classification tree is sparse-to-absent for bacterial targets; `EC_NUMBER` synonyms identified as a viable alternative basis, not built. Pocket-embedding: correctly not attempted (plan's own "experimental at most"). Nothing here blocks or delays G1/Phase 7.

### Phase 7 — Release (weeks 11–13)

**STATUS: readiness assessed, NOT released.** See `PHASE7_RELEASE_READINESS.md` for the full honest checklist. **The ChEMBL date-retention pull was launched and completed** (54,717 documents, `phase2/stage_document_year_backfill.py`) — the temporal holdout and document-level split are now both built and real (two real bugs found and fixed along the way: a mega-document catastrophically skewing the document split, only caught by actually running it). **G1 now has exactly one remaining blocker**: the adjudication study's human reviewer pass. **The locked scaffold test was deliberately NOT opened** — doing so without a resolved authorization process (§9 item 4) would be exactly the kind of irreversible, undisclosed action this program has avoided throughout.

---

## 6. Gates

| Gate | Status / criterion |
|---|---|
| G0b | **PASS.** Closed — see §2. |
| G1 | **Six of seven conditions done; one blocker remains.** Benchmark card ✅; locked scaffold test sealed ✅ (still sealed, not opened); temporal holdout ✅ (built, `phase2/build_temporal_holdout.py`); four split types ✅ (all built, `phase2/build_document_split.py` completed the set); Set C confirmed drug-realistic ✅; full-population coverage numbers ✅; **adjudication study reported as a precision range ❌ — the sole remaining blocker**, needs a human reviewer pass on `phase1/data/adjudication_pilot_output.csv` (`PHASE1_ADJUDICATION_STUDY_STATUS.md`). See `PHASE7_RELEASE_READINESS.md` for full detail. |
| G2 | **PASSED, on the clean subset.** (`PHASE3_G2_LOCKED_TEST_RESULT.md`.) Locked scaffold test opened 2026-09-23 per explicit user authorization, after a fresh-holdout sanity check showed real signal. **Critical finding caught before scoring**: 60% of the 795 locked compounds had already been used in Phase 3 tuning (a real methodological gap — tuning captures sampled from `eval_sets.json` directly, never excluded the locked test's scaffold groups). Reported result uses only the 370 genuinely clean compounds: v2 beats v1 `best_similarity` on any-tier (+0.034, CI [0.012,0.059]) and primary-tier (+0.072, CI [0.021,0.128]), and beats the internal unweighted-10NN baseline on any-tier (+0.031, CI [0.004,0.059]) — all excluding zero. Primary-tier vs. baseline not significant (likely underpowered, n=103). **The locked test is now spent — do not reuse for further tuning.** Temporal-holdout confirmation and per-density-stratum reporting not yet done (a real follow-up, not required to clear this pass's bar). |
| G3 | Calibration: Brier, ECE, marginal + stratum-conditional coverage, 95% CI, n per stratum, all four splits. **Must beat the L-score baseline** (§2) or ship L-score instead. |
| G4 | D3 in/out per like-for-like precision on the 64-target subset. |
| G5 | External comparison, report-only, not pass/fail (even [3]'s own 7-method comparison didn't match predicted-target counts across methods). |
| G6 | Honesty gate: calibrated-event definition, incompleteness range, abstention policy, orthologue provenance — all stated, no bare percentage without its definition. |
| G7 | Reproducibility: clean checkout rebuilds the benchmark and headline table from the frozen index hash + recorded seeds. |

---

## 7. Statistical protocol (rev5 §8, unchanged, restated as a checklist)

- Pre-register one primary metric before Phase 3 (proposal: macro-averaged precision at matched breadth, temporal holdout), secondaries, margins, analysis plan — timestamped into the repo.
- Splits: tuning / validation / locked scaffold test (by scaffold cluster) + renewable temporal holdout.
- Nested CV for every tuned hyperparameter whose gain gets reported (infrastructure already built and validated — `phase0b/analyze_h4*.py`).
- Paired BCa bootstrap, 10,000 resamples, clustered by Bemis–Murcko scaffold group, never per-compound (infrastructure already built and sanity-tested — `phase0b/bca.py`).
- Both averaging conventions reported once per table where it matters; macro is primary thereafter, always labelled.
- Holm–Bonferroni across secondaries; pre-registered primary uncorrected.
- Power calculation before fixing G2's margin, using the formula in rev5 §8.7 — `σ_diff` is now estimable from real Phase 0b paired differences, not assumed.
- Effect sizes with CIs in every table, every phase report. No bare point estimates — this has been the house standard since Phase 0b and should not regress in Phase 1+.
- Sample size printed beside every stratified statistic (the L-score table already does this; keep doing it).

---

## 8. Reusable infrastructure (do not rebuild)

| Path | What it does |
|---|---|
| `target_fishing_v1_freeze/` | Frozen, checksummed v1 snapshot — the fixed reference point, never touched |
| `target_prediction_v2/phase0/` | `fetch_chembl_reference.py`, `build_eval_sets.py`, `build_leakage_fp2.py`, `leakage.py`, `metrics.py`, `run_phase0.py` — eval-set construction, dual-fingerprint leakage control, core metrics |
| `target_prediction_v2/phase0b/capture.py`, `capture_scaffold_strict.py` | Rich per-query neighbour capture (near-dup-only and scaffold-strict leakage variants), schema reusable for Phase 1's full-population scale-up — just raise sample sizes |
| `target_prediction_v2/phase0b/score.py` | Generic re-scoring engine (any k, α, potency/confidence weight function) over captured neighbour data — no re-querying needed for new weighting experiments |
| `target_prediction_v2/phase0b/bca.py` | From-scratch, sanity-tested scaffold-clustered BCa bootstrap |
| `target_prediction_v2/phase0b/analyze_h4*.py` | Nested-CV + BCa harness, reusable for any future retrieval-variant comparison, not just H4 |
| `target_prediction_v2/phase0b/diagnose_A_vs_B.py`, `test_density_stratification.py` | Density/spread/depth diagnostic pattern — directly extensible to Phase 1's decile-resolution confirmatory re-test |
| `target_prediction_v2/phase0b/fetch_mechanism_support.py`, `scope_orthologues.py` | Per-document ChEMBL pull and orthologue-candidate scoping — patterns to scale to full-population in Phase 2 |

All of Phase 1's "entry requirements" (§5, items 8–11) are extensions of code that already exists and already works — this is a scale-up, not a new build.

---

## 9. Decisions still needing explicit sign-off (consolidated from rev5 §12, rev6 §4, rev7 §3)

1. Confidence-filter relaxation for Phase 2 (§5, Phase 2 item 2) — without it, the confidence weight is a dead feature.
2. Conformal coverage level (80% or 90%) and minimum stratum size (proposal: 100).
3. D3 in/out of scope as a gated experiment (rev5 decision #9).
4. **RESOLVED.** Locked scaffold test opened 2026-09-23 per explicit user authorization, after a fresh-holdout sanity check justified it. See `PHASE3_G2_LOCKED_TEST_RESULT.md`. The test is now spent — this decision does not reopen for a second look.
5. Novelty position (rev5 §9) signed off before any external "first" claim is drafted.
6. **RESOLVED.** G1's literal text (§6) requires "temporal holdout built" and "four split types defined," both genuinely blocked on per-pair date/document data that only Phase 2's rebuild can produce — so G1 cannot pass before Phase 2 runs. **Decision: G1 is evaluated once, at the Phase 2→Phase 3 boundary**, not at the Phase 1→Phase 2 boundary — this is forced by G2's own text (§6), which requires candidates be evaluated "primary on the temporal holdout," meaning the temporal holdout must already exist by the time Phase 3's bake-off is gated by G2, which only happens after Phase 2 runs. **Phase 2 starts now**, in parallel with Phase 1's one remaining item (adjudication `final_verdict`, pending human review) — Phase 2's engineering work (potency/confidence weighting, salt/parent collapse, orthologue tier, document/date retention) does not depend on the adjudication result or any other still-open Phase 1 item.

**Not open, do not re-ask:** averaging convention (macro primary — settled), ground-truth recipe's core shape (§4 — settled, modulo item 1 above), mechanistically-supported operational rule (rev5 §4.1's 2-of-3 — adopted and already used), H4b/weighting functional form (closed, §2).

---

## References

Full list: `target_prediction_v2_procedure_rev5.md` items 1–25, `_rev6.md` item 26. No new references introduced in this plan.
