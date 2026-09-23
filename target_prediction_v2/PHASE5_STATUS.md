# Phase 5 — Differentiators: Status

Most of Phase 5's items are product/UI decisions or require a human step already flagged elsewhere in this program — not data-science tasks with a computable answer. Honest status per item, not fabricated completions.

## D1 — reframed reliability-score wording. BLOCKED on the same human step as before.
Wording template already fixed by the plan: *"for queries like yours, about X% of predictions at this evidence level are already-supported targets; adjudication suggests a further Y% are plausible but untested."* **X is now computable** — from `phase4/PHASE4_CALIBRATION.md`'s per-density-band calibration (or, since that failed its gate, from the L-score table directly, `PHASE0B_ADDENDUM.md` §7). **Y is not** — it requires the P1-A adjudication study's human reviewer pass, still pending (`PHASE1_ADJUDICATION_STUDY_STATUS.md`). D1's copy cannot ship complete until that reviewer pass lands.

## D2 — polypharmacology grouping. Dependency broken, needs a new decision.
Originally specified as "de-flooded by 3A's popularity correction" — **popularity correction was tested and rejected** (`phase3/PHASE3A_POPULARITY_CORRECTION.md`, significantly harmful at every tested strength). D2 cannot rely on a mechanism that doesn't exist. **Recommendation** (a genuine design call, not fabricated data): ship D2 without algorithmic de-flooding — rely on the wording safeguard already specified ("high-evidence predicted targets," never "primary") plus a simple display cap (e.g., top-10 shown, matching the ranking depth already used throughout this program's own metrics) rather than any popularity-based score adjustment. This is a scope reduction, not a full solution — a real de-flooding mechanism (if one is later found to work) would need its own tested lever, not a revival of the rejected one.

## D3 — QSAR/docking consensus. Genuinely unresolved — a human decision, not fabricated.
`BUILD_PLAN.md` itself marks this "rev5 decision #9, unresolved." Nothing in this pass changes that — in/out scope for this gated experiment is not something I can decide on the program's behalf. Left open, flagged, not silently dropped or invented.

## D4 — evidence panel. Data ready; UI implementation is separate frontend work.
Every field the panel needs already exists in the pipeline's own data: neighbour structures + Tanimoto (capture's `neighbours` field), potency (`pchembl_value`/`pchembl_value_computed`), assay type (`assay_types_seen`), species (`species_provenance`, orthologue tier), ChEMBL links (derivable from `target_chembl_id`/`molecule_chembl_id` directly). **Not attempted here**: the actual frontend panel component — a UI engineering task outside this session's data/research scope, not a research question with a data answer.

## D5 — abstention, density-driven. DONE.
Already covered in Phase 4's status: abstain (or flag low-confidence) below density=8, the just-refit threshold — supersedes the old density≥12 number and any fixed 0.4 similarity floor.

## D7 — orthologue transparency. Data ready; UI implementation is separate frontend work.
`species_provenance` is already a first-class field on every orthologue-tier pair (`stage7_orthologue_pairs.jsonl`), never merged into native-human rows (item 4's decision, enforced throughout Phase 3A's testing). Surfacing it in the UI when a prediction's evidence includes orthologue-tier data is a frontend task, not a data gap — the data to surface already exists and is already kept separate.

## Novelty position (rev5 §9). NOT attempted — a literature-review task, not fabricated.
"Formal novelty check, not a limited search" requires actually reading primary literature to check specific claims against prior art — a research task with no shortcut through data analysis. Not attempted in this pass; required before any external "first" claim is made, per the plan's own text.
