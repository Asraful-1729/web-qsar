# Adjudication Study — Hand-off to the Reviewer

**What this is for**: the one remaining task before this project is a fully-closed research tool. Everything else — the method, its validation against a real held-out test, a comparison against a published baseline (SEA), reproducibility, and a live interface — is done (`METHODS_AND_VALIDATION.md`, `RESEARCH_HANDOFF_ROADMAP.md`). This is the one piece that genuinely needs your judgment, not more engineering.

## The question

When the model ranks a target highly for a drug, but that (drug, target) pair isn't in the ground-truth benchmark, is that a **real miss** (the model found a spurious structural coincidence) or a **real relationship the benchmark's ground truth is simply missing**?

Both happen in ChEMBL-derived benchmarks. The two most common reasons a real relationship goes "missing" from the ground truth:
- **Salt-form gap**: the drug's salt/hydrate form is registered as a separate ChEMBL ID from its parent, so the query molecule's own record looks sparser than the compound's true evidence.
- **Orthologue/homologue gap**: the drug is annotated against the target in another species, but the index's human-only scope missed a genuinely-supported human ortholog.

The answer matters because it turns "the model has false positives" into a defensible, disclosed number: "X–Y% of the model's apparent false positives are actually real, just unannotated" — a precision range, not a point estimate.

## What's already done for you

153 candidate (drug, target) pairs — the model's top-ranked predictions that fall outside the benchmark's ground truth — have already been through an automated evidence-gathering pass (not a verdict, just evidence triage):

| Category | n | Basis |
|---|---:|---|
| Salt-form gap found | 15 | Direct ChEMBL database fact (sibling molecule ID already has this activity) |
| Orthologue/homologue gap found | 4 | Direct ChEMBL database fact (same drug already active on the non-human ortholog) |
| Uncertain | 35 | A PubMed title/abstract hit exists, but nobody has read it to confirm it's actually about this relationship |
| No evidence found | 99 | None of the automated checks found anything |

**File**: `phase1/data/adjudication_pilot_output.csv` (153 rows). The 5 right-most `ai_*`-prefixed columns are that automated evidence (drug/target names, gene symbols, the actual PubMed IDs for the "uncertain" rows, a suggested verdict). The `final_verdict` column is blank — that's yours.

## What to actually do

You don't need to review all 153 rows the same way:

1. **19 rows** (`ai_salt_form_gap` or `ai_orthologue_gap` = "Yes") — fast. The evidence is already a database fact; just confirm the automated categorization reads correctly and copy it into `final_verdict`.
2. **35 rows** (`ai_suggested_verdict` starts with "uncertain") — the real work. Each has a PubMed ID (`ai_pubmed_example_pmids`) already looked up for you — read the abstract and judge whether it actually supports the predicted drug-target relationship, or is just a coincidental co-occurrence (e.g., both terms appear in an unrelated review article).
3. **99 rows** ("genuine_miss") — spot-check a sample rather than reviewing exhaustively. The automated checks only look at ChEMBL sibling records and PubMed title/abstract — they can't see DrugBank (blocked, no API key available) or full-text literature, so a real relationship could still hide here. Flag any you happen to know about.

For each row, fill in `final_verdict` with one of:
- `real_but_unannotated (salt-form gap)`
- `real_but_unannotated (orthologue/homologue gap)`
- `real_but_unannotated (other)`
- `genuine_miss`
- `uncertain` (you looked, genuinely can't resolve it)

Rev 5's own estimate for a full pass: 2–3 reviewer-days. The triage above should make it faster in practice, especially if you only spot-check the 99.

## Getting the final number

Once you've filled in `final_verdict` (for as many or as few rows as you get to — partial review is fine and clearly reported as such):

```
python3 target_prediction_v2/phase1/compute_adjudication_precision_range.py
```

This reads your verdicts directly from `adjudication_pilot_output.csv` and prints the precision range — no further engineering needed, and it correctly reports "not yet reviewed" for anything you skip rather than guessing.

## Where this result goes

Once you have a range, it slots into `BUILD_PLAN.md`'s G1 gate (the one remaining unmet condition — see `PHASE7_RELEASE_READINESS.md`) and should be added to `METHODS_AND_VALIDATION.md` §6 (currently: "an adjudication study... needs a qualified human reviewer's judgment. This is the single most important open item.").

## Background, if useful

- `PHASE1_ADJUDICATION_STUDY_STATUS.md` — the full account of why this can't be automated, exactly how the evidence-gathering pass works, and its disclosed limitations.
- `METHODS_AND_VALIDATION.md` — what the model actually is and what's already been validated (start here for the big picture).
- `phase1/adjudication_pilot.py` — the script that produced the `ai_*` evidence columns, if you want to see exactly how each check works or re-run it after a ChEMBL update.
