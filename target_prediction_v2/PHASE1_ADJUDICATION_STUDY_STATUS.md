# P1-A Adjudication Study — Status: Prepared, Not Run

**`BUILD_PLAN.md` Phase 1 item 6 (rev5 §4.4).** Rev 5/Rev 7 both call this "the best value-per-hour item in the whole plan." It is also the one Phase 1 item that **cannot be completed by this program autonomously** — it requires genuine human domain-expert judgment, per compound-target pair, against primary literature and database records. This document states that plainly rather than fabricating adjudication outcomes or silently skipping the item.

---

## What the study is

For a sample of the model's top-ranked "false positives" (predictions the current any-annotated ground truth doesn't support), determine — by actually looking up each pair in ChEMBL/UniProt/DrugBank/PubMed — whether the prediction is:

1. A genuine miss (the model found a spurious structural coincidence), or
2. A **real target relationship the ground truth is simply missing**, for one of two specifically named reasons (rev5 §4.4):
   - **Orthologue/homologue annotation gap**: the drug is annotated against target X in one species but the index's human-target scope missed a genuinely-supported human orthologue/paralogue.
   - **Salt-form annotation gap**: the drug's specific salt/hydrate form is registered as a separate ChEMBL molecule ID from the parent, splitting its real annotations across IDs so the query molecule's own record looks sparser than the compound's true evidence base.

This produces a **precision range** ("X% of top false positives are actually real, just unannotated"), not a point estimate — exactly the number `BUILD_PLAN.md` §9's D1 wording is written to use ("adjudication suggests a further Y% are plausible but untested").

## Why this program cannot do the adjudication itself

Determining which of the three categories a given (drug, target) pair falls into requires:
- Reading the actual primary literature or database entries for that specific drug-target pair (not a pattern this program's existing data captures).
- Domain judgment about what counts as a "genuinely supported" orthologue relationship vs. a coincidental structural similarity.
- Recognizing salt/hydrate parent-child relationships that ChEMBL's own molecule IDs don't always make mechanically obvious.

None of this is derivable from the index, the capture files, or any script — it is real reviewer time. Rev 5's own estimate: **2-3 days of a qualified domain reviewer**. Fabricating plausible-sounding verdicts here would be worse than leaving the item open — it would silently corrupt the exact number (`D1`'s incompleteness range) this program plans to make a user-facing claim about.

## What this pass DID complete: the mechanical prep

`phase1/build_adjudication_candidates.py` — fully scriptable, no judgment involved:

1. Scored every candidate target for all 1,600 pooled-sample (Set A+B) captured query drugs using **weighted k-NN** (k=25, α=1.0 — the confirmed baseline-to-beat, `PHASE1_BASELINES.md` §2), from the already-captured per-query neighbour lists.
2. Ranked each drug's candidate targets by score.
3. Selected each drug's top-ranked predictions that fall **outside** its any-annotated ground truth — the most informative population for adjudication, since a miss even under the loosest ground-truth tier is either a real annotation gap or a genuine error, not a threshold artifact.
4. Wrote a ready-to-fill worksheet.

**Output: `phase1/data/adjudication_candidates.csv` — 153 candidate pairs across 54 distinct drugs** (rev5 §4.4 target: 150-200 pairs, 40-50 drugs — within range on pairs, slightly over on drug count because per-drug false-positive yield was uneven; widening to 54 drugs was needed to reach 150 pairs at a cap of 4 candidates/drug). Columns: `molecule_chembl_id`, `smiles`, `set`, `max_phase`, `predicted_target_chembl_id`, `predicted_rank`, `weighted_knn_score`, `max_contributing_neighbour_similarity`, plus five blank columns for the reviewer: `orthologue_or_homologue_annotation_gap`, `salt_form_annotation_gap`, `other_explanation`, `reviewer_notes`, `final_verdict`.

## Update: AI-assisted evidence-gathering pass completed (still not the adjudication itself)

Rather than leave the worksheet fully blank, an evidence-gathering pass was run over all 153 rows (`phase1/adjudication_pilot.py`) to triage the reviewer's workload. Design, deliberately grounded rather than free-text literature guessing:

- **Salt-form gap check (mechanical, database fact)**: queries ChEMBL for sibling molecule IDs sharing the same `molecule_hierarchy.parent_chembl_id` as the query drug, then checks directly whether any sibling already has an annotated activity against the predicted target.
- **Orthologue/homologue gap check (mechanical, database fact)**: looks up the predicted target's ChEMBL gene-symbol synonym, finds other single-protein targets sharing that gene symbol (catches cross-species orthologues ChEMBL names identically), and checks directly whether the *same query drug* already has an annotated activity against that non-human target.
- **PubMed co-occurrence (weak, explicitly unverified)**: only when neither mechanical check fires — searches drug name + target gene symbol/name in PubMed title/abstract, reported as a hit count and PMIDs, never as a verdict.
- **DrugBank**: not used — confirmed via a live connectivity test that its public pages return `403 Forbidden` (bot-blocked) and no API key is available.

**Result, `phase1/data/adjudication_pilot_output.csv` (153/153 rows, all `ai_*`-prefixed columns, the five original human-judgment columns left untouched):**

| Category | n | % | Basis |
|---|---:|---:|---|
| Salt-form gap | 15 | 9.8% | Direct ChEMBL database fact |
| Orthologue/homologue gap | 4 | 2.6% | Direct ChEMBL database fact |
| Uncertain | 35 | 22.9% | PubMed hit found, content unread/unverified |
| No evidence found | 99 | 64.7% | None of the three checks found anything |

**19/153 (12.4%) are database-confirmed real-but-unannotated relationships** — an auditable floor for D1's incompleteness range, not an estimate. This is a genuine, disclosed limitation set, not a complete substitute for review: DrugBank is inaccessible; PubMed matching is title/abstract only and gene-symbol-based (won't catch papers using neither); "no evidence found" for the 99 means the three cheap checks found nothing, not that a human read those cases and confirmed they're genuine misses.

## What happens next (needs a human) — now triaged, not blind

A domain reviewer's workload is reduced from 153 blind lookups to:
1. **19 rows** — fast confirmation (evidence is already a database fact; just check the categorization reads right).
2. **35 rows** — the real work: read the specific PubMed hits listed (`ai_pubmed_example_pmids` in the output file) and judge whether they actually support the predicted relationship.
3. **99 rows** — spot-check a sample rather than exhaustively review; flag any known false negatives from the mechanical checks' blind spots (DrugBank-only relationships, non-PubMed-indexed sources, synonym mismatches).

Once the reviewer fills in `final_verdict` (on `adjudication_candidates.csv` or the AI-evidence-augmented `adjudication_pilot_output.csv`), computing the resulting precision range is mechanical again.

**This item stays open until that reviewer pass happens.** `BUILD_PLAN.md` §9 gate G1 lists "adjudication study reported as a precision range" as a literal gate condition — the AI-assisted pass narrows the work substantially but does not substitute for it, and cannot be marked closed by this program alone.
