"""
Phase 1, BUILD_PLAN.md item 6: MECHANICAL prep for the P1-A adjudication
study (rev5 Section 4.4) -- NOT the study itself.

What this script does (scriptable, no judgment involved):
  1. Score every candidate target for each captured query drug using
     weighted k-NN (k=25, alpha=1.0 -- the H4-confirmed baseline to beat,
     see PHASE1_BASELINES.md Section 2), from the already-captured
     per-query neighbour lists (phase0b/results/capture_scaffold_strict_{A,B}.jsonl).
  2. Rank each drug's candidate targets by score.
  3. Select the drug's top-ranked predictions that are NOT in its
     any-annotated ground-truth set -- i.e. apparent false positives under
     the loosest (most generous) ground-truth tier, which is exactly the
     population rev5 SS4.4 wants adjudicated (if a prediction is wrong even
     under the loosest tier, it's the most informative case: either a real
     annotation gap, or a genuine miss).
  4. Write a candidate worksheet (CSV) with blank columns for the two
     coding categories rev5 SS4.4 specifies, plus a free-text/verdict
     column -- ready for a human domain reviewer to fill in.

What this script explicitly does NOT do: make any adjudication judgment.
Whether a given (drug, target) false positive is actually an orthologue/
homologue annotation gap, a salt-form annotation gap, or a genuine miss
requires literature/database lookup per pair -- that is real domain-expert
work (rev5's own estimate: 2-3 days reviewer time), not something this
script or any script can produce. See PHASE1_ADJUDICATION_STUDY_STATUS.md.

Usage: python3 build_adjudication_candidates.py
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE0B_RESULTS = os.path.join(HERE, "..", "phase0b", "results")
DATA_DIR = os.path.join(HERE, "data")

K = 25
ALPHA = 1.0
TARGET_N_DRUGS = 54
MAX_FP_PER_DRUG = 4
TARGET_N_PAIRS = 180


def load_pooled():
    rows = []
    for letter in ("A", "B"):
        path = os.path.join(PHASE0B_RESULTS, f"capture_scaffold_strict_{letter}.jsonl")
        with open(path) as f:
            for line in f:
                d = json.loads(line)
                d["_set"] = letter
                rows.append(d)
    return rows


def weighted_knn_scores(neighbours, k=K, alpha=ALPHA):
    top_k = neighbours[:k]
    scores = {}
    for _smiles, sim, target_pairs in top_k:
        for tcid, _pchembl in target_pairs:
            scores[tcid] = scores.get(tcid, 0.0) + sim ** alpha
    return scores


def main():
    rows = load_pooled()
    print(f"{len(rows)} captured query drugs (pooled A+B)", flush=True)

    candidates = []
    for d in rows:
        scores = weighted_knn_scores(d["neighbours"])
        if not scores:
            continue
        annotated = set(d["annotated_targets"])
        ranked = sorted(scores.items(), key=lambda kv: -kv[1])
        fps = [(rank + 1, tcid, score) for rank, (tcid, score) in enumerate(ranked) if tcid not in annotated]
        if not fps:
            continue
        for rank, tcid, score in fps[:MAX_FP_PER_DRUG]:
            max_nb_sim = max((sim for _s, sim, tp in d["neighbours"][:K] if any(t == tcid for t, _p in tp)), default=None)
            candidates.append({
                "molecule_chembl_id": d["molecule_chembl_id"],
                "smiles": d["smiles"],
                "set": d["_set"],
                "max_phase": d["max_phase"],
                "predicted_target_chembl_id": tcid,
                "predicted_rank": rank,
                "weighted_knn_score": round(score, 4),
                "max_contributing_neighbour_similarity": max_nb_sim,
                "orthologue_or_homologue_annotation_gap": "",
                "salt_form_annotation_gap": "",
                "other_explanation": "",
                "reviewer_notes": "",
                "final_verdict": "",
            })

    # spread across as many distinct drugs as possible, favoring drugs
    # with a high-scoring (more confident, more interesting) top false
    # positive, up to TARGET_N_DRUGS / TARGET_N_PAIRS
    by_drug = {}
    for c in candidates:
        by_drug.setdefault(c["molecule_chembl_id"], []).append(c)
    drug_order = sorted(by_drug.keys(), key=lambda m: -max(c["weighted_knn_score"] for c in by_drug[m]))

    selected = []
    for mcid in drug_order[:TARGET_N_DRUGS]:
        drug_cands = sorted(by_drug[mcid], key=lambda c: c["predicted_rank"])
        selected.extend(drug_cands)

    n_drugs = len(set(c["molecule_chembl_id"] for c in selected))
    print(f"Selected {len(selected)} candidate pairs from {n_drugs} distinct drugs "
          f"(target: {TARGET_N_PAIRS} pairs / {TARGET_N_DRUGS} drugs, rev5 SS4.4)", flush=True)

    os.makedirs(DATA_DIR, exist_ok=True)
    out_path = os.path.join(DATA_DIR, "adjudication_candidates.csv")
    fieldnames = list(selected[0].keys())
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(selected)
    print(f"Wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
