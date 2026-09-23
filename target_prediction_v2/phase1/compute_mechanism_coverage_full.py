"""
Phase 1, BUILD_PLAN.md item 7 (closing step): apply Rev 5 Section 4.1's
"mechanistically supported" rule -- a (molecule, target) pair counts if it
satisfies AT LEAST 2 OF 3 criteria (primary/intended via drug_mechanism;
>=2 independent documents; binding-assay evidence present) -- to the full
matched population (~4,715 drugs), using phase1/data/mechanism_support_full.json
(the batched, recipe-filtered pull) + eval_sets.json's annotated/primary
targets. Same literal rule Phase 0b applied to the 710-drug sample
(PHASE0B_ADDENDUM.md Section 13); this is the full-population extension,
not a redefinition.

Reports both:
  - pair-level coverage (% of annotated pairs meeting the rule, of pairs
    with pulled data -- matching Section 13's own denominator convention)
  - compound-level coverage (% of drugs with >=1 supported target -- matching
    PHASE0B_REV7_RESPONSE.md's headline numbers: A 50.5%, B 63.5%, D 74.3%
    on the 710-drug sample)

Usage: python3 compute_mechanism_coverage_full.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    with open(os.path.join(HERE, "..", "phase0", "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)
    with open(os.path.join(HERE, "data", "mechanism_support_full.json")) as f:
        support = json.load(f)

    out = {}
    for key, setname in [("set_A_approved", "A"), ("set_B_clinical", "B")]:
        n_drugs = 0
        n_drugs_with_support = 0
        n_pairs_annotated = 0
        n_pairs_with_data = 0
        n_pairs_supported = 0

        for r in eval_sets[key]:
            mcid = r["molecule_chembl_id"]
            annotated = r["annotated_targets"]
            primary = set(r["primary_targets"])
            n_drugs += 1
            n_pairs_annotated += len(annotated)

            sup_for_mol = support.get(mcid, {})
            drug_has_support = False
            for tcid in annotated:
                entry = sup_for_mol.get(tcid)
                if entry is None:
                    continue
                n_pairs_with_data += 1
                criteria_met = sum([
                    tcid in primary,
                    entry["n_documents"] >= 2,
                    entry["has_binding_assay"],
                ])
                if criteria_met >= 2:
                    n_pairs_supported += 1
                    drug_has_support = True
            if drug_has_support:
                n_drugs_with_support += 1

        out[setname] = {
            "n_drugs": n_drugs,
            "compound_level_coverage_pct": round(100 * n_drugs_with_support / n_drugs, 1),
            "n_drugs_with_ge1_supported_target": n_drugs_with_support,
            "n_pairs_annotated": n_pairs_annotated,
            "n_pairs_with_pulled_data": n_pairs_with_data,
            "pct_pairs_with_pulled_data": round(100 * n_pairs_with_data / n_pairs_annotated, 1) if n_pairs_annotated else None,
            "n_pairs_supported": n_pairs_supported,
            "pair_level_coverage_pct_of_pairs_with_data": round(100 * n_pairs_supported / n_pairs_with_data, 1) if n_pairs_with_data else None,
        }

    out_path = os.path.join(HERE, "data", "mechanism_coverage_full.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
