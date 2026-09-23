"""
Phase 1, BUILD_PLAN.md items 1 and 3 (partial): construct the locked
scaffold test set and the random/scaffold split types for the tuning/
validation population.

Scope, disclosed directly:
  - Locked scaffold test (20%, sealed): BUILT HERE, real, deterministic,
    seeded. This is the actual holdout partition -- opened once at G2,
    never touched before then.
  - Random and scaffold split types (for Phase 3 model tuning/validation
    within the remaining 80%): BUILT HERE.
  - Document/assay-campaign and temporal split types: NOT built here --
    genuinely blocked. v1's aggregated index (compounds.csv.gz) retains no
    document_chembl_id or activity date for the ~1.3M-pair FULL reference
    (only for the small, query-scoped mechanism-support pull already done
    in Phase 0b, which covers query drugs, not the reference pool this
    split needs). BUILD_PLAN.md Phase 2 already commits to retaining
    document/assay metadata at the next index rebuild -- these two split
    types are correctly deferred to after that rebuild, not silently
    skipped now.

Population: the FULL matched population (Sets A + B, ~4,715 drugs, not
the 800/set sample used for the density work) -- this is a data-partition
exercise, so it should cover everything Phase 1 forward will actually use,
not just the exploratory sample.

Partition unit: Bemis-Murcko scaffold group (never split a scaffold group
across partitions -- the same leakage-safety principle used throughout
this program).

Usage: python3 build_holdouts.py
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "target_fishing_v1_freeze", "code"))
os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))

import target_fishing as TF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 42
LOCKED_TEST_FRACTION = 0.20
VALIDATION_FRACTION_OF_REMAINDER = 0.25  # of the 80% left after locking the test


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def main():
    phase0_dir = os.path.join(HERE, "..", "phase0")
    with open(os.path.join(phase0_dir, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)

    drugs = []
    for key, setname in [("set_A_approved", "A"), ("set_B_clinical", "B")]:
        for r in eval_sets[key]:
            drugs.append({"smiles": r["smiles"], "molecule_chembl_id": r["molecule_chembl_id"], "set": setname})
    # de-dup by smiles (a handful of compounds can appear with max_phase
    # ambiguity between the two ChEMBL pulls) -- keep first occurrence
    seen = set()
    dedup = []
    for d in drugs:
        if d["smiles"] in seen:
            continue
        seen.add(d["smiles"])
        dedup.append(d)
    drugs = dedup
    print(f"Full matched population: {len(drugs)} distinct drugs (Sets A+B combined)", flush=True)

    _, _, full_df = TF._load()
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))

    groups = {}
    for i, d in enumerate(drugs):
        sk = scaffold_key(d["smiles"], scaffold_by_smiles.get(d["smiles"]))
        groups.setdefault(sk, []).append(i)
    group_keys = list(groups.keys())
    print(f"{len(group_keys)} distinct scaffold groups", flush=True)

    rng = random.Random(SEED)
    rng.shuffle(group_keys)

    n_total = len(drugs)
    target_locked_n = int(round(n_total * LOCKED_TEST_FRACTION))

    locked_groups, locked_idx = [], []
    remainder_groups = []
    running = 0
    for gk in group_keys:
        if running < target_locked_n:
            locked_groups.append(gk)
            locked_idx.extend(groups[gk])
            running += len(groups[gk])
        else:
            remainder_groups.append(gk)

    remainder_idx_by_group = {gk: groups[gk] for gk in remainder_groups}
    rng2 = random.Random(SEED + 1)
    remainder_keys_shuffled = list(remainder_groups)
    rng2.shuffle(remainder_keys_shuffled)
    n_remainder = sum(len(groups[gk]) for gk in remainder_groups)
    target_val_n = int(round(n_remainder * VALIDATION_FRACTION_OF_REMAINDER))

    validation_groups, validation_idx = [], []
    tuning_groups, tuning_idx = [], []
    running = 0
    for gk in remainder_keys_shuffled:
        if running < target_val_n:
            validation_groups.append(gk)
            validation_idx.extend(groups[gk])
            running += len(groups[gk])
        else:
            tuning_groups.append(gk)
            tuning_idx.extend(groups[gk])

    # ---- scaffold split type: tuning vs validation is ALREADY a scaffold
    # split by construction above. Also build a RANDOM split type (same
    # tuning/validation sizes, but by compound not scaffold group) as the
    # second split type, for comparison during Phase 3 model tuning.
    rng3 = random.Random(SEED + 2)
    remainder_compound_idx = list(remainder_idx_by_group.values())
    flat_remainder = [i for sub in remainder_compound_idx for i in sub]
    rng3.shuffle(flat_remainder)
    n_val_random = int(round(len(flat_remainder) * VALIDATION_FRACTION_OF_REMAINDER))
    random_split_validation_idx = flat_remainder[:n_val_random]
    random_split_tuning_idx = flat_remainder[n_val_random:]

    def summarize(idx_list, label):
        subset = [drugs[i] for i in idx_list]
        n_a = sum(1 for d in subset if d["set"] == "A")
        n_b = sum(1 for d in subset if d["set"] == "B")
        return {"label": label, "n": len(subset), "n_from_set_A": n_a, "n_from_set_B": n_b}

    out = {
        "seed": SEED,
        "n_total_drugs": n_total,
        "n_scaffold_groups": len(group_keys),
        "partition_unit": "Bemis-Murcko scaffold group (never split across partitions)",
        "locked_scaffold_test": {
            **summarize(locked_idx, "locked_scaffold_test (SEALED — open once, at G2)"),
            "n_scaffold_groups": len(locked_groups),
            "smiles": [drugs[i]["smiles"] for i in locked_idx],
        },
        "tuning_validation_remainder": {
            "n_scaffold_groups_total": len(remainder_groups),
            "scaffold_split_type": {
                "tuning": {**summarize(tuning_idx, "tuning (scaffold split)"), "n_scaffold_groups": len(tuning_groups)},
                "validation": {**summarize(validation_idx, "validation (scaffold split)"), "n_scaffold_groups": len(validation_groups)},
            },
            "random_split_type": {
                "tuning": summarize(random_split_tuning_idx, "tuning (random split)"),
                "validation": summarize(random_split_validation_idx, "validation (random split)"),
            },
        },
        "document_assay_campaign_split_type": {
            "status": "BLOCKED",
            "reason": ("v1's aggregated index retains no document_chembl_id for the ~1.3M-pair full reference "
                       "(only for Phase 0b's small query-scoped mechanism-support pull, which covers query drugs "
                       "only, not the reference pool this split needs). Deferred to after Phase 2's index rebuild, "
                       "which already commits to retaining document/assay metadata per pair (BUILD_PLAN.md Section 4)."),
        },
        "temporal_split_type": {
            "status": "BLOCKED",
            "reason": ("No first-seen / ChEMBL-release date is currently retained per activity record. Would need "
                       "a fresh pull recording deposition/release dates -- a real, scoped-but-not-yet-done data "
                       "engineering task, not attempted here rather than faked."),
        },
    }

    out_dir = os.path.join(HERE, "data")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "holdouts.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps({k: v for k, v in out.items() if k != "locked_scaffold_test"}, indent=2, default=str), flush=True)
    print(f"\nlocked_scaffold_test: n={out['locked_scaffold_test']['n']} "
          f"(A={out['locked_scaffold_test']['n_from_set_A']}, B={out['locked_scaffold_test']['n_from_set_B']}), "
          f"{out['locked_scaffold_test']['n_scaffold_groups']} scaffold groups", flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
