"""
Phase 3A: orthologue-tier ablation -- does adding orthologue-derived
evidence (species_provenance != Homo sapiens, item 4's separate stratum)
to the k=10 pooled vote improve retrieval over native-human-only pooling?
Tested, not assumed, on the same n=1600 v2 captures.

Self-contained direct numpy/pandas search against orthologue_index
(reading files directly rather than through target_fishing.py's _load(),
which lru_caches a single module-level index directory and can't hold two
different indices open in one process).

Combination rule: combined_score(t) = native_score(t) + orthologue_score(t),
both on the same k=10-unweighted-vote scale (1.0 per qualifying neighbour
compound per target) -- a simple additive union, the natural first thing
to test before any more elaborate down-weighting of the orthologue term.

Usage: python3 test_orthologue_ablation.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd

PHASE0B_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0b")
sys.path.insert(0, PHASE0B_DIR)
import score as S  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2_DIR = os.path.join(HERE, "..", "phase2")
RESULTS_DIR = os.path.join(HERE, "results")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

NEIGHBOR_CAP = 300
ORTHO_TANIMOTO_MIN = 0.4


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_records():
    records = []
    for letter in ("A", "B"):
        with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_v2_{letter}.jsonl")) as f:
            for line in f:
                r = json.loads(line)
                r["_set"] = letter
                records.append(r)
    return records


def tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop):
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)


def main():
    records = load_records()
    print(f"{len(records)} pooled queries loaded", flush=True)

    # query fingerprint lookup: v2_index's own compounds+fingerprints,
    # first occurrence per smiles (same convention leakage.py uses)
    v2_df = pd.read_csv(os.path.join(PHASE2_DIR, "v2_index", "compounds.csv.gz"))
    v2_fps = np.load(os.path.join(PHASE2_DIR, "v2_index", "fingerprints.npz"))["fps"]
    first_idx = v2_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]
    v2_pop = _POPCOUNT_TABLE[v2_fps].sum(axis=1)

    ortho_df = pd.read_csv(os.path.join(HERE, "orthologue_index", "compounds.csv.gz"))
    ortho_fps = np.load(os.path.join(HERE, "orthologue_index", "fingerprints.npz"))["fps"]
    ortho_pop = _POPCOUNT_TABLE[ortho_fps].sum(axis=1)
    print(f"orthologue index: {len(ortho_df)} pairs, {ortho_df['smiles'].nunique()} distinct compounds", flush=True)

    sys.path.insert(0, HERE)
    from fit_density_adaptive_rule_v2 import load_scaffolds_v2
    real_scaffolds = load_scaffolds_v2(records)
    scaffold_groups = [f"{r['_set']}:{sc}" for r, sc in zip(records, real_scaffolds)]
    densities = [density_of(r) for r in records]

    diffs_any_all, diffs_any_dense, groups_dense = [], [], []
    diffs_prim_all, groups_prim_all = [], []
    n_with_ortho_evidence = 0

    for i, r in enumerate(records):
        native_score = S.score_query(r["neighbours"], k=10, alpha=0)

        smi = r["smiles"]
        if smi not in first_idx.index:
            combined_score = native_score
        else:
            row = int(first_idx.loc[smi]) if not hasattr(first_idx.loc[smi], "__len__") else int(first_idx.loc[smi].iloc[0])
            q_packed = v2_fps[row]
            q_pop = int(v2_pop[row])

            tanimoto = tanimoto_to_all(q_packed, ortho_fps, ortho_pop, q_pop)
            tmp = ortho_df.copy()
            tmp["tanimoto"] = tanimoto
            tmp = tmp[tmp["tanimoto"] >= ORTHO_TANIMOTO_MIN]
            by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
            top_n = by_compound.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)

            ortho_score = {}
            if not top_n.empty:
                # unweighted vote: each of the top-N orthologue neighbour
                # COMPOUNDS votes 1.0 for every human target it's linked to
                for nsmi in top_n.index:
                    targets_for_compound = set(tmp[tmp["smiles"] == nsmi]["target_chembl"])
                    for t in targets_for_compound:
                        ortho_score[t] = ortho_score.get(t, 0.0) + 1.0
                if ortho_score:
                    n_with_ortho_evidence += 1

            combined_score = dict(native_score)
            for t, s in ortho_score.items():
                combined_score[t] = combined_score.get(t, 0.0) + s

        native_ranked = S.rank_from_scores(native_score)
        combined_ranked = S.rank_from_scores(combined_score)

        native_rank_any = S.best_rank(native_ranked, r["annotated_targets"])
        combined_rank_any = S.best_rank(combined_ranked, r["annotated_targets"])
        diffs_any_all.append(rr(combined_rank_any) - rr(native_rank_any))

        if r["primary_targets"]:
            native_rank_prim = S.best_rank(native_ranked, r["primary_targets"])
            combined_rank_prim = S.best_rank(combined_ranked, r["primary_targets"])
            diffs_prim_all.append(rr(combined_rank_prim) - rr(native_rank_prim))
            groups_prim_all.append(scaffold_groups[i])

        if densities[i] >= 8:
            diffs_any_dense.append(rr(combined_rank_any) - rr(native_rank_any))
            groups_dense.append(scaffold_groups[i])

        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(records)} processed", flush=True)

    print(f"\n{n_with_ortho_evidence}/{len(records)} queries had >=1 orthologue neighbour "
          f"(Tanimoto>={ORTHO_TANIMOTO_MIN})", flush=True)

    ci_any_all = BCA.scaffold_clustered_bca(diffs_any_all, scaffold_groups)
    ci_any_dense = BCA.scaffold_clustered_bca(diffs_any_dense, groups_dense) if diffs_any_dense else None
    ci_prim_all = BCA.scaffold_clustered_bca(diffs_prim_all, groups_prim_all) if diffs_prim_all else None

    print(f"any-tier all: {ci_any_all}", flush=True)
    print(f"any-tier dense: {ci_any_dense}", flush=True)
    print(f"primary-tier all: {ci_prim_all}", flush=True)

    out = {"n_with_ortho_evidence": n_with_ortho_evidence, "n_total": len(records),
           "ci_any_all": ci_any_all, "ci_any_dense": ci_any_dense, "ci_prim_all": ci_prim_all}
    with open(os.path.join(RESULTS_DIR, "orthologue_ablation_test.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {os.path.join(RESULTS_DIR, 'orthologue_ablation_test.json')}", flush=True)


if __name__ == "__main__":
    main()
