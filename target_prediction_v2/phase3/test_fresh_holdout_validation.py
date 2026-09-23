"""
Phase A step 1 (making v2 usable): the real test. Runs the FULL ASSEMBLED
v2 pipeline -- density-adaptive gate (pool if density>=8 else best-
similarity) + potency weighting (thr=6.0/floor=0.0/steep=5.0) + additive
orthologue term -- against 800 queries that were NEVER used in any prior
Phase 3 tuning round (capture_fresh_holdout_{A,B}.jsonl), and compares it
against:
  (a) v1's actual production ranking (best_similarity)
  (b) the internal baseline every Phase 3A lever was tested against
      (native unweighted k=10 pooling)

This is the first time the ASSEMBLED pipeline (not individual levers) has
been tested on data it wasn't tuned on. Not the locked scaffold test
(that stays sealed) -- a cheap, reversible check of whether the adopted
improvements generalize before deciding whether opening the real locked
test is warranted.

Usage: python3 test_fresh_holdout_validation.py
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

DENSITY_THRESHOLD = 8  # phase3/PHASE3A_DENSITY_REFIT.md
NEIGHBOR_CAP = 300
ORTHO_TANIMOTO_MIN = 0.4
POTENCY_CFG = dict(threshold=6.0, floor=0.0, steepness=5.0)


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_records():
    records = []
    for letter in ("A", "B"):
        with open(os.path.join(RESULTS_DIR, f"capture_fresh_holdout_{letter}.jsonl")) as f:
            for line in f:
                r = json.loads(line)
                r["_set"] = letter
                records.append(r)
    return records


def tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop):
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)


def ortho_score_for(smi, first_idx, v2_fps, v2_pop, ortho_df, ortho_fps, ortho_pop):
    if smi not in first_idx.index:
        return {}
    row = int(first_idx.loc[smi]) if not hasattr(first_idx.loc[smi], "__len__") else int(first_idx.loc[smi].iloc[0])
    q_packed, q_pop = v2_fps[row], int(v2_pop[row])
    tanimoto = tanimoto_to_all(q_packed, ortho_fps, ortho_pop, q_pop)
    tmp = ortho_df.copy()
    tmp["tanimoto"] = tanimoto
    tmp = tmp[tmp["tanimoto"] >= ORTHO_TANIMOTO_MIN]
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    top_n = by_compound.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)
    out = {}
    for nsmi in top_n.index:
        for t in set(tmp[tmp["smiles"] == nsmi]["target_chembl"]):
            out[t] = out.get(t, 0.0) + 1.0
    return out


def bestsim_scores(record):
    return {t: s for t, s in record["bestsim_results"]}


def main():
    records = load_records()
    print(f"{len(records)} fresh (never-tuned-on) queries loaded", flush=True)

    sys.path.insert(0, HERE)
    from fit_density_adaptive_rule_v2 import load_scaffolds_v2
    real_scaffolds = load_scaffolds_v2(records)
    scaffold_groups = [f"{r['_set']}:{sc}" for r, sc in zip(records, real_scaffolds)]
    densities = [density_of(r) for r in records]

    v2_df = pd.read_csv(os.path.join(PHASE2_DIR, "v2_index", "compounds.csv.gz"))
    v2_fps = np.load(os.path.join(PHASE2_DIR, "v2_index", "fingerprints.npz"))["fps"]
    first_idx = v2_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]
    v2_pop = _POPCOUNT_TABLE[v2_fps].sum(axis=1)
    ortho_df = pd.read_csv(os.path.join(HERE, "orthologue_index", "compounds.csv.gz"))
    ortho_fps = np.load(os.path.join(HERE, "orthologue_index", "fingerprints.npz"))["fps"]
    ortho_pop = _POPCOUNT_TABLE[ortho_fps].sum(axis=1)

    def pfn(pchembl):
        return S.potency_weight(pchembl, **POTENCY_CFG)

    print("Scoring: best_similarity (v1), native k=10 pooling (internal baseline), "
          "and the FULL assembled v2 pipeline...", flush=True)
    diffs_v2_vs_bestsim_any, diffs_v2_vs_native_any = [], []
    diffs_v2_vs_bestsim_prim, diffs_v2_vs_native_prim, groups_prim = [], [], []
    n_pooled, n_bestsim_regime = 0, 0

    for i, r in enumerate(records):
        bs_scores = bestsim_scores(r)
        native_scores = S.score_query(r["neighbours"], k=10, alpha=0)

        if densities[i] >= DENSITY_THRESHOLD:
            n_pooled += 1
            v2_scores = S.score_query(r["neighbours"], k=10, alpha=0, potency_fn=pfn)
            ortho = ortho_score_for(r["smiles"], first_idx, v2_fps, v2_pop, ortho_df, ortho_fps, ortho_pop)
            for t, s in ortho.items():
                v2_scores[t] = v2_scores.get(t, 0.0) + s
        else:
            n_bestsim_regime += 1
            v2_scores = dict(bs_scores)  # density-adaptive gate: fall back to best-similarity

        bs_ranked = sorted(bs_scores.items(), key=lambda kv: (-kv[1], kv[0]))
        native_ranked = S.rank_from_scores(native_scores)
        v2_ranked = S.rank_from_scores(v2_scores)

        bs_rank_any = S.best_rank(bs_ranked, r["annotated_targets"])
        native_rank_any = S.best_rank(native_ranked, r["annotated_targets"])
        v2_rank_any = S.best_rank(v2_ranked, r["annotated_targets"])
        diffs_v2_vs_bestsim_any.append(rr(v2_rank_any) - rr(bs_rank_any))
        diffs_v2_vs_native_any.append(rr(v2_rank_any) - rr(native_rank_any))

        if r["primary_targets"]:
            bs_rank_prim = S.best_rank(bs_ranked, r["primary_targets"])
            native_rank_prim = S.best_rank(native_ranked, r["primary_targets"])
            v2_rank_prim = S.best_rank(v2_ranked, r["primary_targets"])
            diffs_v2_vs_bestsim_prim.append(rr(v2_rank_prim) - rr(bs_rank_prim))
            diffs_v2_vs_native_prim.append(rr(v2_rank_prim) - rr(native_rank_prim))
            groups_prim.append(scaffold_groups[i])

    print(f"\n{n_pooled} queries in the pooling regime (density>={DENSITY_THRESHOLD}), "
          f"{n_bestsim_regime} fell back to best-similarity (density<{DENSITY_THRESHOLD})", flush=True)

    ci_bs_any = BCA.scaffold_clustered_bca(diffs_v2_vs_bestsim_any, scaffold_groups)
    ci_native_any = BCA.scaffold_clustered_bca(diffs_v2_vs_native_any, scaffold_groups)
    ci_bs_prim = BCA.scaffold_clustered_bca(diffs_v2_vs_bestsim_prim, groups_prim)
    ci_native_prim = BCA.scaffold_clustered_bca(diffs_v2_vs_native_prim, groups_prim)

    print(f"\nFULL v2 pipeline vs. v1 best_similarity, any-tier:  {ci_bs_any}", flush=True)
    print(f"FULL v2 pipeline vs. native k=10 baseline, any-tier: {ci_native_any}", flush=True)
    print(f"FULL v2 pipeline vs. v1 best_similarity, primary:    {ci_bs_prim}", flush=True)
    print(f"FULL v2 pipeline vs. native k=10 baseline, primary:  {ci_native_prim}", flush=True)

    out = {
        "n_queries": len(records), "n_pooled_regime": n_pooled, "n_bestsim_regime": n_bestsim_regime,
        "v2_vs_bestsim_any": ci_bs_any, "v2_vs_native_any": ci_native_any,
        "v2_vs_bestsim_primary": ci_bs_prim, "v2_vs_native_primary": ci_native_prim,
    }
    with open(os.path.join(RESULTS_DIR, "fresh_holdout_validation.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {os.path.join(RESULTS_DIR, 'fresh_holdout_validation.json')}", flush=True)


if __name__ == "__main__":
    main()
