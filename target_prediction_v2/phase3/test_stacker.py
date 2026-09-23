"""
Phase 3E: logistic-regression stacker over three base scores (native k=10
unweighted pooling; potency-weighted pooling at the adopted config
thr=6.0/floor=0.0/steep=5.0; orthologue additive score) -- pointwise
learning-to-rank, out-of-fold (5-fold scaffold-grouped CV, reusing
analyze_h4.group_kfold), tested against the current best single-lever
(potency-weighted alone) per BUILD_PLAN.md's own gate: GBM/anything fancier
only if logistic regression wins first.

Usage: python3 test_stacker.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

PHASE0B_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0b")
sys.path.insert(0, PHASE0B_DIR)
import score as S  # noqa: E402
import bca as BCA  # noqa: E402
import analyze_h4 as H4  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2_DIR = os.path.join(HERE, "..", "phase2")
RESULTS_DIR = os.path.join(HERE, "results")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

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


def main():
    records = load_records()
    print(f"{len(records)} pooled queries loaded", flush=True)

    sys.path.insert(0, HERE)
    from fit_density_adaptive_rule_v2 import load_scaffolds_v2
    real_scaffolds = load_scaffolds_v2(records)
    scaffold_groups = [f"{r['_set']}:{sc}" for r, sc in zip(records, real_scaffolds)]

    v2_df = pd.read_csv(os.path.join(PHASE2_DIR, "v2_index", "compounds.csv.gz"))
    v2_fps = np.load(os.path.join(PHASE2_DIR, "v2_index", "fingerprints.npz"))["fps"]
    first_idx = v2_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]
    v2_pop = _POPCOUNT_TABLE[v2_fps].sum(axis=1)
    ortho_df = pd.read_csv(os.path.join(HERE, "orthologue_index", "compounds.csv.gz"))
    ortho_fps = np.load(os.path.join(HERE, "orthologue_index", "fingerprints.npz"))["fps"]
    ortho_pop = _POPCOUNT_TABLE[ortho_fps].sum(axis=1)

    def pfn(pchembl):
        return S.potency_weight(pchembl, **POTENCY_CFG)

    print("Computing base scores per query...", flush=True)
    native_scores, potency_scores, ortho_scores = [], [], []
    for i, r in enumerate(records):
        native_scores.append(S.score_query(r["neighbours"], k=10, alpha=0))
        potency_scores.append(S.score_query(r["neighbours"], k=10, alpha=0, potency_fn=pfn))
        ortho_scores.append(ortho_score_for(r["smiles"], first_idx, v2_fps, v2_pop, ortho_df, ortho_fps, ortho_pop))
        if (i + 1) % 400 == 0:
            print(f"  {i+1}/{len(records)}", flush=True)

    print("Building per-(query,candidate) feature table...", flush=True)
    rows = []
    for i, r in enumerate(records):
        candidates = set(native_scores[i]) | set(potency_scores[i]) | set(ortho_scores[i])
        annotated = set(r["annotated_targets"])
        for t in candidates:
            rows.append({
                "query_idx": i, "target": t,
                "f_native": native_scores[i].get(t, 0.0),
                "f_potency": potency_scores[i].get(t, 0.0),
                "f_ortho": ortho_scores[i].get(t, 0.0),
                "label": 1 if t in annotated else 0,
            })
    table = pd.DataFrame(rows)
    print(f"{len(table)} (query,candidate) rows, {table['label'].sum()} positive", flush=True)

    folds = H4.group_kfold(len(records), scaffold_groups, k=5, seed=42)
    query_to_fold = {}
    for fold_i, idxs in enumerate(folds):
        for qi in idxs:
            query_to_fold[qi] = fold_i

    table["fold"] = table["query_idx"].map(query_to_fold)
    feat_cols = ["f_native", "f_potency", "f_ortho"]

    stacked_score = np.zeros(len(table))
    for fold_i in range(5):
        train = table[table["fold"] != fold_i]
        test_mask = (table["fold"] == fold_i).values
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(train[feat_cols], train["label"])
        stacked_score[test_mask] = clf.predict_proba(table.loc[test_mask, feat_cols])[:, 1]
    table["stacked"] = stacked_score

    print("Ranking per query, comparing stacked vs. potency-weighted-alone...", flush=True)
    diffs_any, diffs_prim, groups_prim = [], [], []
    for i, r in enumerate(records):
        sub = table[table["query_idx"] == i]
        stacked_ranked = sorted(zip(sub["target"], sub["stacked"]), key=lambda kv: (-kv[1], kv[0]))
        potency_ranked = S.rank_from_scores(potency_scores[i])

        stacked_rank_any = S.best_rank(stacked_ranked, r["annotated_targets"])
        potency_rank_any = S.best_rank(potency_ranked, r["annotated_targets"])
        diffs_any.append(rr(stacked_rank_any) - rr(potency_rank_any))

        if r["primary_targets"]:
            stacked_rank_prim = S.best_rank(stacked_ranked, r["primary_targets"])
            potency_rank_prim = S.best_rank(potency_ranked, r["primary_targets"])
            diffs_prim.append(rr(stacked_rank_prim) - rr(potency_rank_prim))
            groups_prim.append(scaffold_groups[i])

    ci_any = BCA.scaffold_clustered_bca(diffs_any, scaffold_groups)
    ci_prim = BCA.scaffold_clustered_bca(diffs_prim, groups_prim) if diffs_prim else None
    print(f"any-tier (stacked vs potency-alone): {ci_any}", flush=True)
    print(f"primary-tier: {ci_prim}", flush=True)

    out = {"ci_any": ci_any, "ci_prim": ci_prim, "n_rows": len(table), "n_positive": int(table["label"].sum())}
    with open(os.path.join(RESULTS_DIR, "stacker_test.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {os.path.join(RESULTS_DIR, 'stacker_test.json')}", flush=True)


if __name__ == "__main__":
    main()
