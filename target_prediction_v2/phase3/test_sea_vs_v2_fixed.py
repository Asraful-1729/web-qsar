"""
Corrected re-run of test_sea_vs_v2.py, using the fixed continuous
power-law null model (build_sea_null_distributions_v2.py /
sea_null_distributions_v2.json) instead of the original discrete-bin null
(sea_null_distributions.json), which was found to be broken for
mega-promiscuous targets (up to 195,809 ligands) far outside its largest
bin's ~1640-ligand representative size -- see build_sea_null_distributions_v2.py's
docstring for the full diagnosis. That bug produced absurd z-scores
(up to z=641) for huge real targets regardless of true relevance, causing
them to dominate SEA's ranked output and explaining the implausible
original result (SEA scoring catastrophically worse than v1 best_similarity).

Uses the SAME query subset/seed/methodology as test_sea_vs_v2.py (n=150,
seed=999, fresh-holdout population) so the two results are a clean
apples-to-apples before/after comparison of the null-model fix alone.

Usage: python3 test_sea_vs_v2_fixed.py
"""
import json
import os
import random
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
PHASE2_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase2")
PHASE0B_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0b")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, PHASE0B_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(os.path.join(PHASE2_DIR, "v2_index"))

import target_fishing as TF  # noqa: E402
import leakage as LK  # noqa: E402
import score as S  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

MIN_LIGANDS = 5
NEAR_DUP_CUTOFF = 0.95
N_PER_SET = 75
SEED = 999
DENSITY_THRESHOLD = 8
POTENCY_CFG = dict(threshold=6.0, floor=0.0, steepness=5.0)
ORTHO_TANIMOTO_MIN = 0.4
NEIGHBOR_CAP = 300

EULER_MASCHERONI = 0.5772156649
GUMBEL_STD_FACTOR = 1.2825498  # pi / sqrt(6)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_power_law():
    with open(os.path.join(RESULTS_DIR, "sea_null_distributions_v2.json")) as f:
        d = json.load(f)
    pl = d["power_law"]
    return pl["mean_a"], pl["mean_b"], pl["std_a"], pl["std_b"]


def gumbel_params_for_size(n_ligands, mean_a, mean_b, std_a, std_b):
    mean = mean_a * (n_ligands ** mean_b)
    std = std_a * (n_ligands ** std_b)
    scale = std / GUMBEL_STD_FACTOR
    loc = mean - EULER_MASCHERONI * scale
    return loc, scale


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def load_fresh_records():
    records = []
    for letter in ("A", "B"):
        with open(os.path.join(HERE, "results", f"capture_fresh_holdout_{letter}.jsonl")) as f:
            for line in f:
                r = json.loads(line)
                r["_set"] = letter
                records.append(r)
    return records


def main():
    print("Loading v2_index, leak index, power-law null params...", flush=True)
    full_fps, full_pop, full_df = TF._load()
    fp2 = np.load(os.path.join(PHASE2_DIR, "data", "v2_leakage_fp2_filtered.npz"), allow_pickle=True)
    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2["smiles"], fp2["fps"])
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))
    scaffold_key_by_smiles = {s: scaffold_key(s, sc) for s, sc in scaffold_by_smiles.items()}

    mean_a, mean_b, std_a, std_b = load_power_law()
    print(f"Power law: mean(N)={mean_a:.5f}*N^{mean_b:.4f}, std(N)={std_a:.5f}*N^{std_b:.4f}", flush=True)

    all_fresh = load_fresh_records()
    rng = random.Random(SEED)
    rng.shuffle(all_fresh)
    subset = all_fresh[:N_PER_SET * 2]
    print(f"{len(subset)} queries selected for SEA comparison (subset of fresh-holdout, seed={SEED})", flush=True)

    v2_df = pd.read_csv(os.path.join(PHASE2_DIR, "v2_index", "compounds.csv.gz"))
    v2_fps = np.load(os.path.join(PHASE2_DIR, "v2_index", "fingerprints.npz"))["fps"]
    first_idx_map = v2_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]
    v2_pop_all = _POPCOUNT_TABLE[v2_fps].sum(axis=1)
    ortho_df = pd.read_csv(os.path.join(HERE, "orthologue_index", "compounds.csv.gz"))
    ortho_fps = np.load(os.path.join(HERE, "orthologue_index", "fingerprints.npz"))["fps"]
    ortho_pop = _POPCOUNT_TABLE[ortho_fps].sum(axis=1)

    def pfn(pchembl):
        return S.potency_weight(pchembl, **POTENCY_CFG)

    def tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop):
        inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
        union = q_pop + pop_counts - inter
        return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)

    def ortho_score_for(smi):
        if smi not in first_idx_map.index:
            return {}
        row = int(first_idx_map.loc[smi]) if not hasattr(first_idx_map.loc[smi], "__len__") else int(first_idx_map.loc[smi].iloc[0])
        q_packed, q_pop = v2_fps[row], int(v2_pop_all[row])
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

    diffs_v2_vs_sea_any, diffs_v2_vs_sea_prim, groups_prim = [], [], []
    diffs_sea_vs_bestsim_any = []
    n_done = 0
    t0 = time.time()
    scaffold_groups = []

    for row in subset:
        if n_done >= N_PER_SET * 2:
            break
        smi = row["smiles"]

        reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, NEAR_DUP_CUTOFF)
        q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)
        q_scaffold = scaffold_key_by_smiles.get(smi, f"__no_ring__:{smi}")
        same_scaffold_mask = reduced_df["smiles"].map(
            lambda s: scaffold_key_by_smiles.get(s, f"__no_ring__:{s}")) == q_scaffold
        keep_mask = ~same_scaffold_mask.values
        sc_fps, sc_pop, sc_df = reduced_fps[keep_mask], reduced_pop[keep_mask], reduced_df.loc[keep_mask]

        tanimoto = LK.tanimoto_to_all(q_packed, sc_fps, sc_pop, q_pop)
        tmp = sc_df.copy()
        tmp["tanimoto"] = tanimoto

        per_target = tmp.groupby("target_chembl").agg(
            raw_score=("tanimoto", "sum"),
            n_ligands=("smiles", "nunique"),
        )
        per_target = per_target[per_target["n_ligands"] >= MIN_LIGANDS]

        sea_pvals = {}
        for tcid, rowdata in per_target.iterrows():
            loc, scale = gumbel_params_for_size(int(rowdata["n_ligands"]), mean_a, mean_b, std_a, std_b)
            pval = float(stats.gumbel_r.sf(rowdata["raw_score"], loc=loc, scale=scale))
            sea_pvals[tcid] = pval

        sea_ranked = sorted(sea_pvals.items(), key=lambda kv: kv[1])

        neighbours = []
        by_compound2 = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
        top_n = by_compound2.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)
        if not top_n.empty:
            neighbor_rows = tmp[tmp["smiles"].isin(set(top_n.index))][["smiles", "target_chembl", "pchembl_value"]]
            by_smi = {}
            for r in neighbor_rows.itertuples(index=False):
                pv = None if r.pchembl_value != r.pchembl_value else float(r.pchembl_value)
                by_smi.setdefault(r.smiles, []).append([r.target_chembl, pv])
            for nsmi, sim in top_n["tanimoto"].items():
                neighbours.append([nsmi, round(float(sim), 4), by_smi.get(nsmi, [])])

        density = sum(1 for _, sim, _ in neighbours if sim >= 0.5)
        bs_results = TF._search_against(sc_fps, sc_pop, sc_df, q_packed, q_pop, threshold=0.0, max_compounds_per_target=1)
        bs_scores = {r["target_chembl"]: r["best_similarity"] for r in bs_results}
        bs_ranked = sorted(bs_scores.items(), key=lambda kv: (-kv[1], kv[0]))

        if density >= DENSITY_THRESHOLD:
            v2_scores = S.score_query(neighbours, k=10, alpha=0, potency_fn=pfn)
            ortho = ortho_score_for(smi)
            for t, s in ortho.items():
                v2_scores[t] = v2_scores.get(t, 0.0) + s
        else:
            v2_scores = dict(bs_scores)
        v2_ranked = S.rank_from_scores(v2_scores)

        sea_rank_any = S.best_rank(sea_ranked, row["annotated_targets"])
        v2_rank_any = S.best_rank(v2_ranked, row["annotated_targets"])
        bs_rank_any = S.best_rank(bs_ranked, row["annotated_targets"])
        diffs_v2_vs_sea_any.append(rr(v2_rank_any) - rr(sea_rank_any))
        diffs_sea_vs_bestsim_any.append(rr(sea_rank_any) - rr(bs_rank_any))

        sc = f"{row['_set']}:{q_scaffold}"
        scaffold_groups.append(sc)

        if row["primary_targets"]:
            sea_rank_prim = S.best_rank(sea_ranked, row["primary_targets"])
            v2_rank_prim = S.best_rank(v2_ranked, row["primary_targets"])
            diffs_v2_vs_sea_prim.append(rr(v2_rank_prim) - rr(sea_rank_prim))
            groups_prim.append(sc)

        n_done += 1
        if n_done % 10 == 0:
            elapsed = time.time() - t0
            rate = n_done / elapsed
            remaining = (N_PER_SET * 2 - n_done) / rate if rate > 0 else 0
            print(f"  {n_done}/{N_PER_SET*2} ({elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining)", flush=True)

    print(f"\n{n_done} queries scored (SEA-fixed + v2 + best_similarity)", flush=True)
    ci_v2_vs_sea_any = BCA.scaffold_clustered_bca(diffs_v2_vs_sea_any, scaffold_groups)
    ci_v2_vs_sea_prim = BCA.scaffold_clustered_bca(diffs_v2_vs_sea_prim, groups_prim) if diffs_v2_vs_sea_prim else None
    ci_sea_vs_bestsim_any = BCA.scaffold_clustered_bca(diffs_sea_vs_bestsim_any, scaffold_groups)

    print(f"\nv2 (full pipeline) vs SEA-fixed, any-tier:   {ci_v2_vs_sea_any}", flush=True)
    print(f"v2 (full pipeline) vs SEA-fixed, primary-tier: {ci_v2_vs_sea_prim}", flush=True)
    print(f"SEA-fixed vs v1 best_similarity, any-tier:   {ci_sea_vs_bestsim_any}", flush=True)

    out = {"n": n_done, "v2_vs_sea_any": ci_v2_vs_sea_any, "v2_vs_sea_primary": ci_v2_vs_sea_prim,
           "sea_vs_bestsim_any": ci_sea_vs_bestsim_any, "null_model": "continuous power-law (v2 fix)"}
    with open(os.path.join(RESULTS_DIR, "sea_vs_v2_result_fixed.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {os.path.join(RESULTS_DIR, 'sea_vs_v2_result_fixed.json')}", flush=True)


if __name__ == "__main__":
    main()
