"""
Phase 0 (target_prediction_v2_procedure_rev3.md, Section 4) -- the actual
autopsy run. Loads the FROZEN v1 snapshot read-only, the drug-realistic
eval sets (build_eval_sets.py) and the leakage-control fingerprints
(build_leakage_fp2.py), then runs:

  1. Baseline measurement of v1-as-shipped (ranked by best_similarity, its
     frozen production default) on Sets A (approved), B (clinical), D
     (sparse-target), with leakage control at Tanimoto>=0.95 (primary) --
     any-annotated and primary-target Top-k/MRR/median-rank, plus
     precision/recall/MCC/targets-per-query and the precision-coverage
     curve via post-hoc threshold filtering (see module docstring in
     metrics.py for why one search per query, at threshold=0.0, is
     sufficient to reconstruct every threshold's predicted set).
  2. Leakage sensitivity: the same queries re-run with the looser 0.90
     cutoff, to see how much the leakage boundary itself moves the numbers
     (review point 4).
  3. Potency-filter experiment (H3): re-run a subsample with the reduced
     reference additionally restricted to pchembl_value >= 5.0 (10 uM) and
     >= 6.0 (1 uM), to test whether v1's lack of a stated potency floor is
     inflating or depressing recovery.
  4. Weighted k-NN vs best_similarity (H4): for the same queries, also rank
     targets by score(t) = sum over the k nearest DISTINCT reference
     compounds i (by Tanimoto, post-leakage-removal) whose annotation set
     contains t, of tanimoto_i^alpha -- swept over a small (k, alpha) grid
     -- compared against the frozen best_similarity ranking on identical
     queries.
  5. Popularity baseline: a query-independent sanity floor (rank ALL
     targets by how many compounds support them globally, ignore the query
     molecule entirely).

This is a SAMPLE-based first pass (see SAMPLE_* constants below), not a
full census of every matched drug -- each query costs ~3-4s (two full
857k-compound leakage-fingerprint passes plus a 1.3M-row search), so the
full ~4,700-drug matched set would take several hours. The sample sizes
below were chosen to finish in roughly an hour while still giving a
genuine, honestly-labeled statistical sample (Gate G3 explicitly calls for
"n per stratum with 95% CI," which presupposes sampling is legitimate
methodology here, not a shortcut). Scaling up is mechanical: raise the
SAMPLE_* constants and re-run; no code changes needed.

Usage: python3 run_phase0.py
"""
import json
import os
import random
import statistics
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "target_fishing_v1_freeze", "target_fishing_index"))

import target_fishing as TF  # noqa: E402

import leakage as LK  # noqa: E402
import metrics as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
RESULTS_DIR = os.path.join(HERE, "results")

SEED = 42
SAMPLE_A = 200
SAMPLE_B = 200
LEAKAGE_SENSITIVITY_SUBSAMPLE = 100
POTENCY_SUBSAMPLE = 200
KNN_GRID = [(10, 0.5), (10, 1.0), (25, 1.0), (50, 1.0)]
THRESHOLDS = [0.2, 0.3, 0.4, 0.5, 0.6]
PRIMARY_LEAKAGE_CUTOFF = 0.95
SENSITIVITY_LEAKAGE_CUTOFF = 0.90
POTENCY_FLOORS = [("10uM_5.0", 5.0), ("1uM_6.0", 6.0)]


def wilson_ci(successes, n, z=1.96):
    if n == 0:
        return (None, None)
    p = successes / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)) / denom
    return (round(max(0.0, center - half), 4), round(min(1.0, center + half), 4))


def load_everything():
    t0 = time.time()
    print("Loading v1 frozen index...", flush=True)
    full_fps, full_pop, full_df = TF._load()
    print(f"  {len(full_df)} pairs, {full_df['smiles'].nunique()} compounds ({time.time()-t0:.0f}s)", flush=True)

    print("Loading leakage fp2 (topological)...", flush=True)
    fp2 = np.load(os.path.join(DATA_DIR, "leakage_fp2.npz"), allow_pickle=True)
    fp2_smiles = fp2["smiles"]
    fp2_arr = fp2["fps"]
    print(f"  {len(fp2_smiles)} distinct compounds fingerprinted ({time.time()-t0:.0f}s)", flush=True)

    leak_idx = LK.LeakageIndex(full_fps, full_pop, full_df, fp2_smiles, fp2_arr)

    with open(os.path.join(DATA_DIR, "eval_sets.json")) as f:
        eval_sets = json.load(f)

    return full_df, leak_idx, eval_sets


def annotated_universe(row):
    return set(row["annotated_targets"])


def primary_universe(row):
    return set(row["primary_targets"])


def weighted_knn_rank(reduced_df, tanimoto, true_annotated, true_primary, k, alpha):
    """tanimoto: array aligned to reduced_df rows (compound-target pairs).
       Dedup to one row per DISTINCT compound (max tanimoto if it somehow
       differs -- it won't, same fingerprint per compound, but max is the
       safe reduction), take the k highest-similarity compounds, and have
       each vote tanimoto_i^alpha for EVERY target in ITS OWN full
       annotation set within the reduced reference (not just the target
       being scored) -- this is the review's formula-fix (point 5):
       'each neighbour compound votes for every target it is annotated to.'
    """
    tmp = reduced_df.copy()
    tmp["tanimoto"] = tanimoto
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    top_k = by_compound.sort_values("tanimoto", ascending=False).head(k)
    if top_k.empty:
        return None, None
    neighbor_smiles = set(top_k.index)
    neighbor_targets = tmp[tmp["smiles"].isin(neighbor_smiles)][["smiles", "target_chembl"]].drop_duplicates()
    sim_by_smiles = top_k["tanimoto"].to_dict()
    scores = {}
    for row in neighbor_targets.itertuples(index=False):
        w = sim_by_smiles[row.smiles] ** alpha
        scores[row.target_chembl] = scores.get(row.target_chembl, 0.0) + w
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    rank_map = {t: i + 1 for i, (t, _) in enumerate(ranked)}
    best_annot = min((rank_map[t] for t in true_annotated if t in rank_map), default=None)
    best_prim = min((rank_map[t] for t in true_primary if t in rank_map), default=None)
    return best_annot, best_prim


def run_baseline_and_knn(query_rows, leak_idx, cutoff, label, do_knn=True, do_potency=False, potency_subsample_n=0, rng=None):
    n = len(query_rows)
    t0 = time.time()
    any_ranks, primary_ranks = [], []
    pr_by_threshold = {th: [] for th in THRESHOLDS}
    knn_any = {cfg: [] for cfg in KNN_GRID}
    knn_primary = {cfg: [] for cfg in KNN_GRID}
    potency_any = {fl[0]: [] for fl in POTENCY_FLOORS}
    potency_primary = {fl[0]: [] for fl in POTENCY_FLOORS}
    n_potency_done = 0

    potency_targets = set()
    if do_potency and potency_subsample_n:
        idxs = list(range(n))
        (rng or random).shuffle(idxs)
        potency_targets = set(idxs[:potency_subsample_n])

    for i, row in enumerate(query_rows):
        smi = row["smiles"]
        annotated = annotated_universe(row)
        primary = primary_universe(row)

        reduced_fps, reduced_pop, reduced_df = leak_idx.reduced_reference(smi, cutoff)
        q_packed, q_pop, _, _ = leak_idx.query_fingerprints(smi)

        results = TF._search_against(reduced_fps, reduced_pop, reduced_df, q_packed, q_pop,
                                      threshold=0.0, max_compounds_per_target=1)
        any_rank = M.rank_of_best(results, annotated)
        prim_rank = M.rank_of_best(results, primary) if primary else None
        any_ranks.append(any_rank)
        if primary:
            primary_ranks.append(prim_rank)

        for th in THRESHOLDS:
            predicted = [r for r in results if r["best_similarity"] >= th]
            pr_by_threshold[th].append(M.precision_recall_at_threshold(predicted, annotated))

        if do_knn:
            tanimoto = LK.tanimoto_to_all(q_packed, reduced_fps, reduced_pop, q_pop)
            for cfg in KNN_GRID:
                k, alpha = cfg
                a_rank, p_rank = weighted_knn_rank(reduced_df, tanimoto, annotated, primary, k, alpha)
                knn_any[cfg].append(a_rank)
                if primary:
                    knn_primary[cfg].append(p_rank)

        if do_potency and i in potency_targets:
            n_potency_done += 1
            for label_f, floor in POTENCY_FLOORS:
                pot_mask = (reduced_df["pchembl_value"] >= floor).to_numpy()
                pot_fps = reduced_fps[pot_mask]
                pot_pop = reduced_pop[pot_mask]
                pot_df = reduced_df.loc[pot_mask]
                pot_results = TF._search_against(pot_fps, pot_pop, pot_df, q_packed, q_pop,
                                                  threshold=0.0, max_compounds_per_target=1)
                potency_any[label_f].append(M.rank_of_best(pot_results, annotated))
                if primary:
                    potency_primary[label_f].append(M.rank_of_best(pot_results, primary))

        if (i + 1) % 25 == 0 or (i + 1) == n:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            remaining = (n - (i + 1)) / rate if rate > 0 else 0
            print(f"  [{label}] {i+1}/{n} ({elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining)", flush=True)

    n_targets_universe = int(query_rows and 4658 or 0)  # v1's own frozen target-universe size
    out = {
        "label": label,
        "leakage_cutoff": cutoff,
        "n_queries": n,
        "any_annotated": M.any_annotated_metrics(any_ranks),
        "primary_target": M.primary_metrics(primary_ranks),
        "precision_recall_by_threshold": {
            str(th): M.aggregate_precision_recall(pr_by_threshold[th], n_targets_universe) for th in THRESHOLDS
        },
    }
    out["any_annotated"]["top_1_ci95"] = wilson_ci(sum(1 for r in any_ranks if r is not None and r <= 1), n)
    out["any_annotated"]["top_10_ci95"] = wilson_ci(sum(1 for r in any_ranks if r is not None and r <= 10), n)

    if do_knn:
        out["weighted_knn_vs_best_similarity"] = {}
        for cfg in KNN_GRID:
            k, alpha = cfg
            out["weighted_knn_vs_best_similarity"][f"k={k}_alpha={alpha}"] = {
                "any_annotated": M.any_annotated_metrics(knn_any[cfg]),
                "primary_target": M.primary_metrics(knn_primary[cfg]),
            }
    if do_potency:
        out["potency_filter"] = {"n_subsample": n_potency_done}
        for label_f, _ in POTENCY_FLOORS:
            out["potency_filter"][label_f] = {
                "any_annotated": M.any_annotated_metrics(potency_any[label_f]),
                "primary_target": M.primary_metrics(potency_primary[label_f]),
            }
    return out


def popularity_baseline(full_df, all_query_rows):
    popularity = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size()
    ranked_targets = popularity.sort_values(ascending=False).index.tolist()
    rank_map = {t: i + 1 for i, t in enumerate(ranked_targets)}
    any_ranks, primary_ranks = [], []
    for row in all_query_rows:
        annotated = annotated_universe(row)
        primary = primary_universe(row)
        any_ranks.append(min((rank_map[t] for t in annotated if t in rank_map), default=None))
        if primary:
            primary_ranks.append(min((rank_map[t] for t in primary if t in rank_map), default=None))
    return {
        "any_annotated": M.any_annotated_metrics(any_ranks),
        "primary_target": M.primary_metrics(primary_ranks),
    }


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    rng = random.Random(SEED)
    full_df, leak_idx, eval_sets = load_everything()

    set_a = eval_sets["set_A_approved"]
    set_b = eval_sets["set_B_clinical"]
    set_d = eval_sets["set_D_sparse_target"]

    rng.shuffle(set_a)
    rng.shuffle(set_b)
    sample_a = set_a[:SAMPLE_A]
    sample_b = set_b[:SAMPLE_B]

    report = {"methodology": {
        "seed": SEED,
        "sample_a_n": len(sample_a), "sample_a_of_total": len(set_a),
        "sample_b_n": len(sample_b), "sample_b_of_total": len(set_b),
        "set_d_n": len(set_d),
        "primary_leakage_cutoff": PRIMARY_LEAKAGE_CUTOFF,
        "sensitivity_leakage_cutoff": SENSITIVITY_LEAKAGE_CUTOFF,
        "knn_grid": KNN_GRID,
        "thresholds": THRESHOLDS,
        "potency_floors": POTENCY_FLOORS,
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }}

    print("\n=== Set A (approved drugs) baseline @0.95 leakage cutoff, with kNN + potency ===", flush=True)
    report["set_A"] = run_baseline_and_knn(sample_a, leak_idx, PRIMARY_LEAKAGE_CUTOFF, "SetA",
                                            do_knn=True, do_potency=True,
                                            potency_subsample_n=min(POTENCY_SUBSAMPLE, len(sample_a)), rng=rng)
    _save(report)

    print("\n=== Set B (clinical, phase 1-3) baseline @0.95 leakage cutoff, with kNN ===", flush=True)
    report["set_B"] = run_baseline_and_knn(sample_b, leak_idx, PRIMARY_LEAKAGE_CUTOFF, "SetB", do_knn=True)
    _save(report)

    print("\n=== Set D (sparse-target, ALL matched) baseline @0.95 leakage cutoff ===", flush=True)
    report["set_D"] = run_baseline_and_knn(set_d, leak_idx, PRIMARY_LEAKAGE_CUTOFF, "SetD", do_knn=False)
    _save(report)

    print(f"\n=== Leakage sensitivity: same {LEAKAGE_SENSITIVITY_SUBSAMPLE} queries @0.90 cutoff ===", flush=True)
    sens_sample = (sample_a + sample_b)[:LEAKAGE_SENSITIVITY_SUBSAMPLE]
    report["leakage_sensitivity_0_90"] = run_baseline_and_knn(sens_sample, leak_idx, SENSITIVITY_LEAKAGE_CUTOFF,
                                                               "LeakageSensitivity", do_knn=False)
    _save(report)

    print("\n=== Popularity baseline (query-independent floor) ===", flush=True)
    report["popularity_baseline"] = popularity_baseline(full_df, sample_a + sample_b + set_d)
    _save(report)

    print("\nPhase 0 sample run complete.", flush=True)


def _save(report):
    out_path = os.path.join(RESULTS_DIR, "phase0_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"  (partial results saved -> {out_path})", flush=True)


if __name__ == "__main__":
    main()
