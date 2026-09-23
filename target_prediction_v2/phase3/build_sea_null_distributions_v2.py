"""
Fix for a real bug found in build_sea_null_distributions.py (v1): the top
size bin (1280-999999 ligands) used a single fixed "representative size"
(~1640) for null generation, but real targets in the v2 index range up to
195,809 ligands (63 targets >10,000; 15 targets >50,000; confirmed via
v2_index/compounds.csv.gz). Scoring a 106,614-ligand or 114,949-ligand
target's raw score against a null fit for ~1640 ligands produces absurd
z-scores (z=636, z=641 observed in diagnostics) purely from scale mismatch,
not genuine relevance -- these mega-promiscuous targets then dominate SEA's
ranked output for nearly every query, explaining SEA's implausible
catastrophic-underperformance result in test_sea_vs_v2.py.

Fix: replace discrete capped bins with a CONTINUOUS null model. Null
raw-score mean and std are sampled at many representative sizes spanning
the full real range (up to ~200,000, covering the true max), then fit as
power laws of size: mean(N) = a * N^b, std(N) = c * N^d (log-log linear
regression). At scoring time this gives a principled null for ANY real
target size, not just the sizes that happened to get sampled.

Empirical justification for the power-law form: in the original bin fit,
null_mean/size was ~constant (~0.118) across all 9 bins -- consistent with
b~1 (mean scales linearly with size, as expected for a sum over N roughly
similarly-distributed terms). null_std/sqrt(size) was NOT constant (rising
from 0.048 to 0.558 across the bins) -- consistent with d>0.5, i.e. terms
are positively correlated (chemical-space clustering among reference
compounds), so std grows faster than sqrt(N). A single power-law fit
handles both cleanly without needing to know the mechanism.

Usage: python3 build_sea_null_distributions_v2.py
"""
import json
import os
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2_DIR = os.path.join(HERE, "..", "phase2")
RESULTS_DIR = os.path.join(HERE, "results")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

MIN_LIGANDS = 5
# log-spaced representative sizes spanning the FULL real range observed in
# v2_index (max ligand-set size = 195,809) -- not capped bins.
SIZES = [7, 15, 30, 60, 120, 240, 480, 960, 1920, 3840, 7680, 15360,
          30720, 61440, 122880, 196608]
N_NULL_QUERIES = 80
N_SUBSETS_PER_QUERY_PER_SIZE = 15
SEED = 42


def tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop):
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)


def main():
    print("Loading v2_index...", flush=True)
    df = pd.read_csv(os.path.join(PHASE2_DIR, "v2_index", "compounds.csv.gz"))
    fps = np.load(os.path.join(PHASE2_DIR, "v2_index", "fingerprints.npz"))["fps"]
    pop = _POPCOUNT_TABLE[fps].sum(axis=1)

    first_idx = df.reset_index().drop_duplicates(subset="smiles", keep="first")["index"].to_numpy()
    print(f"{len(df)} total rows, {len(first_idx)} distinct compounds", flush=True)

    real_sizes = df.groupby("target_chembl")["smiles"].nunique()
    print(f"Real target ligand-set sizes: min={real_sizes.min()}, max={real_sizes.max()}, "
          f"n>10000={int((real_sizes>10000).sum())}, n>50000={int((real_sizes>50000).sum())}", flush=True)

    np_rng = np.random.default_rng(SEED)

    null_samples = {s: [] for s in SIZES}
    t0 = time.time()
    query_rows = np_rng.choice(first_idx, size=N_NULL_QUERIES, replace=False)

    for qi, row in enumerate(query_rows):
        q_packed, q_pop = fps[row], int(pop[row])
        sims = tanimoto_to_all(q_packed, fps[first_idx], pop[first_idx], q_pop)
        pool_mask = first_idx != row
        pool_sims = sims[pool_mask]
        n_pool = len(pool_sims)

        for size in SIZES:
            draw_size = min(size, n_pool)
            for _ in range(N_SUBSETS_PER_QUERY_PER_SIZE):
                subset = np_rng.choice(n_pool, size=draw_size, replace=False)
                null_samples[size].append(float(pool_sims[subset].sum()))

        if (qi + 1) % 10 == 0:
            print(f"  {qi+1}/{N_NULL_QUERIES} null queries processed ({time.time()-t0:.0f}s)", flush=True)

    print("\nRaw mean/std per representative size:", flush=True)
    size_arr, mean_arr, std_arr = [], [], []
    for s in SIZES:
        samples = np.array(null_samples[s])
        m, sd = samples.mean(), samples.std()
        size_arr.append(s); mean_arr.append(m); std_arr.append(sd)
        print(f"  size={s}: n={len(samples)}, mean={m:.3f}, std={sd:.3f}, mean/size={m/s:.4f}, std/sqrt(size)={sd/np.sqrt(s):.4f}", flush=True)

    # power-law fit via log-log linear regression: log(y) = log(a) + b*log(size)
    log_size = np.log(np.array(size_arr, dtype=float))
    log_mean = np.log(np.array(mean_arr, dtype=float))
    log_std = np.log(np.array(std_arr, dtype=float))

    b_mean, loga_mean = np.polyfit(log_size, log_mean, 1)
    b_std, loga_std = np.polyfit(log_size, log_std, 1)
    a_mean, a_std = float(np.exp(loga_mean)), float(np.exp(loga_std))

    # goodness of fit (R^2 in log-log space)
    pred_mean = loga_mean + b_mean * log_size
    pred_std = loga_std + b_std * log_size
    r2_mean = 1 - np.sum((log_mean - pred_mean) ** 2) / np.sum((log_mean - log_mean.mean()) ** 2)
    r2_std = 1 - np.sum((log_std - pred_std) ** 2) / np.sum((log_std - log_std.mean()) ** 2)

    print(f"\nPower-law fit: mean(N) = {a_mean:.5f} * N^{b_mean:.4f}  (log-log R^2={r2_mean:.5f})", flush=True)
    print(f"Power-law fit: std(N)  = {a_std:.5f} * N^{b_std:.4f}  (log-log R^2={r2_std:.5f})", flush=True)

    # sanity check against real max size
    real_max = int(real_sizes.max())
    check_mean = a_mean * real_max ** b_mean
    check_std = a_std * real_max ** b_std
    print(f"\nSanity check at real max size N={real_max}: predicted mean={check_mean:.1f}, std={check_std:.1f}", flush=True)

    out = {
        "method": "continuous power-law null model (v2 fix for mega-target scale mismatch bug in v1)",
        "sampled_sizes": SIZES,
        "min_ligands": MIN_LIGANDS,
        "seed": SEED,
        "n_null_queries": N_NULL_QUERIES,
        "n_subsets_per_query_per_size": N_SUBSETS_PER_QUERY_PER_SIZE,
        "raw_fit_points": [
            {"size": s, "mean": float(m), "std": float(sd)}
            for s, m, sd in zip(size_arr, mean_arr, std_arr)
        ],
        "power_law": {
            "mean_a": a_mean, "mean_b": float(b_mean), "mean_r2": float(r2_mean),
            "std_a": a_std, "std_b": float(b_std), "std_r2": float(r2_std),
        },
        "real_size_range_observed": {"min": int(real_sizes.min()), "max": real_max},
        "note": (
            "For a target with real ligand-set size N, compute null_mean = mean_a * N**mean_b, "
            "null_std = std_a * N**std_b, then Gumbel scale = null_std / 1.2825 (pi/sqrt(6)), "
            "loc = null_mean - 0.5772*scale (Euler-Mascheroni). Valid for any N in-range; "
            "replaces the discrete-bin v1 approach that broke down for the largest, most "
            "promiscuous targets (up to 195,809 ligands)."
        ),
    }
    out_path = os.path.join(RESULTS_DIR, "sea_null_distributions_v2.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path} ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
