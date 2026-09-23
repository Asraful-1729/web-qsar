"""
Research-handoff item #2: a genuine SEA (Similarity Ensemble Approach,
Keiser et al. 2007) implementation -- previously marked "not attempted"
in BUILD_PLAN.md (3B), disclosed rather than faked, because a correct
implementation requires real statistical engineering: raw Tanimoto-sum
scores are only meaningful once converted to a significance value against
a target-size-specific NULL distribution (an Extreme Value / Gumbel fit,
the same logic as BLAST E-values) -- SEA is NOT just "sum up similarities
and rank," that part alone (no null correction) is close to what
popularity_correct() tried to fix and what got rejected
(PHASE3A_POPULARITY_CORRECTION.md) for being too crude.

Raw score rule (disclosed choice, no arbitrary inclusion threshold):
raw_score(Q, T) = sum of Tanimoto(Q, l) over every distinct ligand l
indexed against target T (v2_index, matching the "active set" definition
used everywhere else in this program). Fully vectorized (fast packed-bit
Tanimoto), no thresholding -- simpler and avoids an unvalidated threshold
choice; the null-distribution fit uses the identical rule, so this is
self-consistent regardless of whether it matches the original paper's
exact scoring convention.

Null distribution: for each target-SIZE bin (log-spaced; SEA excludes
targets with <5 ligands "by construction", per BUILD_PLAN.md's own note),
draw many random (query, random-same-size-ligand-set) pairs from the same
reference pool and compute the same raw score -- this is the actual
background model SEA needs. Fit scipy.stats.gumbel_r per bin.

Usage: python3 build_sea_null_distributions.py
"""
import json
import os
import random
import time

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2_DIR = os.path.join(HERE, "..", "phase2")
RESULTS_DIR = os.path.join(HERE, "results")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

MIN_LIGANDS = 5  # SEA's own exclusion rule
SIZE_BINS = [(5, 10), (10, 20), (20, 40), (40, 80), (80, 160), (160, 320),
             (320, 640), (640, 1280), (1280, 999999)]
N_NULL_QUERIES = 80
N_SUBSETS_PER_QUERY_PER_BIN = 15  # random same-size ligand-set draws per null query, per bin
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

    # distinct-compound view (first occurrence per smiles) -- the reference
    # pool for both null-query sampling and null ligand-set sampling
    first_idx = df.reset_index().drop_duplicates(subset="smiles", keep="first")["index"].to_numpy()
    print(f"{len(df)} total rows, {len(first_idx)} distinct compounds", flush=True)

    rng = random.Random(SEED)
    np_rng = np.random.default_rng(SEED)

    null_samples = {b: [] for b in SIZE_BINS}
    t0 = time.time()
    query_rows = np_rng.choice(first_idx, size=N_NULL_QUERIES, replace=False)

    for qi, row in enumerate(query_rows):
        q_packed, q_pop = fps[row], int(pop[row])
        sims = tanimoto_to_all(q_packed, fps[first_idx], pop[first_idx], q_pop)
        # exclude the query's own trivial self-similarity (=1.0) from the pool it draws from
        pool_mask = first_idx != row
        pool_sims = sims[pool_mask]

        for (lo, hi) in SIZE_BINS:
            size = min(int((lo + min(hi, 2000)) / 2), len(pool_sims))  # representative size for this bin
            size = max(size, lo)
            for _ in range(N_SUBSETS_PER_QUERY_PER_BIN):
                subset = np_rng.choice(len(pool_sims), size=min(size, len(pool_sims)), replace=False)
                null_samples[(lo, hi)].append(float(pool_sims[subset].sum()))

        if (qi + 1) % 20 == 0:
            print(f"  {qi+1}/{N_NULL_QUERIES} null queries processed ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nFitting Gumbel EVD per size bin...", flush=True)
    fitted = {}
    for b, samples in null_samples.items():
        samples = np.array(samples)
        loc, scale = stats.gumbel_r.fit(samples)
        fitted[f"{b[0]}-{b[1]}"] = {
            "n_null_samples": len(samples), "loc": float(loc), "scale": float(scale),
            "null_mean": float(samples.mean()), "null_std": float(samples.std()),
        }
        print(f"  bin {b}: n={len(samples)}, loc={loc:.3f}, scale={scale:.3f}, "
              f"mean={samples.mean():.3f}, std={samples.std():.3f}", flush=True)

    out = {"size_bins": SIZE_BINS, "min_ligands": MIN_LIGANDS, "seed": SEED,
           "n_null_queries": N_NULL_QUERIES, "fitted_gumbel_params": fitted}
    out_path = os.path.join(RESULTS_DIR, "sea_null_distributions.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path} ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
