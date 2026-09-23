"""
Phase 1 entry requirement #4 (rev7 Section 2.2 item 4 / BUILD_PLAN.md Phase 1
item 11): fit an actual continuous density-adaptive retrieval rule from the
n=800/set confirmatory data, replacing the hand-set quartile thresholds the
Phase 0b exploration used to locate the shape.

Per PHASE1_DENSITY_CONFIRMATORY.md's finding: density (not which eval set a
query came from) is the driver, and the effect is present in both sets --
so this pools Set A + Set B (1600 queries) rather than fitting two separate
per-set rules, matching that finding directly rather than re-introducing a
set-specific split the confirmatory pass just argued against.

Two fits, deliberately kept separate:

1. ISOTONIC regression (sklearn, monotonic non-decreasing) of the paired
   outcome (reciprocal-rank of unweighted-10NN minus reciprocal-rank of
   best-similarity, any-annotated tier -- the tier with full n=1600) on
   density. This is the rigorous, assumption-light nonparametric fit --
   it assumes ONLY that expected benefit does not decrease as density
   increases (a real assumption, weaker than a specific functional form,
   and consistent with every confirmed finding so far now that the one
   non-monotonic claim, the Q2 valley, has been retracted). A scaffold-
   clustered bootstrap (reusing bca.py's cluster-resampling machinery,
   refitting isotonic regression on each resample) gives a confidence
   band on the fitted curve itself, not just point estimates.

2. A smooth LOGISTIC mixing weight w(density) = 1 / (1 + exp(-(density -
   d0) / s)), least-squares fit to the raw per-query outcomes. This is
   the practically usable formula for a production mixing rule (Rev 7's
   own suggested form: "the pooling/best-similarity mixing weight as a
   smooth function of local density") -- an explicit, disclosed
   parametric simplification of the isotonic fit, not a replacement for
   it; reported alongside the isotonic curve so a reader can see how well
   the smooth approximation tracks the rigorous one.

From the isotonic fit + its bootstrap band, a simple threshold rule is
also derived: the smallest density at which the fitted curve's lower
confidence bound first crosses (and stays above) zero -- the density
level above which pooling can be recommended with actual statistical
support, not just a positive point estimate.

Usage: python3 fit_density_adaptive_rule.py
"""
import json
import os
import sys

import numpy as np
from sklearn.isotonic import IsotonicRegression
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_h4 as H4  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

DENSITY_GRID = [0, 1, 2, 3, 5, 8, 12, 18, 25, 35, 50, 70, 100, 150, 200, 300]
N_BOOT = 2000  # isotonic refit per resample is more expensive than a mean -- kept below the 10k default


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_set(setname):
    with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_{setname}.jsonl")) as f:
        records = [json.loads(line) for line in f]
    scaffolds = H4.load_scaffolds(records)
    tenn = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs = H4.bestsim_ranks(records)
    rows = []
    for i, r in enumerate(records):
        rows.append({
            "density": density_of(r),
            "outcome": rr(tenn[i]["rank_any"]) - rr(bs[i]["rank_any"]),
            "scaffold": f"{setname}:{scaffolds[i]}",  # keep sets' scaffold groups distinct when pooling
        })
    return rows


def isotonic_fit(densities, outcomes):
    ir = IsotonicRegression(out_of_bounds="clip", increasing=True)
    ir.fit(densities, outcomes)
    return ir


def main():
    print("Loading n=800 captures for Set A and Set B, pooling (density is the confirmed driver)...", flush=True)
    rows = load_set("A") + load_set("B")
    densities = np.array([r["density"] for r in rows], dtype=float)
    outcomes = np.array([r["outcome"] for r in rows], dtype=float)
    scaffolds = [r["scaffold"] for r in rows]
    n = len(rows)
    print(f"  pooled n={n} queries, {len(set(scaffolds))} distinct scaffold groups", flush=True)

    # ---- 1. Isotonic fit + scaffold-clustered bootstrap band ----
    print("Fitting isotonic regression + scaffold-clustered bootstrap band...", flush=True)
    point_fit = isotonic_fit(densities, outcomes)
    point_curve = point_fit.predict(DENSITY_GRID)

    groups = {}
    for i, g in enumerate(scaffolds):
        groups.setdefault(g, []).append(i)
    group_keys = list(groups.keys())
    n_groups = len(group_keys)

    import random
    rng = random.Random(42)
    boot_curves = []
    for b in range(N_BOOT):
        resample_groups = [group_keys[rng.randrange(n_groups)] for _ in range(n_groups)]
        idx = []
        for g in resample_groups:
            idx.extend(groups[g])
        idx = np.array(idx)
        try:
            ir_b = isotonic_fit(densities[idx], outcomes[idx])
            boot_curves.append(ir_b.predict(DENSITY_GRID))
        except Exception:
            continue
        if (b + 1) % 500 == 0:
            print(f"  bootstrap {b+1}/{N_BOOT}", flush=True)
    boot_curves = np.array(boot_curves)
    ci_low = np.percentile(boot_curves, 2.5, axis=0)
    ci_high = np.percentile(boot_curves, 97.5, axis=0)

    isotonic_table = []
    threshold_density = None
    for d, pt, lo, hi in zip(DENSITY_GRID, point_curve, ci_low, ci_high):
        sig = lo > 0
        isotonic_table.append({"density": d, "fitted_benefit": round(float(pt), 4),
                                "ci_low": round(float(lo), 4), "ci_high": round(float(hi), 4),
                                "ci_excludes_zero": bool(sig)})
        if sig and threshold_density is None:
            threshold_density = d
    # confirm the threshold "sticks" (every grid point at or above it is also significant) --
    # report the raw table regardless so a reader can see if it doesn't.
    stays_significant = all(row["ci_excludes_zero"] for row in isotonic_table
                             if row["density"] >= (threshold_density or 1e9))

    # ---- 2. Smooth logistic mixing-weight approximation ----
    print("Fitting smooth logistic mixing-weight approximation...", flush=True)

    def logistic(x, d0, s, lo, hi):
        return lo + (hi - lo) / (1.0 + np.exp(-(x - d0) / s))

    # initial guesses: midpoint near the isotonic threshold, floor/ceiling from observed outcome range
    p0 = [threshold_density or 20.0, 10.0, float(np.min(point_curve)), float(np.max(point_curve))]
    try:
        popt, _ = curve_fit(logistic, densities, outcomes, p0=p0, maxfev=20000)
        d0, s, lo, hi = popt
        logistic_table = [{"density": d, "fitted": round(float(logistic(d, *popt)), 4)} for d in DENSITY_GRID]
        logistic_params = {"d0_midpoint_density": round(float(d0), 2), "s_scale": round(float(s), 2),
                            "floor": round(float(lo), 4), "ceiling": round(float(hi), 4)}
    except Exception as e:
        logistic_table = None
        logistic_params = {"error": str(e)}

    out = {
        "n_pooled_queries": n, "n_scaffold_groups": len(group_keys),
        "density_grid": DENSITY_GRID,
        "isotonic_fit_with_scaffold_clustered_95pct_band": isotonic_table,
        "recommended_pooling_threshold_density": threshold_density,
        "threshold_holds_for_all_higher_grid_points": stays_significant,
        "logistic_mixing_weight_fit": {"formula": "w(density) = floor + (ceiling - floor) / (1 + exp(-(density - d0_midpoint_density) / s_scale))",
                                        "params": logistic_params, "table": logistic_table},
    }
    out_path = os.path.join(RESULTS_DIR, "density_adaptive_rule.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str), flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
