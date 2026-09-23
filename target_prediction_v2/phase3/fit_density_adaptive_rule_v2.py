"""
Phase 3A, the re-fit BUILD_PLAN.md's own text requires before the
density-adaptive rule can be used for anything against the Phase 2
rebuilt index: identical methodology to phase0b/fit_density_adaptive_rule.py
(isotonic regression + scaffold-clustered bootstrap band, plus the smooth
logistic mixing-weight approximation), run against the new
capture_scaffold_strict_v2_{A,B}.jsonl captures (n=800/set, against
v2_index, 3,973,438 native-human pairs vs. v1's 1,312,849).

Only difference from the original script: scaffold lookup is pointed at
phase2/v2_index instead of v1's frozen index (analyze_h4.py's
load_scaffolds() hardcodes the v1 path -- reimplemented here rather than
patched in place, to avoid touching a script other Phase 0b results still
depend on). fixed_baseline_ranks/bestsim_ranks are pure functions of each
record's own embedded fields (no index dependency) -- reused unchanged.

Usage: python3 fit_density_adaptive_rule_v2.py
"""
import json
import os
import sys

import numpy as np
from sklearn.isotonic import IsotonicRegression
from scipy.optimize import curve_fit

PHASE0B_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0b")
PHASE2_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase2")
sys.path.insert(0, PHASE0B_DIR)
sys.path.insert(0, os.path.join(PHASE0B_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ["TARGET_FISHING_INDEX_DIR"] = os.path.abspath(os.path.join(PHASE2_DIR, "v2_index"))

import analyze_h4 as H4  # noqa: E402
import bca as BCA  # noqa: E402
import target_fishing as TF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

DENSITY_GRID = [0, 1, 2, 3, 5, 8, 12, 18, 25, 35, 50, 70, 100, 150, 200, 300]
N_BOOT = 2000


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_scaffolds_v2(records):
    _, _, full_df = TF._load()
    needed = set(r["smiles"] for r in records)
    sub = full_df[full_df["smiles"].isin(needed)][["smiles", "murcko_scaffold"]].drop_duplicates(subset="smiles")
    scaffold_by_smiles = dict(zip(sub["smiles"], sub["murcko_scaffold"]))
    out = []
    for r in records:
        sc = scaffold_by_smiles.get(r["smiles"])
        if not isinstance(sc, str) or not sc:
            sc = f"__no_ring__:{r['smiles']}"
        out.append(sc)
    return out


def load_set(setname):
    with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_v2_{setname}.jsonl")) as f:
        records = [json.loads(line) for line in f]
    scaffolds = load_scaffolds_v2(records)
    tenn = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs = H4.bestsim_ranks(records)
    rows = []
    for i, r in enumerate(records):
        rows.append({
            "density": density_of(r),
            "outcome": rr(tenn[i]["rank_any"]) - rr(bs[i]["rank_any"]),
            "scaffold": f"{setname}:{scaffolds[i]}",
        })
    return rows


def isotonic_fit(densities, outcomes):
    ir = IsotonicRegression(out_of_bounds="clip", increasing=True)
    ir.fit(densities, outcomes)
    return ir


def main():
    print("Loading n=800 v2 captures for Set A and Set B, pooling...", flush=True)
    rows = load_set("A") + load_set("B")
    densities = np.array([r["density"] for r in rows], dtype=float)
    outcomes = np.array([r["outcome"] for r in rows], dtype=float)
    scaffolds = [r["scaffold"] for r in rows]
    n = len(rows)
    print(f"  pooled n={n} queries, {len(set(scaffolds))} distinct scaffold groups", flush=True)
    print(f"  density: min={densities.min():.0f} max={densities.max():.0f} "
          f"median={np.median(densities):.0f} mean={densities.mean():.1f}", flush=True)

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
    stays_significant = all(row["ci_excludes_zero"] for row in isotonic_table
                             if row["density"] >= (threshold_density or 1e9))

    print("Fitting smooth logistic mixing-weight approximation...", flush=True)

    def logistic(x, d0, s, lo, hi):
        return lo + (hi - lo) / (1.0 + np.exp(-(x - d0) / s))

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
        "index": "Phase 2 rebuilt (v2_index, 3,973,438 native-human pairs)",
        "n_pooled_queries": n, "n_scaffold_groups": len(group_keys),
        "density_grid": DENSITY_GRID,
        "isotonic_fit_with_scaffold_clustered_95pct_band": isotonic_table,
        "recommended_pooling_threshold_density": threshold_density,
        "threshold_holds_for_all_higher_grid_points": stays_significant,
        "logistic_mixing_weight_fit": {"formula": "w(density) = floor + (ceiling - floor) / (1 + exp(-(density - d0_midpoint_density) / s_scale))",
                                        "params": logistic_params, "table": logistic_table},
    }
    out_path = os.path.join(RESULTS_DIR, "density_adaptive_rule_v2.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str), flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
