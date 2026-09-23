"""
Rev 7 response, flags 2 and 3 -- PRELIMINARY exploratory checks only, on the
same 200-query Set A sample already captured. Rev 7 explicitly scopes the
confirmatory versions of both to Phase 1's full benchmark (decision #18);
this script does NOT attempt that -- it is the cheap same-day check Rev 7
itself suggested ("before committing to 'density alone'... this is a
same-day addition... not a new data collection effort"), to give early
signal before Phase 1 commits to a functional form.

Flag-2: within-pool similarity SPREAD (shape), not just density (count).
  spread_gap  = similarity(1st neighbour) - similarity(10th neighbour)
  spread_var  = variance of the top-10 neighbours' similarities
  Checked against the per-query paired outcome (reciprocal-rank of
  unweighted-10NN minus reciprocal-rank of best-similarity, any-annotated
  tier) via simple correlation, both overall and restricted to the Q2
  density band where the valley was found.

Flag-3: joint density + target-reference-depth model.
  OLS (via numpy least-squares, no statsmodels dependency) of the same
  paired outcome on [density, mean_target_reference_depth, their
  interaction], to check whether depth's coefficient survives once density
  is in the model -- the univariate quartile check alone (diagnose_A_vs_B.py)
  cannot distinguish "depth matters independently" from "depth is just
  correlated with density."

Both use capture_scaffold_strict_A.jsonl (the leakage control the H4a/H4b
conclusions rest on) and the same v1-index reference-depth lookup as
diagnose_A_vs_B.py.

Usage: python3 test_rev7_flags.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))

import target_fishing as TF  # noqa: E402
import score as S  # noqa: E402
import analyze_h4 as H4  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def ols_with_stats(X, y):
    """Plain OLS via least squares, with a manual t-stat/p-value
       computation (two-sided, using the normal approximation -- adequate
       for an exploratory n=200 check; not claiming this replaces a
       properly powered Phase 1 regression)."""
    n, p = X.shape
    Xd = np.column_stack([np.ones(n), X])
    beta, _, _, _ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = n - Xd.shape[1]
    sigma2 = (resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(Xd.T @ Xd)
    se = np.sqrt(np.diag(cov))
    tstats = beta / se
    from math import erf, sqrt
    pvals = [2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2)))) for t in tstats]
    r2 = 1 - (resid @ resid) / ((y - y.mean()) @ (y - y.mean()))
    return {"coef": beta.tolist(), "se": se.tolist(), "t": tstats.tolist(), "p": pvals, "r2": float(r2), "n": n}


def main():
    with open(os.path.join(RESULTS_DIR, "capture_scaffold_strict_A.jsonl")) as f:
        records = [json.loads(line) for line in f]

    _, _, full_df = TF._load()
    depth_by_target = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size().to_dict()

    tenn_ranks = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
    bs_ranks = H4.bestsim_ranks(records)

    rows = []
    for i, r in enumerate(records):
        neighbours = r["neighbours"]
        sims_sorted = sorted((sim for _, sim, _ in neighbours), reverse=True)
        top10 = sims_sorted[:10]
        density_05 = sum(1 for s in sims_sorted if s >= 0.5)
        spread_gap = (top10[0] - top10[-1]) if len(top10) >= 10 else None
        spread_var = float(np.var(top10)) if len(top10) >= 2 else None

        annotated = r["annotated_targets"]
        depths = [depth_by_target.get(t, 0) for t in annotated]
        mean_depth = float(np.mean(depths)) if depths else None

        outcome = rr(tenn_ranks[i]["rank_any"]) - rr(bs_ranks[i]["rank_any"])

        rows.append({
            "density_0.5": density_05, "spread_gap": spread_gap, "spread_var": spread_var,
            "mean_depth": mean_depth, "outcome_rr_diff": outcome,
        })

    # ---- Flag 2: spread vs outcome, overall and within the Q2 density band ----
    valid = [row for row in rows if row["spread_gap"] is not None]
    dens_sorted = sorted(row["density_0.5"] for row in valid)
    n = len(dens_sorted)
    q2_lo, q2_hi = dens_sorted[n // 4], dens_sorted[n // 2]
    q2_rows = [row for row in valid if q2_lo <= row["density_0.5"] <= q2_hi]

    def corr(key_x, key_y, subset):
        xs = np.array([row[key_x] for row in subset])
        ys = np.array([row[key_y] for row in subset])
        if len(xs) < 3 or xs.std() == 0 or ys.std() == 0:
            return None
        return float(np.corrcoef(xs, ys)[0, 1])

    flag2 = {
        "n_overall": len(valid), "n_q2_band": len(q2_rows),
        "q2_density_range": [q2_lo, q2_hi],
        "corr_density_vs_outcome_overall": corr("density_0.5", "outcome_rr_diff", valid),
        "corr_spread_gap_vs_outcome_overall": corr("spread_gap", "outcome_rr_diff", valid),
        "corr_spread_var_vs_outcome_overall": corr("spread_var", "outcome_rr_diff", valid),
        "corr_spread_gap_vs_outcome_within_q2": corr("spread_gap", "outcome_rr_diff", q2_rows),
        "corr_spread_var_vs_outcome_within_q2": corr("spread_var", "outcome_rr_diff", q2_rows),
        "mean_spread_gap_q2": float(np.mean([row["spread_gap"] for row in q2_rows])),
        "mean_outcome_q2": float(np.mean([row["outcome_rr_diff"] for row in q2_rows])),
    }

    # ---- Flag 3: joint density + depth model ----
    valid3 = [row for row in rows if row["mean_depth"] is not None]
    X = np.array([[row["density_0.5"], row["mean_depth"], row["density_0.5"] * row["mean_depth"]] for row in valid3])
    # standardize columns for interpretable, comparable coefficients
    X_std = (X - X.mean(axis=0)) / X.std(axis=0)
    y = np.array([row["outcome_rr_diff"] for row in valid3])
    model_full = ols_with_stats(X_std, y)

    X_density_only = X_std[:, [0]]
    model_density_only = ols_with_stats(X_density_only, y)

    flag3 = {
        "n": len(valid3),
        "model_density_plus_depth_plus_interaction": {
            "labels": ["intercept", "density", "depth", "density_x_depth"],
            "coef": model_full["coef"], "p": model_full["p"], "r2": model_full["r2"],
        },
        "model_density_only": {
            "labels": ["intercept", "density"],
            "coef": model_density_only["coef"], "p": model_density_only["p"], "r2": model_density_only["r2"],
        },
        "interpretation_note": ("depth's p-value in the full model, controlling for density, is the direct test "
                                 "of whether depth adds anything beyond density; r2 improvement from "
                                 "density_only -> full model quantifies how much."),
    }

    out = {"flag2_similarity_spread": flag2, "flag3_joint_density_depth_model": flag3}
    out_path = os.path.join(RESULTS_DIR, "rev7_flags_report.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
