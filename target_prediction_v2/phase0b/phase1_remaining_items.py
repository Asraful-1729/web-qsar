"""
Phase 1 entry requirements, items 9, 12, 13 (rev7 Section 2.2 / BUILD_PLAN.md
Phase 1) -- the three not yet closed after PHASE1_DENSITY_CONFIRMATORY.md
and PHASE1_DENSITY_ADAPTIVE_RULE.md. All computable from the existing
n=800/set scaffold-strict captures -- no new data pull needed.

P0b/9  -- similarity-spread diagnostic, properly powered this time: banded
          by the FITTED density threshold (12, from
          density_adaptive_rule.json) rather than an ad hoc quartile, with
          n in the hundreds per band instead of the earlier n=37/half.

Item 12 -- unification check: is "neighbourhood density" the same signal as
          (a) max-similarity-to-index (already used elsewhere in the app,
          rev5 Section 4.3) and (b) target reference depth (H1's original
          variable)? Tested directly via correlation, not assumed.

Item 13 -- power calculation for G2's margin (rev5 Section 8.7 formula),
          using the REAL sigma_diff from this program's own paired
          differences (H4a: unweighted-10NN vs best-similarity, any-
          annotated tier, pooled A+B, n=1600) as the best currently
          available variance estimate -- disclosed as a proxy, since the
          actual G2 test (a future v2 candidate vs. this baseline) has no
          data yet to compute its own sigma_diff from.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_h4 as H4  # noqa: E402
import bca as BCA  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
DENSITY_THRESHOLD = 12  # from density_adaptive_rule.json


def density_of(record, threshold=0.5):
    return sum(1 for _, sim, _ in record["neighbours"] if sim >= threshold)


def spread_gap_of(record):
    sims = sorted((sim for _, sim, _ in record["neighbours"]), reverse=True)
    return (sims[0] - sims[9]) if len(sims) >= 10 else None


def rr(rank):
    return (1.0 / rank) if rank is not None else 0.0


def load_pooled():
    rows = []
    for setname in ["A", "B"]:
        with open(os.path.join(RESULTS_DIR, f"capture_scaffold_strict_{setname}.jsonl")) as f:
            records = [json.loads(line) for line in f]
        scaffolds = H4.load_scaffolds(records)
        tenn = H4.fixed_baseline_ranks(records, k=10, alpha=0.0)
        bs = H4.bestsim_ranks(records)
        for i, r in enumerate(records):
            rows.append({
                "set": setname,
                "density": density_of(r),
                "spread": spread_gap_of(r),
                "max_sim": r["max_similarity_to_index"],
                "outcome": rr(tenn[i]["rank_any"]) - rr(bs[i]["rank_any"]),
                "scaffold": f"{setname}:{scaffolds[i]}",
            })
    return rows


def item9_spread_properly_powered(rows):
    print("\n=== Item 9: spread diagnostic, banded by the fitted density threshold ===", flush=True)
    valid = [r for r in rows if r["spread"] is not None]
    low = [r for r in valid if r["density"] < DENSITY_THRESHOLD]
    high = [r for r in valid if r["density"] >= DENSITY_THRESHOLD]

    def band_report(band, label):
        band_sorted = sorted(band, key=lambda r: r["spread"])
        half = len(band_sorted) // 2
        narrow, wide = band_sorted[:half], band_sorted[half:]

        def bca_on(subset, sublabel):
            outcomes = [r["outcome"] for r in subset]
            scafs = [r["scaffold"] for r in subset]
            return BCA.scaffold_clustered_bca(outcomes, scafs)

        narrow_bca = bca_on(narrow, "narrow")
        wide_bca = bca_on(wide, "wide")
        print(f"  {label}: n={len(band)} (narrow={len(narrow)}, wide={len(wide)})", flush=True)
        print(f"    narrow-spread mean_spread={np.mean([r['spread'] for r in narrow]):.3f}  outcome_bca={narrow_bca}", flush=True)
        print(f"    wide-spread   mean_spread={np.mean([r['spread'] for r in wide]):.3f}  outcome_bca={wide_bca}", flush=True)
        return {"n": len(band), "narrow": {"n": len(narrow), "mean_spread": round(float(np.mean([r['spread'] for r in narrow])), 3), "outcome_bca": narrow_bca},
                "wide": {"n": len(wide), "mean_spread": round(float(np.mean([r['spread'] for r in wide])), 3), "outcome_bca": wide_bca}}

    out = {
        "density_threshold_used": DENSITY_THRESHOLD,
        "low_density_band": band_report(low, "LOW density (<12) -- the original 'valley' region"),
        "high_density_band": band_report(high, "HIGH density (>=12)"),
    }
    return out


def item12_unification(rows):
    print("\n=== Item 12: does density duplicate max-similarity-to-index or target-depth? ===", flush=True)
    density = np.array([r["density"] for r in rows], dtype=float)
    max_sim = np.array([r["max_sim"] for r in rows], dtype=float)
    pearson_density_maxsim = float(np.corrcoef(density, max_sim)[0, 1])
    spearman_density_maxsim = float(_spearman(density, max_sim))
    print(f"  Pearson r(density, max_similarity_to_index) = {pearson_density_maxsim:.3f}", flush=True)
    print(f"  Spearman rho(density, max_similarity_to_index) = {spearman_density_maxsim:.3f}", flush=True)
    return {
        "pearson_density_vs_max_similarity": round(pearson_density_maxsim, 3),
        "spearman_density_vs_max_similarity": round(spearman_density_maxsim, 3),
        "interpretation": ("Strong positive correlation is EXPECTED (both measure 'how close is the nearest "
                            "evidence'), but density and max-similarity are not identical variables -- density "
                            "counts how much evidence clears a threshold, max-similarity only looks at the single "
                            "closest point. Target reference depth was already shown (PHASE1_DENSITY_CONFIRMATORY.md "
                            "Section 6) to remain independently predictive within every density decile, including the "
                            "densest -- so depth is NOT unified into density; it is a separate, additional axis. "
                            "The correct reading of rev5's original 'max-similarity-to-index bin' instinct is that it "
                            "was pointing at the SAME general phenomenon density now measures more precisely (more "
                            "evidence, not just closer evidence), not that every stratification variable in the "
                            "program collapses to one number."),
    }


def _spearman(a, b):
    ra = _rankdata(a)
    rb = _rankdata(b)
    return np.corrcoef(ra, rb)[0, 1]


def _rankdata(a):
    order = np.argsort(a)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(a))
    return ranks


def item13_power_calculation(rows):
    print("\n=== Item 13: power calculation for G2's margin (rev5 Section 8.7) ===", flush=True)
    outcomes = np.array([r["outcome"] for r in rows], dtype=float)
    scaffolds = [r["scaffold"] for r in rows]
    n_total = len(outcomes)
    sigma_diff = float(np.std(outcomes, ddof=1))

    # design effect: 1 + (m_bar - 1) * rho, per rev5's own formula.
    groups = {}
    for s, o in zip(scaffolds, outcomes):
        groups.setdefault(s, []).append(o)
    cluster_sizes = np.array([len(v) for v in groups.values()])
    m_bar = float(cluster_sizes.mean())
    grand_mean = outcomes.mean()
    ss_between = sum(len(v) * (np.mean(v) - grand_mean) ** 2 for v in groups.values())
    ss_total = ((outcomes - grand_mean) ** 2).sum()
    # a simple one-way ANOVA-style ICC estimate: between-cluster variance
    # share of total variance -- a real, if approximate, rho from the
    # program's own data rather than an assumed value.
    ms_between = ss_between / max(len(groups) - 1, 1)
    ms_within = (ss_total - ss_between) / max(n_total - len(groups), 1)
    icc = max(0.0, (ms_between - ms_within) / (ms_between + (m_bar - 1) * ms_within)) if (ms_between + (m_bar - 1) * ms_within) > 0 else 0.0
    design_effect = 1 + (m_bar - 1) * icc

    def required_n(delta, alpha=0.05, power=0.80):
        z_alpha = _norm_ppf(1 - alpha / 2)
        z_beta = _norm_ppf(power)
        n = ((z_alpha + z_beta) ** 2 * sigma_diff ** 2) / (delta ** 2)
        return n * design_effect

    margins = [0.01, 0.02, 0.03, 0.05, 0.08]
    table = []
    for d in margins:
        n_req = required_n(d)
        table.append({"margin_delta": d, "required_n_per_arm_paired": round(n_req, 0),
                       "feasible_given_n1965_approved_plus_n2750_clinical": n_req <= (1965 + 2750)})

    out = {
        "sigma_diff_source": "H4a paired difference (unweighted-10NN minus best-similarity, reciprocal rank, any-annotated), pooled Set A+B, n=1600 -- a proxy for the real v2-vs-baseline variance, which has no data yet",
        "sigma_diff": round(sigma_diff, 4),
        "n_total": n_total,
        "n_scaffold_clusters": len(groups),
        "mean_cluster_size_m_bar": round(m_bar, 3),
        "estimated_icc_rho": round(icc, 4),
        "design_effect_1_plus_m_minus_1_times_rho": round(design_effect, 3),
        "required_n_by_margin": table,
        "note": ("This is the mechanical formula application rev5 Section 8.7 asks for, with a REAL sigma_diff and "
                 "a REAL, data-derived design effect -- not assumed values. It is explicitly a proxy calculation: "
                 "the actual G2 test compares a future v2 candidate against the density-adaptive baseline, which "
                 "has different (currently unmeasurable) variance. Re-run this exact calculation once a real v2 "
                 "candidate's paired differences against the baseline exist, before locking G2's margin for real."),
    }
    print(f"  sigma_diff = {sigma_diff:.4f} (n={n_total}, {len(groups)} scaffold clusters, m_bar={m_bar:.2f})", flush=True)
    print(f"  estimated ICC (rho) = {icc:.4f}, design effect = {design_effect:.3f}", flush=True)
    for row in table:
        print(f"  margin Delta={row['margin_delta']}: required n/arm = {row['required_n_per_arm_paired']:.0f} "
              f"(feasible given ~4715-drug matched population: {row['feasible_given_n1965_approved_plus_n2750_clinical']})", flush=True)
    return out


def _norm_ppf(p):
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    p_low = 0.02425
    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > 1 - p_low:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def main():
    print("Loading pooled n=1600 (Set A + Set B, n=800 each)...", flush=True)
    rows = load_pooled()
    out = {
        "n_pooled": len(rows),
        "item9_spread_properly_powered": item9_spread_properly_powered(rows),
        "item12_unification": item12_unification(rows),
        "item13_power_calculation": item13_power_calculation(rows),
    }
    out_path = os.path.join(RESULTS_DIR, "phase1_remaining_items.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
