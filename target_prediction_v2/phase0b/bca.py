"""
Phase 0b -- scaffold-clustered paired BCa (bias-corrected and accelerated)
bootstrap, from scratch (rev5 Section 8.4: "Paired BCa bootstrap, 10,000
resamples, clustered by Bemis-Murcko scaffold. Per-compound resampling is
anti-conservative with analogue series present.").

scipy.stats.bootstrap does BCa but not cluster resampling out of the box;
implementing cluster resampling manually here (resample scaffold GROUPS
with replacement, not individual queries) and applying the standard BCa
bias/acceleration correction (Efron & Tibshirani) to the resulting
distribution.
"""
import math
import random
import statistics


def _norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2)))


def _norm_ppf(p):
    """Inverse normal CDF (Acklam's rational approximation) -- avoids a
       scipy.stats.norm dependency for this one call."""
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
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
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > 1 - p_low:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def scaffold_clustered_bca(per_query_values, scaffold_of_query, n_resamples=10000, seed=42, alpha=0.05):
    """per_query_values: list of floats (e.g. paired per-query differences,
       or a single method's per-query metric). scaffold_of_query: parallel
       list of scaffold-group keys, same length.

       Returns {point_estimate, ci_low, ci_high, n_resamples, n_groups,
       n_queries}. Point estimate = pooled mean over all queries."""
    n = len(per_query_values)
    assert n == len(scaffold_of_query)
    groups = {}
    for v, g in zip(per_query_values, scaffold_of_query):
        groups.setdefault(g, []).append(v)
    group_keys = list(groups.keys())
    n_groups = len(group_keys)

    def pooled_mean(selected_group_keys):
        vals = []
        for gk in selected_group_keys:
            vals.extend(groups[gk])
        return statistics.mean(vals) if vals else float("nan")

    observed = pooled_mean(group_keys)

    rng = random.Random(seed)
    boot_stats = []
    for _ in range(n_resamples):
        resample = [group_keys[rng.randrange(n_groups)] for _ in range(n_groups)]
        boot_stats.append(pooled_mean(resample))
    boot_stats.sort()

    # bias-correction z0
    n_less = sum(1 for b in boot_stats if b < observed)
    prop = max(min(n_less / len(boot_stats), 1 - 1e-9), 1e-9)
    z0 = _norm_ppf(prop)

    # acceleration a, via group-level jackknife (leave-one-scaffold-group-out)
    jack_stats = []
    for i in range(n_groups):
        remaining = group_keys[:i] + group_keys[i + 1:]
        jack_stats.append(pooled_mean(remaining))
    jack_mean = statistics.mean(jack_stats)
    num = sum((jack_mean - j) ** 3 for j in jack_stats)
    den = 6.0 * (sum((jack_mean - j) ** 2 for j in jack_stats) ** 1.5)
    a = (num / den) if den != 0 else 0.0

    z_lo = _norm_ppf(alpha / 2)
    z_hi = _norm_ppf(1 - alpha / 2)

    def bca_percentile(z):
        num = z0 + z
        denom = 1 - a * (z0 + z)
        adj_z = z0 + num / denom if denom != 0 else z0
        return _norm_cdf(adj_z)

    p_lo = bca_percentile(z_lo)
    p_hi = bca_percentile(z_hi)
    p_lo = min(max(p_lo, 0.0), 1.0)
    p_hi = min(max(p_hi, 0.0), 1.0)

    idx_lo = max(0, min(len(boot_stats) - 1, int(round(p_lo * (len(boot_stats) - 1)))))
    idx_hi = max(0, min(len(boot_stats) - 1, int(round(p_hi * (len(boot_stats) - 1)))))
    ci_low, ci_high = boot_stats[idx_lo], boot_stats[idx_hi]
    if ci_low > ci_high:
        ci_low, ci_high = ci_high, ci_low

    return {
        "point_estimate": round(observed, 4),
        "ci_low": round(ci_low, 4),
        "ci_high": round(ci_high, 4),
        "ci_excludes_zero": bool(ci_low > 0 or ci_high < 0),
        "n_resamples": n_resamples,
        "n_groups": n_groups,
        "n_queries": n,
    }
