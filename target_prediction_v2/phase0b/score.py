"""
Phase 0b -- generic re-scoring engine operating on capture.py's saved
per-query neighbour lists. Every weighted-kNN variant (baseline B: plain
10-NN vote; baseline C: similarity+potency-weighted; the nested-CV grid;
the 2x2 potency-weighting factorial) is one call into `score_query()` with
different parameters -- no re-querying the index.

Potency weight function g(pchembl) [Design, per rev5 Section 3]: a smooth,
monotone weight, NOT a hard filter (hard filters cost recall -- [Phase 0]
Section 6, [3]). Implemented as a logistic ramp centered at the potency
threshold, saturating to 1.0 well above it and to a configurable floor well
below it -- this is an explicit, disclosed functional-form choice (rev5
Section 3 flags the functional form as "unvalidated -- tune it"; nested CV
in analyze_h4.py tunes its steepness alongside k/alpha).
"""
import math


def potency_weight(pchembl, threshold=6.0, floor=0.2, steepness=2.0):
    """g(pchembl): logistic ramp, g(threshold)=0.5*(1+floor), saturating to
       1.0 above and `floor` (never zero -- a weight, not a filter) below.
       pchembl=None (no potency annotation for this pair) -> neutral 1.0,
       since absence of a potency value is not evidence of weak potency."""
    if pchembl is None:
        return 1.0
    x = steepness * (pchembl - threshold)
    sig = 1.0 / (1.0 + math.exp(-x))
    return floor + (1.0 - floor) * sig


def score_query(neighbours, k, alpha, potency_fn=None):
    """neighbours: capture.py's ["neighbours"] field --
       [[smiles, tanimoto, [[target_chembl, pchembl_or_null], ...]], ...],
       already sorted by tanimoto descending (capture.py's own sort order).

       Returns {target_chembl: score}, using the review's formula-fix
       (rev3 Section 1 point 5 / rev5 Section 3): each of the k nearest
       neighbour COMPOUNDS votes tanimoto_i^alpha (times potency_fn(pchembl)
       if given) for EVERY target in its own annotation set, not just the
       target being scored.

       alpha=0, potency_fn=None reproduces baseline B (plain unweighted
       10-NN vote: every one of the k neighbours contributes exactly 1.0
       per target, so the resulting integer score IS the "l" in rev5's
       L-score l/10)."""
    scores = {}
    for smi, sim, targets in neighbours[:k]:
        base_w = (sim ** alpha) if alpha != 0 else 1.0
        for tcid, pchembl in targets:
            w = base_w
            if potency_fn is not None:
                w = w * potency_fn(pchembl)
            scores[tcid] = scores.get(tcid, 0.0) + w
    return scores


def rank_from_scores(scores):
    """{target: score} -> [(target, score), ...] sorted descending, ties
       broken by target_chembl id (stable, arbitrary but deterministic --
       recorded so reruns are reproducible per Project Rule 6)."""
    return sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))


def best_rank(ranked, target_set):
    if not target_set:
        return None
    rank_map = {t: i + 1 for i, (t, _) in enumerate(ranked)}
    ranks = [rank_map[t] for t in target_set if t in rank_map]
    return min(ranks) if ranks else None


def popularity_correct(scores, popularity, eps=1.0):
    """Observed/expected vote-mass correction (rev5 Section 3: 'same logic
       as SEA's background model' -- corrects for a target simply having
       many indexed compounds, independent of query-specific similarity).
       `popularity`: {target_chembl: n_indexed_compounds_globally}."""
    return {t: s / (popularity.get(t, 0) + eps) for t, s in scores.items()}
