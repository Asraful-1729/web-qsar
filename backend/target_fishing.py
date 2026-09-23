"""
A2 — "compound -> target" mode: given a compound, which protein targets
is it likely to hit? Complements the disease -> target flow TargetBrowser
already offers (see docking/recommend.py) — researchers enter from either
direction depending on the project.

Method: ligand-based target prediction via similarity to KNOWN bioactive
compounds — the "guilt by association" principle SwissTargetPrediction and
similar tools use, NOT a trained multi-label classifier (none exists in
this app). Results are "similar known actives found for these targets,"
never a calibrated probability, and are presented to the caller that way.

Reference pool: a broad ChEMBL bioactivity pull (human targets, IC50/Ki/
Kd/EC50, assay_confidence_score>=8, pchembl_value present — see
scripts/fetch_chembl_activities.py / build_target_fishing_index.py for the
exact filter and provenance, recorded in the index's own manifest.json),
independent of this app's own ~64 QSAR-modeled targets. An earlier version
of this module reused models/curated/*.csv (this app's own QSAR training
data, repurposed) as the search pool — that meant a query could only ever
be routed back to one of those same 64 targets, no matter how it actually
scored against the broader universe of real ChEMBL bioactivity; fixed by
switching the whole reference pool to the broader pull above. Most targets
in this pool will have no corresponding target_id in this app's own
docking_registry.json — that's expected (see target_id=None handling
below), not a bug: the app doesn't have a docking/QSAR setup for every
protein ChEMBL has bioactivity data for, and target fishing's job is to
surface candidates from the WHOLE space, not just the ones already wired
into this app's other tabs.

Reuses the same packed-bit-popcount Tanimoto approach as similarity.py
(measured there: >100x faster than constructing individual RDKit
bitvector objects per compound) — fingerprints are precomputed OFFLINE
into backend/target_fishing_index/ (see build_target_fishing_index.py),
never fingerprinted live per query.

Evidence scoring (search()'s evidence_score, per target): a transparent,
documented, RULE-BASED combination of four signals — never call it or
present it as a probability, it has not been statistically calibrated or
validated against a held-out benchmark:

    evidence_score = 0.6 * best_similarity
                    + 0.2 * scaffold_bonus   (0..1: min(n_scaffolds, 5) / 5)
                    + 0.2 * potency_bonus    (0..1: (best_pchembl - 5) / 5, clamped)

  - best_similarity dominates (0.6 weight) because it's the primary,
    directly-measured structural evidence linking the query to a known
    active — everything else is corroborating context, not a substitute.
  - scaffold_bonus rewards CONVERGENT evidence: several independent
    chemotypes hitting the same target is much stronger support than many
    near-identical analogs of one scaffold (which is really one data
    point, weighted to look like several) — saturates at 5 distinct
    Murcko scaffolds so one promiscuous scaffold family can't dominate.
  - potency_bonus rewards evidence coming from a genuinely potent known
    active (pChEMBL 10 = 0.1 nM) over a barely-active one (pChEMBL 5 =
    10 uM, this pipeline's own activity floor) — zero contribution when a
    target has no potency-annotated support at all.
  - Weights (0.6/0.2/0.2) are a documented, deliberate judgment call, NOT
    fit to data — there is no benchmark yet to fit them against (see the
    module-level TODO: build a scaffold-split held-out validation set and
    report Top-K recovery / MRR / enrichment before trusting this ranking
    beyond "a reasonable starting heuristic").
"""
import functools
import json
import os

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

INDEX_DIR = os.environ.get("TARGET_FISHING_INDEX_DIR", os.path.join(os.path.dirname(__file__), "target_fishing_index"))
_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

MIN_PCHEMBL_FLOOR = 5.0    # matches the index build's own activity floor (10 uM) — potency_bonus's zero point
MAX_PCHEMBL_CAP = 10.0     # 0.1 nM — potency_bonus's saturation point
MAX_SCAFFOLDS_CAP = 5      # distinct Murcko scaffolds at which scaffold_bonus saturates


def available():
    return os.path.exists(os.path.join(INDEX_DIR, "fingerprints.npz")) and os.path.exists(
        os.path.join(INDEX_DIR, "compounds.csv.gz")
    )


def _packed_fingerprint(mol):
    fp = _morgan.GetFingerprint(mol)
    bits = np.zeros(2048, dtype=np.uint8)
    on = list(fp.GetOnBits())
    if on:
        bits[on] = 1
    return np.packbits(bits)


@functools.lru_cache(maxsize=1)
def _load():
    """Loads the offline-precomputed index (see module docstring) — no
       fingerprinting happens here, just an npz/gzip-csv read, so this is
       fast (well under a second) rather than the many minutes it took to
       fingerprint the pool at build time."""
    fps_arr = np.load(os.path.join(INDEX_DIR, "fingerprints.npz"))["fps"]
    df = pd.read_csv(os.path.join(INDEX_DIR, "compounds.csv.gz"))
    pop_counts = _POPCOUNT_TABLE[fps_arr].sum(axis=1)
    return fps_arr, pop_counts, df


def _evidence_score(best_sim, n_scaffolds, best_pchembl):
    scaffold_bonus = min(n_scaffolds, MAX_SCAFFOLDS_CAP) / MAX_SCAFFOLDS_CAP
    if best_pchembl is None:
        potency_bonus = 0.0
    else:
        potency_bonus = (best_pchembl - MIN_PCHEMBL_FLOOR) / (MAX_PCHEMBL_CAP - MIN_PCHEMBL_FLOOR)
        potency_bonus = max(0.0, min(1.0, potency_bonus))
    return round(0.6 * best_sim + 0.2 * scaffold_bonus + 0.2 * potency_bonus, 3)


def _search_against(fps_arr, pop_counts, df, q_packed, q_pop, threshold, max_compounds_per_target):
    """Core vectorized Tanimoto + per-target aggregation, factored out of
       search() so it can run against ARBITRARY (fps_arr, pop_counts, df)
       arrays, not just the production index _load() returns. The only
       other caller is scripts/target_fishing_benchmark.py, which runs
       this same function against a held-out-scaffold-REDUCED copy of the
       index (a compound and every compound sharing its scaffold removed)
       — reusing this instead of a separate reimplementation means the
       benchmark measures the exact ranking logic production users get,
       not a hand-rolled approximation that could silently drift out of
       sync with real changes here.

       Aggregation is deliberately all vectorized pandas groupby/agg, NOT
       a Python-level `for _, r in hits.iterrows(): ...` loop over hit rows
       (an earlier version did exactly that) — iterrows() builds a fresh
       pandas Series per row, fine at the few-thousand-hit scale a normal
       threshold (e.g. 0.4) produces, but at threshold=0.0 (which the
       benchmark needs, to check whether the true target appears ANYWHERE
       in the ranking) EVERY indexed row passes, and iterating that many
       Series objects in Python measured ~20s for a single query against
       the full 1.3M-row index — a real production footgun too (any user
       setting a low similarity threshold hit this same stall). The
       groupby/agg + a bounded `.head(N)` per group below do the
       equivalent aggregation in C, without materializing a Python object
       per hit row.

       Returns the `results` list only (sorted by evidence_score
       descending) — search() wraps this with the pool-level summary
       fields; the benchmark doesn't need those per query."""
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    tanimoto = np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)

    mask = tanimoto >= threshold
    hits = df.loc[mask].copy()
    if hits.empty:
        return []
    hits["tanimoto"] = tanimoto[mask]
    # Sorted ONCE, globally, by similarity descending — everything below
    # relies on this order: groupby(...).head(N) takes each target's first
    # N rows AS ENCOUNTERED, which (because of this sort) are exactly its
    # N most-similar supporting compounds, with no separate per-group sort.
    hits = hits.sort_values("tanimoto", ascending=False)

    grouped = hits.groupby("target_chembl", sort=False)
    agg = grouped.agg(
        target_pref_name=("target_pref_name", "first"),
        target_id=("target_id", "first"),
        n_similar_actives=("tanimoto", "size"),
        best_similarity=("tanimoto", "max"),
        mean_similarity=("tanimoto", "mean"),
        best_pchembl=("pchembl_value", "max"),
        mean_pchembl=("pchembl_value", "mean"),
        n_scaffolds=("murcko_scaffold", "nunique"),
    )

    top_compounds = hits.groupby("target_chembl", sort=False).head(max_compounds_per_target)
    compounds_by_target = {}
    for row in top_compounds.itertuples(index=False):
        compounds_by_target.setdefault(row.target_chembl, []).append({
            "smiles": row.smiles, "tanimoto": round(float(row.tanimoto), 3),
            "pchembl_value": (float(row.pchembl_value) if pd.notna(row.pchembl_value) else None),
        })

    results = []
    for target_chembl, row in agg.iterrows():
        best_sim = float(row["best_similarity"])
        best_pchembl = row["best_pchembl"]
        best_pchembl = float(best_pchembl) if pd.notna(best_pchembl) else None
        mean_pchembl = row["mean_pchembl"]
        n_scaffolds = int(row["n_scaffolds"]) if pd.notna(row["n_scaffolds"]) else 0
        results.append({
            "target_chembl": target_chembl,
            "target_pref_name": (row["target_pref_name"] if pd.notna(row["target_pref_name"]) else None),
            "target_id": (row["target_id"] if pd.notna(row["target_id"]) else None),
            "n_similar_actives": int(row["n_similar_actives"]),
            "best_similarity": round(best_sim, 3),
            "mean_similarity": round(float(row["mean_similarity"]), 3),
            "best_pchembl": (round(best_pchembl, 2) if best_pchembl is not None else None),
            "mean_pchembl": (round(float(mean_pchembl), 2) if pd.notna(mean_pchembl) else None),
            "n_scaffolds": n_scaffolds,
            "evidence_score": _evidence_score(best_sim, n_scaffolds or 1, best_pchembl),
            "compounds": compounds_by_target.get(target_chembl, []),
        })

    # Ranked by best_similarity, NOT evidence_score — the scaffold-split
    # validation benchmark (see TARGET_FISHING_BENCHMARK.md, 2026-09-20;
    # 7,295 novel-scaffold evaluation instances) found evidence_score
    # consistently UNDERPERFORMS raw similarity for actual target recovery
    # (Top-10: 91.1% vs 92.2%; MRR: 0.693 vs 0.705 — small but consistent
    # across every Top-K cut and both the primary and secondary buckets,
    # not noise). evidence_score is still computed and returned on every
    # result (a real, differently-purposed corroboration/robustness
    # signal), just no longer used to decide default ordering — do not
    # switch this back without a benchmark result showing it actually
    # improves recovery.
    results.sort(key=lambda r: -r["best_similarity"])
    return results


def search(query_smiles, threshold=0.4, max_compounds_per_target=5):
    """Returns a result dict, or raises ValueError (bad SMILES) /
       FileNotFoundError (no index built yet). See _search_against() for
       the actual similarity/aggregation logic."""
    if not available():
        raise FileNotFoundError("target-fishing index not built")
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    q_packed = _packed_fingerprint(mol)
    q_pop = int(_POPCOUNT_TABLE[q_packed].sum())

    fps_arr, pop_counts, df = _load()
    results = _search_against(fps_arr, pop_counts, df, q_packed, q_pop, threshold, max_compounds_per_target)

    return {
        "n_indexed_compound_target_pairs": len(df),
        "n_targets_searched": int(df["target_chembl"].nunique()),
        "n_targets_matched": len(results),
        "results": results,
    }


def suggest_compounds(query, limit=8):
    """As-you-type suggestions from the SAME indexed pool search() ranks
       against — no fingerprinting, just a plain substring match, so this
       is cheap enough to call on every keystroke. Matches either a SMILES
       fragment (e.g. a ring system just pasted in) or a target's ChEMBL
       id / name against the pool, deduplicated by SMILES so a compound
       tested against many targets only shows once (with its best-potency
       target as context).

       Returns a list of {smiles, target_chembl, target_pref_name,
       target_id, pchembl_value}, or [] for an empty/too-short query —
       never raises for a bad query, only FileNotFoundError if the index
       itself is missing."""
    if not available():
        raise FileNotFoundError("target-fishing index not built")
    query = (query or "").strip()
    if len(query) < 2:
        return []
    limit = max(1, min(int(limit), 20))

    _, _, df = _load()
    q = query.lower()
    mask = (
        df["smiles"].str.lower().str.contains(q, regex=False)
        | df["target_chembl"].str.lower().str.contains(q, regex=False)
        | df["target_pref_name"].str.lower().str.contains(q, regex=False, na=False)
    )
    hits = df[mask]
    if hits.empty:
        return []
    # one row per compound: keep its highest-pchembl (most-potent, most
    # informative) target when the same SMILES appears against several.
    hits = hits.sort_values("pchembl_value", ascending=False).drop_duplicates(subset="smiles").head(limit)

    out = []
    for _, r in hits.iterrows():
        pchembl = r.get("pchembl_value")
        out.append({
            "smiles": r["smiles"],
            "target_chembl": r["target_chembl"],
            "target_pref_name": (r.get("target_pref_name") if pd.notna(r.get("target_pref_name")) else None),
            "target_id": (r.get("target_id") if pd.notna(r.get("target_id")) else None),
            "pchembl_value": (float(pchembl) if pd.notna(pchembl) else None),
        })
    return out
