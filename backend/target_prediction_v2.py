"""
Target Prediction v2 -- the live serving module for the validated v2
method. See target_prediction_v2/METHODS_AND_VALIDATION.md (repo root)
for the full account of what this method is, what was tested and
rejected, and its validated performance; this docstring only covers what
is specific to SERVING it (research-handoff roadmap item #5: "SMILES in,
ranked targets + evidence + confidence out").

Method (validated -- G2-passed against a genuinely held-out, never-
tuned-on locked scaffold test; PHASE3_G2_LOCKED_TEST_RESULT.md):
density-adaptive retrieval.
  - Count neighbour COMPOUNDS in v2_index at Tanimoto>=0.5 -- "density".
  - density>=8: pool votes from the k=10 nearest-neighbour compounds,
    each voting for every target it is annotated against, weighted by a
    potency logistic ramp (threshold=6.0, floor=0.0, steepness=5.0),
    PLUS an additive orthologue term (non-human-species evidence,
    Tanimoto>=0.4 against a separate orthologue index -- adopted because
    it specifically improves primary/mechanism-tier recovery).
  - density<8: fall back to v1's original best-similarity ranking
    (single nearest-neighbour compound per target, computed over the
    FULL index, not the pooling regime's capped neighbour list) --
    re-fit and re-validated against the rebuilt v2 index, not carried
    over unchanged from v1.

NOT shipped here (tested and rejected/deferred -- see
METHODS_AND_VALIDATION.md Section 4 for the full disclosure):
  - Popularity correction: significantly HARMFUL at every strength tested.
  - Density-stratified calibration: failed its own pre-registered Brier/
    ECE gate against the raw L-score -- so the raw L-score ships as the
    confidence signal, not a "calibrated" probability.
  - The logistic-regression stacker: a real but small gain (PHASE3A_STACKER.md),
    never trained on the full dataset and serialized for production use,
    so it is not called here.

Confidence signal: the raw L-score (pooled regime -- a potency+orthologue-
weighted vote count, roughly 0-10+) or best_similarity (fallback regime).
NEITHER is a calibrated probability -- say so in every result, since Phase 4
measured calibration and it failed to beat the raw score.

CLI usage: python3 target_prediction_v2.py "<SMILES>" [top_k]
"""
import functools
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
V2_INDEX_DIR = os.environ.get(
    "TARGET_PREDICTION_V2_INDEX_DIR",
    os.path.join(HERE, "..", "target_prediction_v2", "phase2", "v2_index"),
)
ORTHO_INDEX_DIR = os.environ.get(
    "TARGET_PREDICTION_V2_ORTHO_DIR",
    os.path.join(HERE, "..", "target_prediction_v2", "phase3", "orthologue_index"),
)
TARGET_NAME_MAP_PATH = os.environ.get(
    "TARGET_PREDICTION_V2_NAME_MAP",
    os.path.join(HERE, "..", "target_prediction_v2", "phase2", "data", "target_name_map.json"),
)

_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

DENSITY_THRESHOLD = 8       # phase3/PHASE3A_DENSITY_REFIT.md -- re-fit against the rebuilt v2 index
DENSITY_TANIMOTO_CUTOFF = 0.5
POOL_K = 10
NEIGHBOR_CAP = 300          # matches every Phase 3 evaluation script -- pooling only ever uses k=10 of these
ORTHO_TANIMOTO_MIN = 0.4
POTENCY_CFG = dict(threshold=6.0, floor=0.0, steepness=5.0)


def available():
    return (
        os.path.exists(os.path.join(V2_INDEX_DIR, "fingerprints.npz"))
        and os.path.exists(os.path.join(V2_INDEX_DIR, "compounds.csv.gz"))
    )


def _packed_fingerprint(mol):
    fp = _morgan.GetFingerprint(mol)
    bits = np.zeros(2048, dtype=np.uint8)
    on = list(fp.GetOnBits())
    if on:
        bits[on] = 1
    return np.packbits(bits)


def _tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop):
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)


def _potency_weight(pchembl, threshold=6.0, floor=0.0, steepness=5.0):
    """g(pchembl): logistic ramp -- absence of a potency value (None) is
       treated as neutral (1.0), not as evidence of weak potency."""
    if pchembl is None:
        return 1.0
    x = steepness * (pchembl - threshold)
    sig = 1.0 / (1.0 + math.exp(-x))
    return floor + (1.0 - floor) * sig


@functools.lru_cache(maxsize=1)
def _load():
    fps = np.load(os.path.join(V2_INDEX_DIR, "fingerprints.npz"))["fps"]
    df = pd.read_csv(os.path.join(V2_INDEX_DIR, "compounds.csv.gz"))
    pop = _POPCOUNT_TABLE[fps].sum(axis=1)
    return fps, pop, df


@functools.lru_cache(maxsize=1)
def _load_ortho():
    if not os.path.exists(os.path.join(ORTHO_INDEX_DIR, "fingerprints.npz")):
        return None
    fps = np.load(os.path.join(ORTHO_INDEX_DIR, "fingerprints.npz"))["fps"]
    df = pd.read_csv(os.path.join(ORTHO_INDEX_DIR, "compounds.csv.gz"))
    pop = _POPCOUNT_TABLE[fps].sum(axis=1)
    return fps, pop, df


@functools.lru_cache(maxsize=1)
def _load_target_names():
    if not os.path.exists(TARGET_NAME_MAP_PATH):
        return {}
    with open(TARGET_NAME_MAP_PATH) as f:
        return json.load(f)


def _target_name(target_chembl):
    return _load_target_names().get(target_chembl)


def _query_tanimoto_df(fps, pop, df, q_packed, q_pop):
    """One vectorized Tanimoto pass against the FULL index -- both the
       pooling-regime neighbour list and the fallback-regime full
       best-similarity ranking are derived from this SAME dataframe, so
       the query is only ever scored against the index once."""
    tanimoto = _tanimoto_to_all(q_packed, fps, pop, q_pop)
    tmp = df.copy()
    tmp["tanimoto"] = tanimoto
    return tmp


def _density_and_neighbours(tmp, cap=NEIGHBOR_CAP):
    """density: distinct compounds at Tanimoto>=DENSITY_TANIMOTO_CUTOFF,
       computed over the FULL index (not capped) -- matching how density
       was computed in every Phase 3 fit/validation script. neighbours:
       top-`cap` distinct compounds by similarity, the pooling input."""
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    density = int((by_compound["tanimoto"] >= DENSITY_TANIMOTO_CUTOFF).sum())
    top_n = by_compound.sort_values("tanimoto", ascending=False).head(cap)
    if top_n.empty:
        return density, []
    neighbor_rows = tmp[tmp["smiles"].isin(set(top_n.index))][["smiles", "target_chembl", "pchembl_value"]]
    by_smi = {}
    for r in neighbor_rows.itertuples(index=False):
        pv = None if r.pchembl_value != r.pchembl_value else float(r.pchembl_value)
        by_smi.setdefault(r.smiles, []).append((r.target_chembl, pv))
    neighbours = [(smi, float(sim), by_smi.get(smi, [])) for smi, sim in top_n["tanimoto"].items()]
    return density, neighbours


def _pooled_scores_with_evidence(neighbours, k=POOL_K):
    """Returns ({target: score}, {target: evidence-list}) -- the evidence
       list records exactly which neighbour compounds/potency values
       contributed, for per-prediction interpretability (research-handoff
       roadmap #6), not just an aggregate number."""
    scores, evidence = {}, {}
    for smi, sim, targets in neighbours[:k]:
        for tcid, pchembl in targets:
            w = _potency_weight(pchembl, **POTENCY_CFG)
            scores[tcid] = scores.get(tcid, 0.0) + w
            evidence.setdefault(tcid, []).append({
                "smiles": smi, "tanimoto": round(sim, 4),
                "pchembl_value": pchembl, "potency_weight": round(w, 4),
            })
    return scores, evidence


def _ortho_scores_with_evidence(q_packed, q_pop):
    ortho = _load_ortho()
    if ortho is None:
        return {}, {}
    ortho_fps, ortho_pop, ortho_df = ortho
    tanimoto = _tanimoto_to_all(q_packed, ortho_fps, ortho_pop, q_pop)
    tmp = ortho_df.copy()
    tmp["tanimoto"] = tanimoto
    tmp = tmp[tmp["tanimoto"] >= ORTHO_TANIMOTO_MIN]
    if tmp.empty:
        return {}, {}
    by_compound = tmp.groupby("smiles").agg(tanimoto=("tanimoto", "max"))
    top_n = by_compound.sort_values("tanimoto", ascending=False).head(NEIGHBOR_CAP)
    scores, evidence = {}, {}
    has_species = "species_provenance" in tmp.columns
    for nsmi in top_n.index:
        rows = tmp[tmp["smiles"] == nsmi]
        sim = float(top_n.loc[nsmi, "tanimoto"])
        species_col = rows["species_provenance"] if has_species else [None] * len(rows)
        for t, species in zip(rows["target_chembl"], species_col):
            scores[t] = scores.get(t, 0.0) + 1.0
            evidence.setdefault(t, []).append({
                "smiles": nsmi, "tanimoto": round(sim, 4), "species_provenance": species,
            })
    return scores, evidence


def _best_similarity_ranking(tmp):
    """v1's original method, re-validated against the rebuilt v2 index --
       single nearest-neighbour compound per target, computed over the
       FULL (uncapped) index, matching backend/target_fishing.py's
       _search_against(threshold=0.0) logic."""
    ordered = tmp.sort_values("tanimoto", ascending=False)
    agg = ordered.groupby("target_chembl", sort=False).agg(
        best_similarity=("tanimoto", "max"),
        best_compound_smiles=("smiles", "first"),
        best_compound_pchembl=("pchembl_value", "first"),
    )
    return agg


def predict(query_smiles, top_k=25):
    """Returns a result dict, or raises ValueError (bad SMILES) /
       FileNotFoundError (index not built).

       {
         "regime": "pooled" | "best_similarity_fallback",
         "density": int, "density_threshold": int,
         "n_neighbours_considered": int, "n_targets_indexed": int,
         "results": [
            {"rank", "target_chembl", "target_pref_name", "score",
             "confidence_label", "evidence": {...}}, ...
         ],
       }

       evidence in the pooled regime: {"native_neighbours": [...],
       "orthologue_neighbours": [...]} -- each entry a supporting compound
       with its similarity/potency/(species, for orthologue). evidence in
       the fallback regime: {"best_compound_smiles", "best_compound_pchembl"}.
    """
    if not available():
        raise FileNotFoundError("target_prediction_v2 index not built")
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    q_packed = _packed_fingerprint(mol)
    q_pop = int(_POPCOUNT_TABLE[q_packed].sum())

    fps, pop, df = _load()
    tmp = _query_tanimoto_df(fps, pop, df, q_packed, q_pop)
    density, neighbours = _density_and_neighbours(tmp)

    results = []
    if density >= DENSITY_THRESHOLD:
        regime = "pooled"
        scores, native_evidence = _pooled_scores_with_evidence(neighbours, k=POOL_K)
        ortho_scores, ortho_evidence = _ortho_scores_with_evidence(q_packed, q_pop)
        for t, s in ortho_scores.items():
            scores[t] = scores.get(t, 0.0) + s
        ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
        for i, (tcid, score) in enumerate(ranked[:top_k]):
            results.append({
                "rank": i + 1,
                "target_chembl": tcid,
                "target_pref_name": _target_name(tcid),
                "score": round(score, 4),
                "confidence_label": ("L-score (uncalibrated pooled-vote score; potency- and "
                                      "orthologue-weighted; NOT a probability)"),
                "evidence": {
                    "native_neighbours": sorted(native_evidence.get(tcid, []), key=lambda e: -e["tanimoto"]),
                    "orthologue_neighbours": sorted(ortho_evidence.get(tcid, []), key=lambda e: -e["tanimoto"]),
                },
            })
    else:
        regime = "best_similarity_fallback"
        agg = _best_similarity_ranking(tmp)
        ranked = agg.sort_values("best_similarity", ascending=False)
        for i, (tcid, row) in enumerate(ranked.head(top_k).iterrows()):
            bp = row["best_compound_pchembl"]
            results.append({
                "rank": i + 1,
                "target_chembl": tcid,
                "target_pref_name": _target_name(tcid),
                "score": round(float(row["best_similarity"]), 4),
                "confidence_label": (f"best-similarity (uncalibrated; single nearest-neighbour match -- "
                                     f"density-adaptive gate fell back here because fewer than "
                                     f"{DENSITY_THRESHOLD} neighbours were found at Tanimoto>="
                                     f"{DENSITY_TANIMOTO_CUTOFF})"),
                "evidence": {
                    "best_compound_smiles": row["best_compound_smiles"],
                    "best_compound_pchembl": (float(bp) if pd.notna(bp) else None),
                },
            })

    return {
        "regime": regime,
        "density": density,
        "density_threshold": DENSITY_THRESHOLD,
        "n_neighbours_considered": len(neighbours),
        "n_targets_indexed": int(df["target_chembl"].nunique()),
        "results": results,
    }


def suggest_compounds(query, limit=8):
    """As-you-type suggestions, matched against this SAME v2 index --
       ported from the removed v1 engine's identically-shaped function
       when v1 was deleted (query substring against a compound's SMILES,
       its target's ChEMBL id, or its target's preferred name; no
       fingerprinting, cheap enough for every keystroke). v2's own index
       (see _load()) carries no target_pref_name column of its own --
       names come from the separate _load_target_names() map -- and no
       app-registry target_id mapping at all (unlike v1's index, which
       had one baked in for the subset of targets this app also docks/
       models); target_id is therefore always None here, same as any
       ChEMBL target v1 itself had no app mapping for.

       Returns a list of {smiles, target_chembl, target_pref_name,
       target_id, pchembl_value}, or [] for an empty/too-short query --
       never raises for a bad query, only FileNotFoundError if the index
       itself is missing."""
    if not available():
        raise FileNotFoundError("target_prediction_v2 index not built")
    query = (query or "").strip()
    if len(query) < 2:
        return []
    limit = max(1, min(int(limit), 20))

    _, _, df = _load()
    names = _load_target_names()
    q = query.lower()
    name_matches = {tc for tc, nm in names.items() if nm and q in nm.lower()}
    mask = (
        df["smiles"].str.lower().str.contains(q, regex=False)
        | df["target_chembl"].str.lower().str.contains(q, regex=False)
        | df["target_chembl"].isin(name_matches)
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
        tc = r["target_chembl"]
        out.append({
            "smiles": r["smiles"],
            "target_chembl": tc,
            "target_pref_name": names.get(tc),
            "target_id": None,
            "pchembl_value": (float(pchembl) if pd.notna(pchembl) else None),
        })
    return out


def _cli():
    if len(sys.argv) < 2:
        print("Usage: python3 target_prediction_v2.py \"<SMILES>\" [top_k]", file=sys.stderr)
        sys.exit(1)
    smiles = sys.argv[1]
    top_k = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    try:
        result = predict(smiles, top_k=top_k)
    except (ValueError, FileNotFoundError) as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    _cli()
