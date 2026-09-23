"""
Per-compound comparison against a target's known active/decoy distribution.
[orchestration WRITTEN; docking runs UNVALIDATED — same status as the rest
of this package's subprocess-backed code]

fresh_decoy_validation(): EXPENSIVE (~50 extra Vina runs). For one compound
worth deeper scrutiny (a shortlisted hit), generates NEW decoys
property-matched + topologically-dissimilar to THAT SPECIFIC compound
(reusing scripts/generate_decoys.py's exact DUD-E-style selection against a
single-molecule query instead of a target's curated actives) and docks all
of them with the same receptor/box/engine, so the percentile reflects this
molecule's own chemical neighborhood.
"""
import os


def percentile_rank(score, reference_scores):
    """% of reference_scores this score beats. Lower Vina score = better
       binding, so 'beats' means strictly lower; ties count as half-beaten
       (standard percentile-of-score convention)."""
    if not reference_scores:
        return None
    worse = sum(1 for s in reference_scores if s > score)
    tied = sum(1 for s in reference_scores if s == score)
    return round(100.0 * (worse + 0.5 * tied) / len(reference_scores), 1)


def discrimination_label(percentile):
    if percentile is None:
        return None
    if percentile >= 90:
        return "Strong"
    if percentile >= 65:
        return "Moderate"
    return "Weak"


def fresh_decoy_validation(target_id, smiles, profile, engine=None, n_decoys=50, seed=None, progress_cb=None):
    """~(n_decoys + 1) Vina runs. Returns a dict with compound_score,
       percentile, discrimination label, and the per-decoy scores — or
       {'error': ...} if the compound/pool can't support the test.

       Decoys are docked CONCURRENTLY (a thread pool, one Vina subprocess
       per worker) — sequentially this was ~50 back-to-back Vina calls,
       10+ minutes wall-clock with 31 of 32 cores sitting idle the whole
       time, since a single exhaustiveness=8 Vina run only keeps ~8 cores
       busy. Workers get an explicit --cpu cap (docking/engines.py) so N
       concurrent Vina processes divide the machine's cores instead of each
       one independently grabbing all of them and thrashing.

       progress_cb(done, total), if given, is called after every dock
       (compound first, then each decoy as its own thread finishes) so a
       caller can show live progress instead of one static "generating..."
       message for the whole run — the single biggest reason this used to
       look hung: no visible movement for many minutes."""
    import os
    import pandas as pd
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from .pipeline import dock_compound
    from .engines import VinaEngine, NullRescorer
    from scripts.generate_decoys import build_pool, select_decoys

    engine = engine or VinaEngine()
    pool = build_pool(exclude_target_id=target_id)
    actives_df = pd.DataFrame([{"smiles": smiles}])
    decoys = select_decoys(actives_df, pool, n_per_active=n_decoys, seed=seed if seed is not None else 42)
    if not decoys:
        return {"error": "No property-matched, topologically-dissimilar decoys found for this "
                         "compound in the candidate pool (models/curated/*.csv) — it may be "
                         "structurally unusual relative to every other target's curated compounds."}

    total = len(decoys) + 1
    done = 0
    if progress_cb:
        progress_cb(done, total)

    compound_res = dock_compound(profile, smiles, engine=engine, rescorer=NullRescorer())
    compound_pose = compound_res.get("consensus_pose")
    compound_score = compound_pose["score"] if compound_pose else None
    done += 1
    if progress_cb:
        progress_cb(done, total)
    if compound_score is None:
        return {"error": "This compound did not produce a valid (PoseBusters-passing) pose against "
                         "this receptor — nothing to rank against decoys."}

    n_workers = min(8, max(1, os.cpu_count() or 4))
    cpu_per_worker = max(1, (os.cpu_count() or n_workers) // n_workers)

    def _dock_one(d):
        worker_engine = VinaEngine(binary=engine.binary, exhaustiveness=engine.exhaustiveness, cpu=cpu_per_worker)
        r = dock_compound(profile, d["smiles"], engine=worker_engine, rescorer=NullRescorer())
        pose = r.get("consensus_pose")
        return {"name": d["name"], "smiles": d["smiles"], "source_target": d.get("source_target"),
               "score": pose["score"] if pose else None}

    decoy_rows = []
    with ThreadPoolExecutor(max_workers=n_workers) as pool_exec:
        futures = {pool_exec.submit(_dock_one, d): d for d in decoys}
        for fut in as_completed(futures):
            decoy_rows.append(fut.result())
            done += 1
            if progress_cb:
                progress_cb(done, total)

    valid_scores = [d["score"] for d in decoy_rows if d["score"] is not None]
    pct = percentile_rank(compound_score, valid_scores) if valid_scores else None
    import numpy as np
    decoy_stats = None
    if valid_scores:
        arr = np.array(valid_scores, dtype=float)
        decoy_stats = {"mean": round(float(arr.mean()), 3), "median": round(float(np.median(arr)), 3),
                       "sd": round(float(arr.std(ddof=1)), 3) if len(arr) > 1 else 0.0,
                       "min": round(float(arr.min()), 3), "max": round(float(arr.max()), 3)}
    return {
        "compound_score": compound_score,
        "n_decoys_generated": len(decoys), "n_decoys_docked": len(valid_scores),
        "n_decoys_failed": len(decoy_rows) - len(valid_scores),
        "percentile": pct, "discrimination": discrimination_label(pct),
        "decoy_stats": decoy_stats,
        "decoys": decoy_rows,
        # B7's "show the run settings used for THIS validation" — same
        # transparency B6 applies to a live docking submission, applied
        # here too so the percentile above can actually be interpreted.
        "run_settings": {
            "center": profile.get("center"), "box_size": profile.get("box_size"),
            "docking_mode": "blind" if profile.get("site_source") == "blind_whole_protein" else "site_specific",
            "pdb_source": profile.get("pdb_source"), "exhaustiveness": engine.exhaustiveness,
        },
    }
