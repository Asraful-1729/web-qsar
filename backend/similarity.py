"""
A3 — natural-product similarity search against the COCONUT-derived index
(see scripts/build_similarity_index.py for how it's built). This is a
SINGLE shared resource (not per-target like models/docking data), fetched
on demand via /api/similarity/download rather than the per-target
download-gate machinery in downloads.py.

Required methods per the upgrade proposal: Morgan/ECFP4 fingerprints +
Tanimoto similarity (the primary ranking), Murcko scaffold matching, and
maximum common substructure (MCS) — computed only for the top handful of
hits since per-pair MCS is comparatively expensive.

Tanimoto is computed via bitwise AND/popcount directly on the PACKED
fingerprint bytes (a lookup table maps each byte 0-255 to its bit count),
not by constructing 738K individual RDKit ExplicitBitVect objects —
that approach measured ~95s just to load the index; this one measures
well under 1s to load AND run the first query (both correctness-checked
against each other: identical Tanimoto scores to 3 decimal places).
"""
import functools
import hashlib
import json
import os
import shutil
import tempfile
import threading

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

INDEX_DIR = os.environ.get("SIMILARITY_INDEX_DIR", "similarity_index")
_morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)

# The same public, zero-credential R2 bucket the per-target model/docking
# downloads already use (see downloads.py) — this is ONE shared resource
# (not per-target), so it gets its own small download-job system instead
# of downloads.py's per-target manifest-driven one.
_DEFAULT_DOWNLOAD_BASE_URL = "https://pub-8ffe9174838b492cab28d50839e49dd7.r2.dev"
DOWNLOAD_BASE_URL = os.environ.get("DOWNLOAD_BASE_URL", _DEFAULT_DOWNLOAD_BASE_URL).rstrip("/")
INDEX_ZIP_URL = f"{DOWNLOAD_BASE_URL}/similarity_index/similarity_index.zip"

_dl_jobs_lock = threading.Lock()
_dl_jobs = {}


class _DownloadCancelled(Exception):
    pass


def available():
    return all(os.path.exists(os.path.join(INDEX_DIR, f)) for f in ("fingerprints.npy", "compounds.csv", "manifest.json"))


def start_download():
    import httpx  # local import matches this project's existing lazy-import convention for optional deps
    import uuid
    job_id = uuid.uuid4().hex
    with _dl_jobs_lock:
        _dl_jobs[job_id] = {"state": "starting", "done": 0, "total": 0, "error": None}

    def run():
        tmp_zip = None
        try:
            os.makedirs(INDEX_DIR, exist_ok=True)
            fd, tmp_zip = tempfile.mkstemp(prefix="similarity_index.", suffix=".zip.part", dir=INDEX_DIR)
            os.close(fd)
            with _dl_jobs_lock:
                _dl_jobs[job_id]["state"] = "downloading"
            with httpx.stream("GET", INDEX_ZIP_URL, timeout=60) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length") or 0)
                with _dl_jobs_lock:
                    _dl_jobs[job_id]["total"] = total
                done = 0
                with open(tmp_zip, "wb") as f:
                    for chunk in r.iter_bytes(chunk_size=1024 * 1024):
                        with _dl_jobs_lock:
                            if _dl_jobs[job_id].get("cancel_requested"):
                                raise _DownloadCancelled()
                        f.write(chunk)
                        done += len(chunk)
                        with _dl_jobs_lock:
                            _dl_jobs[job_id]["done"] = done
            with _dl_jobs_lock:
                _dl_jobs[job_id]["state"] = "extracting"
            with tempfile.TemporaryDirectory(dir=INDEX_DIR) as tmp_extract:
                shutil.unpack_archive(tmp_zip, tmp_extract, format="zip")
                for name in ("fingerprints.npy", "compounds.csv", "manifest.json"):
                    src = os.path.join(tmp_extract, name)
                    if os.path.exists(src):
                        shutil.move(src, os.path.join(INDEX_DIR, name + ".new"))
                for name in ("fingerprints.npy", "compounds.csv", "manifest.json"):
                    new_path = os.path.join(INDEX_DIR, name + ".new")
                    if os.path.exists(new_path):
                        os.replace(new_path, os.path.join(INDEX_DIR, name))
            _load.cache_clear()
            with _dl_jobs_lock:
                _dl_jobs[job_id]["state"] = "done"
        except _DownloadCancelled:
            with _dl_jobs_lock:
                _dl_jobs[job_id]["state"] = "cancelled"
        except Exception as e:
            with _dl_jobs_lock:
                _dl_jobs[job_id]["state"] = "error"; _dl_jobs[job_id]["error"] = str(e)
        finally:
            if tmp_zip and os.path.exists(tmp_zip):
                os.remove(tmp_zip)

    threading.Thread(target=run, daemon=True).start()
    return job_id


def download_progress(job_id):
    with _dl_jobs_lock:
        return dict(_dl_jobs[job_id]) if job_id in _dl_jobs else None


def cancel_download(job_id):
    with _dl_jobs_lock:
        if job_id not in _dl_jobs:
            return False
        _dl_jobs[job_id]["cancel_requested"] = True
        return True


@functools.lru_cache(maxsize=1)
def _load():
    fps_packed = np.load(os.path.join(INDEX_DIR, "fingerprints.npy"))  # (N, 256) uint8
    pop_counts = _POPCOUNT_TABLE[fps_packed].sum(axis=1)  # (N,) — precomputed once, reused every query
    df = pd.read_csv(os.path.join(INDEX_DIR, "compounds.csv"))
    with open(os.path.join(INDEX_DIR, "manifest.json")) as f:
        manifest = json.load(f)
    return fps_packed, pop_counts, df, manifest


def _packed_fingerprint(mol):
    fp = _morgan.GetFingerprint(mol)
    bits = np.zeros(2048, dtype=np.uint8)
    on = list(fp.GetOnBits())
    if on:
        bits[on] = 1
    return np.packbits(bits)


def _murcko_smiles(mol):
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def search(query_smiles, threshold=0.4, top_n=50, mcs_top_n=10):
    """Returns a result dict, or raises ValueError (bad SMILES) /
       FileNotFoundError (index not downloaded yet)."""
    if not available():
        raise FileNotFoundError("similarity index not downloaded yet")
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    query_scaffold = _murcko_smiles(mol)
    q_packed = _packed_fingerprint(mol)
    q_pop = int(_POPCOUNT_TABLE[q_packed].sum())

    fps_packed, pop_counts, df, manifest = _load()
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_packed, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    tanimoto = np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)
    order = np.argsort(-tanimoto)

    results = []
    for idx in order:
        score = float(tanimoto[idx])
        if score < threshold or len(results) >= top_n:
            break
        row = df.iloc[int(idx)]
        murcko = row.get("murcko_framework")
        results.append({
            "id": row["id"], "smiles": row["smiles"], "name": (row.get("name") if pd.notna(row.get("name")) else None),
            "tanimoto": round(score, 3),
            "scaffold_match": bool(query_scaffold is not None and pd.notna(murcko) and murcko == query_scaffold),
            "molecular_weight": (float(row["molecular_weight"]) if pd.notna(row.get("molecular_weight")) else None),
            "chemical_super_class": (row.get("chemical_super_class") if pd.notna(row.get("chemical_super_class")) else None),
            "np_classifier_class": (row.get("np_classifier_class") if pd.notna(row.get("np_classifier_class")) else None),
        })

    # MCS only for the top few — FindMCS is comparatively expensive per
    # pair (up to its own timeout), not worth running on all top_n.
    for r in results[:mcs_top_n]:
        try:
            cand_mol = Chem.MolFromSmiles(r["smiles"])
            mcs_res = rdFMCS.FindMCS([mol, cand_mol], timeout=2)
            r["mcs_smarts"] = mcs_res.smartsString
            r["mcs_n_atoms"] = mcs_res.numAtoms
        except Exception:
            pass

    return {
        "query_scaffold": query_scaffold,
        "n_indexed": len(df),
        "n_results": len(results),
        "manifest": manifest,
        "results": results,
    }
