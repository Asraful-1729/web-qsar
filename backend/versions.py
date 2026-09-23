"""
Software-version snapshot for reproducibility metadata (A6) — the
"Compound database version / Software version / ... " fields a run's
export package's metadata.json carries. Cached after the first call since
none of this changes while the app is running, and the Vina subprocess
call has real (if small) latency.
"""
import shutil
import subprocess

_cache = None


def snapshot():
    global _cache
    if _cache is not None:
        return _cache
    out = {}
    try:
        import rdkit
        out["rdkit"] = rdkit.__version__
    except Exception:
        pass
    try:
        import meeko
        out["meeko"] = getattr(meeko, "__version__", None)
    except Exception:
        pass
    try:
        import admet_ai
        out["admet_ai"] = getattr(admet_ai, "__version__", None)
    except Exception:
        pass
    vina_bin = shutil.which("vina")
    if vina_bin:
        try:
            r = subprocess.run([vina_bin, "--version"], capture_output=True, text=True, timeout=5)
            out["vina"] = (r.stdout or r.stderr).strip()
        except Exception:
            pass
    _cache = out
    return out
