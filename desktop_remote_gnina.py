"""
Desktop-build-only patch for GNINA rescoring — used ONLY by desktop.py,
never imported by the plain web deployment (backend/app.py stays exactly
as it is, unmodified, when run standalone via `uvicorn app:app`).

Why this exists: GNINA has no official Windows build at all (checked
every GitHub release, v1.0 through v1.3.3 — every asset is a Linux ELF
binary; qsar-desktop's own prior desktop build effort independently
found and documented the same thing). So a Windows desktop build can
never run GNINA locally. Instead, when no local `gnina` binary is found
on PATH, this patches docking.engines.GninaRescorer to call out to a
remote GNINA service (backend/remote_gnina_service.py) running on a real
GPU machine, optionally — if that machine isn't reachable (offline, no
network), rescoring is simply unavailable and docking still works fine
with Vina alone (GNINA has always been an optional "second opinion",
never required for docking to activate).

HOW THE PATCH IS APPLIED (see desktop.py): this module only DEFINES the
replacement behaviour — desktop.py assigns these onto the actual
docking.engines.GninaRescorer CLASS OBJECT at startup:

    engines.GninaRescorer.available = patched_available
    engines.GninaRescorer.rescore = make_patched_rescore(original_rescore)

Patching the class object's own attributes (not rebinding a module-level
name) means every existing call site (docking/pipeline.py's `from .engines
import GninaRescorer` at its own module level, app.py's lazy per-request
import) picks this up automatically, regardless of import order — they
all hold a reference to the SAME class object, and mutating its methods
is visible through every one of those references. Zero lines of any
existing file change.

Local gnina (if ever present, e.g. someone runs this desktop build on
Linux) is always preferred and used completely unchanged via
`original_rescore` — this module only changes behaviour in the "no local
binary" case, which is the expected, common case on Windows.

Config (env vars, all optional):
  REMOTE_GNINA_URL     base URL of the remote service (e.g. http://host:8300)
  REMOTE_GNINA_TOKEN    bearer token matching that server's REMOTE_GNINA_TOKEN
  REMOTE_GNINA_TIMEOUT  seconds, default 45 (a real docking rescoring call,
                        not instant -- generous but bounded)

If REMOTE_GNINA_URL is unset, remote rescoring is simply disabled (not an
error) -- available() returns whatever the local shutil.which() check
says, matching the ORIGINAL class's own behaviour exactly.
"""
import json
import os
import shutil
import urllib.error
import urllib.request

REMOTE_URL = os.environ.get("REMOTE_GNINA_URL", "").rstrip("/")
REMOTE_TOKEN = os.environ.get("REMOTE_GNINA_TOKEN", "")
TIMEOUT = float(os.environ.get("REMOTE_GNINA_TIMEOUT", "45"))

_health_cache = {"value": None, "checked": False}


def _remote_configured():
    return bool(REMOTE_URL)


def _remote_health():
    """Cached for the process lifetime -- a health probe on every single
       available() call (called often, e.g. to decide whether to show
       the GNINA toggle in Advanced Settings at all) would be wasteful
       and slow down the UI on every docking-tab render."""
    if _health_cache["checked"]:
        return _health_cache["value"]
    ok = False
    if _remote_configured():
        try:
            req = urllib.request.Request(f"{REMOTE_URL}/health")
            with urllib.request.urlopen(req, timeout=5) as resp:
                ok = bool(json.loads(resp.read()).get("available", False))
        except Exception:
            ok = False
    _health_cache["value"] = ok
    _health_cache["checked"] = True
    return ok


def patched_available(self):
    """self is a GninaRescorer instance -- self.binary is its configured
       binary name (default "gnina"), matching the original method."""
    if shutil.which(self.binary) is not None:
        return True
    return _remote_health()


def make_patched_rescore(original_rescore):
    """original_rescore: the UNPATCHED GninaRescorer.rescore function,
       captured by desktop.py before patching -- called unchanged
       whenever a local binary genuinely is present, so local-gnina
       behaviour (e.g. testing this desktop build on Linux) is
       byte-for-byte identical to the plain web deployment."""

    def patched_rescore(self, receptor_pdb, pose_mol):
        if shutil.which(self.binary) is not None:
            return original_rescore(self, receptor_pdb, pose_mol)
        if not _remote_configured():
            return None
        return _remote_rescore(receptor_pdb, pose_mol)

    return patched_rescore


def _remote_rescore(receptor_pdb_path, pose_mol):
    import tempfile

    from rdkit import Chem

    with open(receptor_pdb_path) as f:
        receptor_text = f.read()

    # Match the ORIGINAL local method's own ligand-serialization exactly
    # (Chem.SDWriter to a real .sdf file) rather than a MolBlock string,
    # so the remote server's Chem.SDMolSupplier parses it identically.
    with tempfile.TemporaryDirectory() as td:
        lig_path = os.path.join(td, "pose.sdf")
        w = Chem.SDWriter(lig_path)
        w.write(pose_mol)
        w.close()
        with open(lig_path) as f:
            ligand_sdf = f.read()

    body = json.dumps({"receptor_pdb": receptor_text, "ligand_sdf": ligand_sdf}).encode()
    req = urllib.request.Request(
        f"{REMOTE_URL}/rescore",
        data=body,
        headers={"Content-Type": "application/json", "Authorization": f"Bearer {REMOTE_TOKEN}"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        return {"error": f"remote GNINA service returned {e.code}: {e.read().decode(errors='replace')[:300]}"}
    except Exception as e:
        return {"error": f"could not reach remote GNINA service: {e}"}
