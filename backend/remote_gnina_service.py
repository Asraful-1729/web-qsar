"""
Remote GNINA rescoring service — runs on a machine that actually has a
working GPU + gnina (this server), so the Windows DESKTOP build (which
has no official GNINA binary to bundle at all — confirmed empirically:
every GNINA GitHub release, v1.0 through v1.3.3, ships Linux-only ELF
binaries, and qsar-desktop's own prior build effort independently found
and documented the same thing) can still offer real CNN rescoring as an
OPTIONAL network call, instead of skipping GNINA on Windows entirely.

Deliberately a SEPARATE process/file from backend/app.py, not a new route
added there — this is desktop-conversion infrastructure, not a change to
the existing served API surface. It reuses docking.engines.GninaRescorer
UNMODIFIED (imported, not copied/reimplemented) — this file is a thin
HTTP wrapper around that existing, already-validated class.

Auth: a single shared bearer token (REMOTE_GNINA_TOKEN env var, generated
once — see desktop/REMOTE_GNINA_SETUP.md) — this endpoint runs arbitrary
docking computations on a real GPU server, so it must not be left open to
the internet unauthenticated even though the payloads themselves are
just chemistry data.

Run:  REMOTE_GNINA_TOKEN=<token> uvicorn remote_gnina_service:app --app-dir backend --host 0.0.0.0 --port 8300
"""
import os
import tempfile

from fastapi import FastAPI, Header, HTTPException
from pydantic import BaseModel

TOKEN = os.environ.get("REMOTE_GNINA_TOKEN")

app = FastAPI(title="Remote GNINA rescoring")


class RescoreBody(BaseModel):
    receptor_pdb: str   # full PDB file text
    ligand_sdf: str      # full SDF file text (single pose)


def _check_auth(authorization):
    if not TOKEN:
        raise HTTPException(500, "REMOTE_GNINA_TOKEN not set on this server — refusing to serve unauthenticated.")
    expected = f"Bearer {TOKEN}"
    if authorization != expected:
        raise HTTPException(401, "invalid or missing token")


@app.get("/health")
def health():
    """No auth needed -- just confirms the service is up and whether
       gnina is actually usable here (the desktop client's own
       availability check calls this, not /rescore, to avoid spending a
       real docking computation just to probe reachability)."""
    from docking.engines import GninaRescorer
    return {"available": GninaRescorer().available()}


@app.post("/rescore")
def rescore(body: RescoreBody, authorization: str = Header(None)):
    _check_auth(authorization)
    from docking.engines import GninaRescorer  # reused unmodified
    from rdkit import Chem

    rescorer = GninaRescorer()
    if not rescorer.available():
        raise HTTPException(503, "gnina is not available on this server right now")

    with tempfile.TemporaryDirectory() as td:
        receptor_path = os.path.join(td, "receptor.pdb")
        ligand_path = os.path.join(td, "pose.sdf")
        with open(receptor_path, "w") as f:
            f.write(body.receptor_pdb)
        with open(ligand_path, "w") as f:
            f.write(body.ligand_sdf)

        supplier = Chem.SDMolSupplier(ligand_path)
        mol = next(iter(supplier), None)
        if mol is None:
            raise HTTPException(400, "could not parse ligand_sdf as a valid molecule")

        result = rescorer.rescore(receptor_path, mol)
        if result is None:
            raise HTTPException(503, "gnina rescoring failed on this server")
        return result
