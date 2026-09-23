"""
============================================================
  PhytoScreen — Desktop serving backend (app.py)
============================================================
  Prediction-first serving app. NO training here (CLAUDE.md §1). Loads
  the real target buckets in models/<target_id>/ (AutoGluon + Chemprop,
  see serving/model_adapter.py) — the app never touches a registry.json
  or raw pickle asset.

  One chosen model per target — no multi-model comparison UI (CLAUDE.md
  §13 guardrail).

  Run:
      uvicorn app:app --host 127.0.0.1 --port 8000
      open http://localhost:8000/

  Tabs served by one single-page UI: Predict / ADMET / Compare / Docking
  / Screen / Target Info.
============================================================
"""
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field
from typing import Optional, List

from serving import model_adapter as MA
import analysis as A
import admet as ADMET
import factory_browser
import downloads

app = FastAPI(title="PhytoScreen", version="3.0")
app.include_router(factory_browser.router)
app.include_router(downloads.router)


class JobCancelled(Exception):
    """Raised from inside a background job's progress callback / per-item
       loop once its job dict's 'cancel_requested' flag has been set (via
       a POST .../cancel/{jid} endpoint) — the owning work() function
       catches this SEPARATELY from a real error, marking the job
       'cancelled' rather than 'error'. Python threads can't be killed
       safely; this cooperative check at natural per-item checkpoints
       (one compound, one decoy, one downloaded chunk) is the standard
       safe way to stop a long-running background job partway through."""


def _check_cancelled(job):
    if job.get("cancel_requested"):
        raise JobCancelled()

# The UI (frontend/) is a separate app/process — dev server (Vite) or a
# built static bundle served by anything, no longer served by this backend
# (see the removed desktop.py / static-mount). No cookies or session auth
# exist here, so a permissive default is safe; ALLOWED_ORIGINS lets a real
# deployment lock it down to its actual frontend origin(s).
_allowed_origins = [o.strip() for o in os.environ.get("ALLOWED_ORIGINS", "*").split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_allowed_origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

DISCLAIMER = ("Prioritisation aid, not a substitute for assays. Trust predictions only for "
              "in-domain molecules; treat the top of the list as a shortlist.")


def _num(v):
    try:
        if v is None or pd.isna(v):
            return None
        return round(float(v), 3)
    except Exception:
        return None


def _rows(df):
    """One dict per row. Out-of-domain rows NEVER carry a potency number —
       AD gating is enforced here, not just in the UI (CLAUDE.md §2)."""
    out = []
    for _, r in df.iterrows():
        in_ad = bool(r["In_AD"])
        out.append({
            "input_smiles": r.get("Input_SMILES"),
            "smiles": r.get("Standardised_SMILES"),
            "parsed_ok": bool(r.get("Parsed_OK")),
            "predicted_pIC50": _num(r.get("Predicted_pIC50")) if in_ad else None,
            "in_domain": in_ad,
            "ad_z": _num(r.get("AD_z")),
            "confidence": r.get("Confidence"),
            "confidence_label": r.get("Confidence_Label"),
            "confidence_basis": r.get("Confidence_Basis"),
        })
    return out


def _predict(target_id, smiles):
    if not smiles:
        raise HTTPException(400, "No SMILES provided.")
    try:
        target = MA.load_target(target_id)
    except MA.BucketError as e:
        raise HTTPException(409, str(e))
    except KeyError:
        raise HTTPException(404, f"Unknown target '{target_id}'")

    df = target.predict_smiles(smiles)
    unparsed = df[~df["Parsed_OK"]]
    parsed = df[df["Parsed_OK"]]
    in_dom = parsed[parsed["In_AD"]].sort_values("Predicted_pIC50", ascending=False).reset_index(drop=True)
    in_dom.insert(0, "Rank", range(1, len(in_dom) + 1))
    out_dom = parsed[~parsed["In_AD"]].reset_index(drop=True)

    in_rows = _rows(in_dom)
    for i, row in enumerate(in_rows):
        row["rank"] = i + 1

    return {
        "target": {"id": target.target_id, "name": target.name},
        "model": target.metrics.get("Best_Model") or target.metrics.get("best_model"),
        "model_metrics": {
            "test_r2": target.metrics.get("R2_Test"),
            "test_rmse": target.metrics.get("RMSE_Test"),
            "pearson_r": target.metrics.get("Pearson_r"),
            "ad_coverage_pct": target.metrics.get("AD_Coverage_pct"),
            "tropsha_pass": target.metrics.get("Tropsha_Pass"),
            "y_random_delta_r2": target.metrics.get("Y_Random_DeltaR2"),
        },
        "counts": {"in_domain": int(len(in_dom)), "out_of_domain": int(len(out_dom)),
                   "skipped": int(len(unparsed)), "submitted": int(len(smiles))},
        "in_domain": in_rows,
        "out_of_domain": _rows(out_dom),
        "skipped": [str(s) for s in unparsed["Input_SMILES"].tolist()],
        "ranked_by": "predicted_pIC50",
        "disclaimer": DISCLAIMER,
    }


# ---------------- models ----------------
class PredictBody(BaseModel):
    target_id: str
    smiles: List[str]


class MultiBody(BaseModel):
    smiles: List[str]
    target_ids: Optional[List[str]] = None


@app.get("/api/health")
def health():
    ids = MA.list_target_ids()
    try:
        dock_ready = bool(DOCK_AVAIL and DOCK_AVAIL.status()["ready"])
    except Exception:
        dock_ready = False
    return {"targets_in_bucket_dir": len(ids), "targets_dir": MA.TARGETS_DIR,
            "admet_ai": ADMET.learned_status()["available"],
            "docking": "ready" if dock_ready else "not_ready", "disclaimer": DISCLAIMER}


@app.get("/api/health/ml_backends")
def health_ml_backends():
    """Whether AutoGluon's base-learner libraries can actually load their
       native compute engine (xgboost.dll/lib_lightgbm.dll/etc), not just
       their Python module. `import xgboost` alone triggers xgboost's
       ctypes CDLL() load of that native library at module-import time
       (see xgboost/libpath.py's find_lib_path()) — invisible to
       PyInstaller's static analysis (a plain --hidden-import bundles the
       .py files but not the runtime-probed .dll), so a frozen build can
       import fine and /api/health can report OK while every prediction
       against an XGBoost-based AutoGluon ensemble member crashes.
       WHICH real target exercises this depends on AutoGluon's per-target
       model composition (not every target's chosen ensemble includes
       every base-learner type), so a single smoke-tested /api/predict
       call can pass while this is still broken — see BUILD_WINDOWS.md's
       "Why --collect-all lightgbm/catboost/xgboost". CI calls this
       directly instead of hoping the right target gets tested."""
    out = {}
    for lib in ("xgboost", "lightgbm", "catboost"):
        try:
            __import__(lib)
            out[lib] = {"ok": True}
        except Exception as e:
            out[lib] = {"ok": False, "error": str(e)}
    return out


@app.get("/api/targets")
def targets():
    out = []
    for meta in MA.list_targets_meta():
        m = meta["metrics"]
        out.append({"target_id": meta["target_id"], "name": meta["target_id"],
                    "best_model": m.get("Best_Model") or m.get("best_model"),
                    "n_compounds": m.get("Total_N") or (m.get("counts") or {}).get("total"),
                    "test_r2": m.get("R2_Test"), "test_rmse": m.get("RMSE_Test"),
                    "ad_coverage_pct": m.get("AD_Coverage_pct"), "tropsha_pass": m.get("Tropsha_Pass")})
    out.sort(key=lambda x: x["target_id"])
    return {"targets": out}


@app.post("/api/predict")
def predict(body: PredictBody):
    return _predict(body.target_id, [s.strip() for s in body.smiles if s and s.strip()])


@app.post("/api/predict_csv")
async def predict_csv(target_id: str = Form(...), file: UploadFile = File(...)):
    import io
    df = pd.read_csv(io.BytesIO(await file.read()))
    col = "SMILES" if "SMILES" in df.columns else ("smiles" if "smiles" in df.columns else df.columns[0])
    smiles = [str(s).strip() for s in df[col].dropna().tolist() if str(s).strip()]
    return _predict(target_id, smiles)


@app.post("/api/parse_sdf")
async def parse_sdf(file: UploadFile = File(...)):
    """SDF -> SMILES list, for the SDF input mode on every input tab (Screen/
       Predict/ADMET/Compare/Docking) — RDKit is far more reliable at reading
       an SDF's bond/stereo perception server-side than any client-side JS
       parser would be, so every tab funnels SDF uploads through here and
       then treats the result exactly like a pasted SMILES list."""
    import io
    from rdkit import Chem
    data = await file.read()
    if not data:
        raise HTTPException(400, "Empty file.")
    suppl = Chem.ForwardSDMolSupplier(io.BytesIO(data), removeHs=False, sanitize=True)
    smiles, names, n_skipped = [], [], 0
    for mol in suppl:
        if mol is None:
            n_skipped += 1
            continue
        try:
            smiles.append(Chem.MolToSmiles(mol))
        except Exception:
            n_skipped += 1
            continue
        names.append(mol.GetProp("_Name") if mol.HasProp("_Name") else None)
    if not smiles:
        raise HTTPException(400, "No parsable molecules found in this SDF file.")
    return {"smiles": smiles, "names": names, "n_parsed": len(smiles), "n_skipped": n_skipped}


@app.post("/api/predict_multi")
def predict_multi(body: MultiBody):
    smiles = [s.strip() for s in body.smiles if s and s.strip()]
    if not smiles:
        raise HTTPException(400, "No SMILES provided.")
    if not body.target_ids:
        raise HTTPException(400, "Provide target_ids.")
    target_ids = body.target_ids
    result = A.analyse(smiles, target_ids)
    if not result["targets"]:
        raise HTTPException(409, f"No usable targets among {target_ids}.")
    return result


ADMET_SYNC_MAX = 50          # lists this size or smaller run synchronously
_ADMET_JOBS = {}             # main-app job store: id -> {det, worker_job, ...}


@app.post("/api/admet")
def admet(body: PredictBody):
    import uuid
    smiles = [s.strip() for s in body.smiles if s and s.strip()]
    if not smiles:
        raise HTTPException(400, "No SMILES provided.")
    det = ADMET.deterministic_profiles(smiles)
    status = ADMET.learned_status()
    std = [p["standardised_smiles"] for p in det if p.get("parsed_ok")]

    # no worker, or nothing valid -> deterministic only, immediately
    if not status["available"] or not std:
        for p in det:
            p["learned"] = {"available": False, "note": status.get("note")}
        return {"mode": "result", "profiles": det, "learned": status, "disclaimer": DISCLAIMER}

    # small -> synchronous
    if len(std) <= ADMET_SYNC_MAX:
        preds = ADMET.worker_profile(std)
        if preds is None:
            for p in det:
                p["learned"] = {"available": False, "note": "ADMET-AI worker did not respond."}
        else:
            ADMET.attach_learned(det, preds)
        return {"mode": "result", "profiles": det, "learned": status, "disclaimer": DISCLAIMER}

    # large -> async job (isolated worker), UI polls
    wj = ADMET.worker_submit(std)
    if wj is None:
        for p in det:
            p["learned"] = {"available": False, "note": "ADMET-AI worker did not accept the job."}
        return {"mode": "result", "profiles": det, "learned": status, "disclaimer": DISCLAIMER}
    jid = uuid.uuid4().hex[:12]
    _ADMET_JOBS[jid] = {"det": det, "worker_job": wj}
    return {"mode": "job", "job_id": jid, "total": len(std), "disclaimer": DISCLAIMER}


@app.get("/api/admet/job/{jid}")
def admet_job(jid: str):
    job = _ADMET_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    st = ADMET.worker_poll(job["worker_job"])
    if st is None:
        raise HTTPException(502, "ADMET-AI worker unreachable")
    if st["status"] == "done":
        ADMET.attach_learned(job["det"], st.get("predictions") or {})
        _ADMET_JOBS.pop(jid, None)
        return {"status": "done", "profiles": job["det"], "learned": {"available": True, "source": "ADMET-AI"},
                "disclaimer": DISCLAIMER}
    if st["status"] == "error":
        _ADMET_JOBS.pop(jid, None)
        return {"status": "error", "error": st.get("error")}
    return {"status": st["status"], "done": st.get("done", 0), "total": st.get("total", 0)}


# ---------------- docking (availability-gated) ----------------
try:
    from docking import availability as DOCK_AVAIL
    from docking import pipeline as DOCK_PIPE
    from docking import profile as DOCK_PROFILE
    from docking import recommend as DOCK_RECOMMEND
    _DOCK_IMPORT_ERR = None
except Exception as e:
    DOCK_AVAIL = DOCK_PIPE = DOCK_PROFILE = DOCK_RECOMMEND = None
    _DOCK_IMPORT_ERR = str(e)

_DOCK_JOBS = {}


@app.get("/api/docking/status")
def docking_status():
    if DOCK_AVAIL is None:
        return {"ready": False, "import_error": _DOCK_IMPORT_ERR,
                "note": "Docking package not importable — install its dependencies (rdkit, meeko, posebusters)."}
    st = DOCK_AVAIL.status()
    st["planned"] = ["AutoDock Vina docking + PoseBusters physical-validity gate",
                     "GNINA CNN rescoring (second opinion, optional)",
                     "2D protein-ligand interaction diagram"]
    try:
        reg = DOCK_PROFILE.load_registry()
        st["docking_targets"] = list(reg.keys())
        st["target_details"] = [{
            "target_id": t, "name": r.get("name", t),
            "site_source": r.get("site_source"),
        } for t, r in reg.items()]
    except Exception:
        st["docking_targets"] = []
        st["target_details"] = []
    return st


@app.get("/api/docking/targets")
def docking_targets():
    if DOCK_PROFILE is None:
        return {"targets": []}
    reg = DOCK_PROFILE.load_registry()
    return {"targets": [{"target_id": t, "name": r.get("name", t)} for t, r in reg.items()]}


@app.get("/api/diseases")
def diseases():
    """Diseases (from panel_results_v2.csv) that have >=1 associated target
       among our QSAR-modeled targets — most of the CSV's diseases fall
       away here since it covers ~669 genes and we only model ~52."""
    if DOCK_RECOMMEND is None:
        return {"diseases": []}
    return {"diseases": DOCK_RECOMMEND.list_diseases()}


@app.get("/api/diseases/{disease_id}/targets")
def disease_targets(disease_id: str):
    """Our targets associated with disease_id, ranked by disease-association
       score, annotated with real QSAR quality — nothing here is invented."""
    if DOCK_RECOMMEND is None:
        return {"targets": []}
    return {"targets": DOCK_RECOMMEND.targets_for_disease(disease_id)}


@app.get("/api/targets/{target_id}/recommendation")
def target_recommendation(target_id: str):
    """The 'Why this?' evidence bundle for one target's recommended docking
       structure: our own empirical validation plus panel_results_v2.csv's
       crystallographic context."""
    if DOCK_RECOMMEND is None:
        raise HTTPException(503, "Docking package not available")
    rec = DOCK_RECOMMEND.recommendation(target_id)
    if rec is None:
        raise HTTPException(404, f"'{target_id}' is not one of our QSAR-modeled targets")
    return rec


@app.get("/api/docking/receptor/{target_id}")
def docking_receptor(target_id: str):
    """Raw receptor PDB text, for the 3D pose viewer (protein context around
       the docked ligand). Same file docking already uses for PoseBusters/PLIP."""
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    profile = DOCK_PROFILE.load_profile(target_id)
    receptor_pdb = profile.get("receptor_pdb")
    if not receptor_pdb or not os.path.exists(receptor_pdb):
        raise HTTPException(404, f"no receptor PDB for '{target_id}'")
    return PlainTextResponse(open(receptor_pdb).read())


@app.get("/api/docking/receptor_file")
def docking_receptor_file(path: str):
    """Serves a receptor PDB by its own file path rather than by target_id —
       used for the binding-site 3D preview and the dock-results pose viewer
       whenever the receptor in play isn't (or might not be) the registry's
       automatic default: an Advanced Settings custom structure (any
       target, not just no-QSAR-model ones — its receptor_pdb is an
       absolute path already known client-side from the custom-receptor
       job's response, not a registry target_id) or the exact receptor a
       completed docking job actually used. Restricted to docking_targets/
       — never an arbitrary filesystem path — so this can't be used to read
       anything else on the box."""
    allowed_root = os.path.abspath("docking_targets")
    real = os.path.abspath(path)
    if not (real == allowed_root or real.startswith(allowed_root + os.sep)) or not os.path.exists(real):
        raise HTTPException(403, "path not allowed")
    return PlainTextResponse(open(real).read())


def _structure_candidates_for_gene(gene, default_pdb=None):
    """Shared by both structure_candidates() (QSAR-modeled targets) and
       gene_structure_candidates() (any protein in the Version 2 CSV, model
       or not) — everything here only needs a gene symbol, never a real
       QSAR target_id."""
    from scripts import panel_candidates as PC
    df = PC._panel_df()
    if df is None or gene not in df.index:
        return {"gene": gene, "default_pdb_id": default_pdb, "candidates": [], "n_qualifying_structures": 0}
    row = df.loc[gene]
    ranked_ids = [p for p in str(row.get("all_pdb_ids_ranked") or "").split(";") if p]
    cands = PC._parse_top5(row.get("top5_pdb_summary"))
    import re as _re
    quality_re = _re.compile(r"^(\S+)\s+\(res=([\d.]+),\s*RSCC=([\d.]+),\s*RSR=([\d.]+)\)")
    quality_by_pdb = {}
    for chunk in str(row.get("top5_pdb_summary") or "").split(" | "):
        m = quality_re.match(chunk.strip())
        if m:
            quality_by_pdb[m.group(1)] = {"resolution": float(m.group(2)), "ligand_RSCC": float(m.group(3)),
                                          "ligand_RSR": float(m.group(4))}
    out = []
    for c in cands:
        try:
            rank = ranked_ids.index(c["pdb_id"]) + 1
        except ValueError:
            rank = None
        q = quality_by_pdb.get(c["pdb_id"], {})
        out.append({"pdb_id": c["pdb_id"], "resname": c["resname"], "csv_rank": rank,
                    "is_current_default": c["pdb_id"] == default_pdb,
                    "resolution": q.get("resolution"), "ligand_RSCC": q.get("ligand_RSCC"),
                    "ligand_RSR": q.get("ligand_RSR")})
    note = "Crystallographic quality ranking only — not a guarantee of pose accuracy for docking."
    # _parse_top5 already drops any candidate whose ONLY ligand is a
    # buffer/cryoprotectant/ion (ADDITIVE_BLACKLIST) or under the MW floor
    # — for ~70 genes across panel_results_v2.csv (e.g. PLEC/2ODV, whose
    # sole HETATM is S-1,2-propanediol) that leaves `out` empty even though
    # a real, good-quality structure exists (best_pdb_id/best_chain were
    # already picked by quality alone, independent of ligand validity).
    # Surface it as ONE blind_only candidate (resname omitted) rather than
    # leaving the picker empty with no path forward — build_receptor()'s
    # no-ligand fallback (see docking/receptor_prep.py) can build a real
    # receptor from it for Blind docking, just never a site-specific one.
    best_pdb = row.get("best_pdb_id")
    if not out and best_pdb and str(best_pdb).strip():
        best_pdb = str(best_pdb).strip().upper()
        try:
            rank = ranked_ids.index(best_pdb) + 1
        except ValueError:
            rank = None
        q = quality_by_pdb.get(best_pdb, {})
        out.append({"pdb_id": best_pdb, "resname": None, "csv_rank": rank,
                    "is_current_default": best_pdb == default_pdb,
                    "resolution": q.get("resolution"), "ligand_RSCC": q.get("ligand_RSCC"),
                    "ligand_RSR": q.get("ligand_RSR"), "blind_only": True})
        note = ("No candidate here has a validated small-molecule ligand (only buffers/cryoprotectants/ions, "
                "if anything was co-crystallized) — showing the best-quality structure on record for Blind "
                "docking only; it won't get a site-specific default.")
    return {"gene": gene, "default_pdb_id": default_pdb, "candidates": out,
           "n_qualifying_structures": row.get("n_qualifying_structures"),
           "note": note}


_PLOT_MEDIA_TYPES = {"svg": "image/svg+xml", "png": "image/png", "tif": "image/tiff", "tiff": "image/tiff", "pdf": "application/pdf"}


@app.get("/api/targets/{target_id}/structure_candidates")
def structure_candidates(target_id: str):
    """Every qualifying structure for this target's gene from the Version 2
       CSV (panel_results_v2.csv), for Advanced Settings' manual structure
       picker — resolution/RSCC/RSR/ligand per candidate, plus which one is
       the current automatic default, so 'View Evidence' works per-row, not
       just for the top pick."""
    if DOCK_RECOMMEND is None:
        raise HTTPException(503, "Docking package not available")
    t2g, _ = DOCK_RECOMMEND._target_gene_map()
    gene = t2g.get(target_id)
    if gene is None:
        raise HTTPException(404, f"'{target_id}' is not one of our QSAR-modeled targets")
    default_pdb = None
    reg = DOCK_PROFILE.load_registry()
    src = (reg.get(target_id) or {}).get("pdb_source")
    if src:
        default_pdb = src.split("_raw")[0].split(".")[0].upper()
    return {"target_id": target_id, **_structure_candidates_for_gene(gene, default_pdb)}


@app.get("/api/genes/{gene_symbol}/structure_candidates")
def gene_structure_candidates(gene_symbol: str):
    """Same as structure_candidates(), but for ANY protein in the Version 2
       CSV — including the ~620 genes that have disease association evidence
       but no trained QSAR model. Powers 'docking only, no QSAR model'
       targets: the Screen tab's disease->target list no longer hides these,
       it routes them here instead so a receptor can still be built and
       docked against on demand (see /api/docking/receptor/custom, which
       already accepts any string as target_id — a QSAR model was never a
       structural requirement for it, only a UI convention)."""
    if DOCK_RECOMMEND is None:
        raise HTTPException(503, "Docking package not available")
    default_pdb = None
    reg = DOCK_PROFILE.load_registry()
    src = (reg.get(f"GENE_{gene_symbol.upper()}") or {}).get("pdb_source")
    if src:
        default_pdb = src.split("_raw")[0].split(".")[0].upper()
    return _structure_candidates_for_gene(gene_symbol.upper(), default_pdb)


@app.get("/api/targets/{target_id}/binding_site")
def binding_site(target_id: str):
    """Binding-site evidence for the AUTOMATIC (registry-default) profile:
       the pocket residues within 5 A of the reference ligand, plus the
       center/box_size already stored. Computed fresh from files already on
       disk (raw PDB + cleaned receptor) every call — cheap (no PDBFixer/
       Vina involved) and works retroactively for every registry target
       without needing to rebuild anything."""
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    try:
        profile = DOCK_PROFILE.load_profile(target_id)
    except Exception as e:
        raise HTTPException(404, str(e))
    from docking import receptor_prep as RP
    resname = profile.get("reference_ligand_resname")
    raw_path = os.path.join(DOCK_PROFILE.DOCKING_TARGETS_DIR, target_id, profile.get("pdb_source") or "")
    residues, error = [], None
    if not resname:
        # site_source == "none_validated": no real co-crystallized ligand
        # exists for this target (see scripts/audit_ligand_selection.py) —
        # distinct from the file-missing case below, and a different
        # message so the UI doesn't imply a download would fix it.
        error = "no validated small-molecule binding site for this target — use Blind mode"
    elif os.path.exists(raw_path) and profile.get("receptor_pdb") and os.path.exists(profile["receptor_pdb"]):
        try:
            ref_coords = RP.locate_ligand_near(raw_path, resname, profile["center"])
            residues = RP.pocket_residues(profile["receptor_pdb"], ref_coords, cutoff=5.0)
        except Exception as e:
            error = str(e)
    else:
        error = "raw structure or cleaned receptor no longer on disk"
    crystal_sdf = os.path.join(DOCK_PROFILE.DOCKING_TARGETS_DIR, target_id, "crystal_ligand.sdf")
    blind_center = blind_box_size = None
    all_res = []
    if profile.get("receptor_pdb") and os.path.exists(profile["receptor_pdb"]):
        try:
            blind_center, blind_box_size = RP.box_from_receptor(profile["receptor_pdb"])
        except Exception:
            pass   # blind mode just won't be offered for this target; site-specific evidence above is unaffected
        try:
            # The full receptor residue list — for MANUAL site definition,
            # which must not be limited to the automatically-detected
            # pocket neighborhood the same way `residues` above is.
            all_res = RP.all_residues(profile["receptor_pdb"])
        except Exception:
            pass   # manual mode just won't have a full list to pick from for this target
    return {"target_id": target_id, "center": profile.get("center"), "box_size": profile.get("box_size"),
           "reference_ligand_resname": resname, "pocket_residues": residues, "n_pocket_residues": len(residues),
           "all_residues": all_res, "has_reference_ligand_mol": os.path.exists(crystal_sdf),
           "blind_center": blind_center, "blind_box_size": blind_box_size,
           "error": error if not residues else None}


@app.get("/api/targets/{target_id}/reference_ligand.sdf")
def reference_ligand_sdf(target_id: str):
    """The crystal (experimentally-observed) reference-ligand pose, real
       bond orders included — only exists for targets that went through the
       full validate_target.py pipeline (not Advanced Settings' on-demand
       custom structures, which skip this step)."""
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    path = os.path.join(DOCK_PROFILE.DOCKING_TARGETS_DIR, target_id, "crystal_ligand.sdf")
    if not os.path.exists(path):
        raise HTTPException(404, f"no reference-ligand pose on file for '{target_id}'")
    return PlainTextResponse(open(path).read())


    # NOTE: redocking-pose validation is no longer a separate user action —
    # it now runs automatically at docking/screen submission time, using
    # whichever PDB/exhaustiveness/n_poses the user already picked for
    # their compounds (see docking/pipeline.py's redock_reference_for_
    # profile, wired into _run_docking_job below and serving/screen.py's
    # run()). docking_receptor_custom (below) builds the crystal reference
    # pose needed for that alongside the receptor itself, for whichever
    # ligand the user explicitly picks from /api/docking/pdb_ligands —
    # there is no automatic "main ligand" selection any more.


class BoxFromResiduesBody(BaseModel):
    target_id: str
    residues: List[dict]     # [{chain, resnum}, ...] — from a binding_site/custom_profile response
    receptor_pdb: Optional[str] = None   # pass a custom_profile's receptor_pdb to target that structure instead
    padding: float = 8.0


@app.post("/api/docking/box_from_residues")
def box_from_residues(body: BoxFromResiduesBody):
    """Advanced Settings' 'define the binding site from residues' path:
       user checks residues off the pocket list, this recomputes center/
       box_size from their coordinates (same math as the automatic ligand-
       centered box) — an alternative to typing raw XYZ numbers."""
    if not body.residues:
        raise HTTPException(400, "Pick at least one residue.")
    from docking import receptor_prep as RP
    receptor_pdb = body.receptor_pdb
    if not receptor_pdb:
        if DOCK_PROFILE is None:
            raise HTTPException(503, "Docking package not available")
        try:
            receptor_pdb = DOCK_PROFILE.load_profile(body.target_id)["receptor_pdb"]
        except Exception as e:
            raise HTTPException(404, str(e))
    if not os.path.exists(receptor_pdb):
        raise HTTPException(404, f"receptor file not found: {receptor_pdb}")
    try:
        center, box_size = RP.box_from_residues(receptor_pdb, body.residues, padding=body.padding)
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"center": center, "box_size": box_size}


_CUSTOM_RECEPTOR_JOBS = {}


def _build_receptor_profile(raw_pdb, work_id, out_dir, target_id, name, ref_resname, chain, resnum=None, progress=None):
    """Shared core of 'build a receptor centered on ONE specific ligand':
       used by docking_receptor_custom below (which fetches raw_pdb from
       RCSB first) and by docking_alternate_ligand_build (which already has
       raw_pdb sitting on disk from the ORIGINAL job, no fetch needed — see
       that endpoint's docstring). Never touches docking_registry.json —
       this is always a per-request profile, not a promotion to the
       registry default (see docking_receptor_custom's own docstring for
       why that matters)."""
    import time
    from docking import receptor_prep as RP

    t_build0 = time.time()
    profile = RP.build_receptor(raw_pdb, work_id, name=name, ref_resname=ref_resname, chain=chain, resnum=resnum,
                                out_dir=out_dir, progress=progress)
    profile["target_id"] = target_id   # advertise the REAL target_id to the caller, not the scratch work_id
    # Original, unmodified structure — lets the UI show a before/after
    # comparison against the stripped/repaired receptor_pdb.
    profile["raw_pdb_path"] = os.path.abspath(raw_pdb)
    build_seconds = round(time.time() - t_build0, 1)
    profile["timing"] = {"build_seconds": build_seconds}
    return profile


def _cached_raw_pdb(pdb_id):
    """Shared per-PDB cache (docking_targets/_pdb_cache/<PDB_ID>/raw.pdb) —
       a PDB entry's own atoms/ligands never depend on which target it's
       being considered for, so /api/docking/pdb_ligands and
       docking_receptor_custom share one fetch instead of each target
       re-downloading the same structure from RCSB."""
    return os.path.join("docking_targets", "_pdb_cache", pdb_id.upper(), "raw.pdb")


_PDB_LIGANDS_JOBS = {}


class PdbLigandsBody(BaseModel):
    pdb_id: str


@app.post("/api/docking/pdb_ligands")
def docking_pdb_ligands(body: PdbLigandsBody):
    """Every real co-crystallized ligand in a PDB entry (see
       receptor_prep.list_ligands's own docstring) — lets the user pick
       WHICH one to build a receptor/box around, instead of the app
       silently choosing "the largest non-additive HETATM group" for them.
       There is deliberately no "main"/pre-selected ligand in the response
       — which one matters depends on the biology of that specific
       structure, not something a generic heuristic should decide."""
    import threading, uuid
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    pdb_id = body.pdb_id.strip().upper()
    if not pdb_id:
        raise HTTPException(400, "pdb_id is required")
    jid = uuid.uuid4().hex[:12]
    _PDB_LIGANDS_JOBS[jid] = {"status": "queued", "step": "Queued", "ligands": None, "error": None}

    def work():
        job = _PDB_LIGANDS_JOBS[jid]
        job["status"] = "running"
        try:
            from docking import receptor_prep as RP
            raw_pdb = _cached_raw_pdb(pdb_id)
            if not os.path.exists(raw_pdb):
                from scripts.pdb_fetch import fetch_pdb
                os.makedirs(os.path.dirname(raw_pdb), exist_ok=True)
                job["step"] = f"Fetching {pdb_id} from RCSB"
                fetch_pdb(pdb_id, raw_pdb)
            job["step"] = "Listing co-crystallized ligands"
            job["ligands"] = RP.list_ligands(raw_pdb)
            job["status"] = "done"
        except Exception as e:
            job["status"] = "error"; job["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"job_id": jid}


@app.get("/api/docking/pdb_ligands/job/{jid}")
def docking_pdb_ligands_job(jid: str):
    job = _PDB_LIGANDS_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    r = {"status": job["status"], "step": job.get("step")}
    if job["status"] == "done":
        r["ligands"] = job["ligands"]
    if job["status"] == "error":
        r["error"] = job.get("error")
    return r


class CustomReceptorBody(BaseModel):
    target_id: str
    pdb_id: str
    chain: Optional[str] = None
    ligand_resname: Optional[str] = None


@app.post("/api/docking/receptor/custom")
def docking_receptor_custom(body: CustomReceptorBody):
    """Advanced Settings' manual structure override: build (strip/repair/
       PDBQT) a receptor for an expert-chosen PDB entry AND ligand on
       demand. Never touches docking_registry.json regardless of the
       outcome — this is still a per-request profile the caller passes
       back as advanced.custom_profile on submit, so it can't silently
       overwrite the automatic default (or race a concurrent
       batch_validate.py run's registry writes).

       Also builds a bond-order-correct crystal reference pose for
       body.ligand_resname (profile["crystal_sdf"]/["reference_smiles"]) —
       needed by docking/pipeline.py's redock_reference_for_profile, which
       runs automatically at submission time using whichever ligand was
       picked here. Building it is non-fatal if it fails (e.g. RCSB has no
       SMILES on record for that ligand): compound docking still proceeds,
       just without a redocking-validation number for this structure.

       Cached per (target_id, pdb_id, ligand_resname) — not rebuilt on
       every pick of the same structure+ligand, since PDBFixer's own
       hydrogenation step is NOT deterministic (two independent builds
       from the identical raw PDB place atoms at measurably different
       coordinates), which would otherwise make redocking-validation RMSD
       drift between repeated runs of the same pick."""
    import threading, uuid
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")

    jid = uuid.uuid4().hex[:12]
    _CUSTOM_RECEPTOR_JOBS[jid] = {"status": "queued", "step": "Queued", "profile": None, "error": None}

    def work():
        job = _CUSTOM_RECEPTOR_JOBS[jid]
        job["status"] = "running"
        def _step(label):
            job["step"] = label
        try:
            import json, shutil, time
            out_dir = os.path.join("docking_targets", "_custom")
            work_id = f"{body.target_id}__{body.pdb_id}__{body.ligand_resname or 'blind'}"
            work_dir = os.path.join(out_dir, work_id)
            profile_json = os.path.join(work_dir, "profile.json")
            if os.path.exists(profile_json):
                _step(f"Reusing previously-built {body.pdb_id} receptor")
                with open(profile_json) as f:
                    profile = json.load(f)
            else:
                from scripts.detect_chain import chain_for_ligand
                os.makedirs(work_dir, exist_ok=True)
                raw_pdb = os.path.join(work_dir, f"{body.pdb_id}_raw.pdb")
                t0 = time.time()
                cached = _cached_raw_pdb(body.pdb_id)
                if os.path.exists(cached):
                    # Already fetched by /api/docking/pdb_ligands when the
                    # user listed this PDB's ligands — no need to hit RCSB
                    # again for the exact same structure.
                    _step(f"Reusing already-fetched {body.pdb_id}")
                    shutil.copyfile(cached, raw_pdb)
                else:
                    from scripts.pdb_fetch import fetch_pdb
                    _step(f"Fetching {body.pdb_id} from RCSB")
                    fetch_pdb(body.pdb_id, raw_pdb)
                fetch_seconds = round(time.time() - t0, 1)
                chain = body.chain
                if not chain and body.ligand_resname:
                    _step("Locating chain containing the reference ligand")
                    chain = chain_for_ligand(raw_pdb, body.ligand_resname)
                    if chain is None:
                        raise RuntimeError(f"ligand '{body.ligand_resname}' not found in any chain of {body.pdb_id}")
                profile = _build_receptor_profile(raw_pdb, work_id, out_dir, body.target_id,
                                                  f"{body.target_id} ({body.pdb_id}, manual)",
                                                  body.ligand_resname, chain, progress=_step)
                # fetch_pdb has no retry and a per-request 30s timeout; a slow/
                # degraded RCSB connection (or an mmCIF-only entry needing the
                # fetch-.pdb-404-then-fetch-.cif fallback) can silently cost
                # real seconds here with no other visibility — recorded
                # alongside build_receptor's own per-phase timing (see
                # _build_receptor_profile) so a slow run is diagnosable
                # with real numbers instead of one opaque "Preparing…" message.
                profile.setdefault("timing", {})["fetch_pdb_seconds"] = fetch_seconds
                if body.ligand_resname:
                    try:
                        from scripts.validate_target import fetch_ligand_smiles, make_crystal_sdf
                        _step("Fetching reference ligand SMILES from RCSB")
                        template_smiles = fetch_ligand_smiles(body.ligand_resname)
                        _step("Building crystal reference pose (correct bond orders)")
                        crystal_sdf = os.path.join(work_dir, "crystal_ligand.sdf")
                        make_crystal_sdf(raw_pdb, body.ligand_resname, crystal_sdf, template_smiles, chain=chain)
                        profile["crystal_sdf"] = os.path.abspath(crystal_sdf)
                        profile["reference_smiles"] = template_smiles
                    except Exception as e:
                        profile["crystal_sdf_error"] = str(e)
                with open(profile_json, "w") as f:
                    json.dump(profile, f)
            job["step"] = "Done"
            job["profile"] = profile
            job["status"] = "done"
        except Exception as e:
            job["status"] = "error"; job["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"job_id": jid}


@app.get("/api/docking/receptor/custom/job/{jid}")
def docking_receptor_custom_job(jid: str):
    job = _CUSTOM_RECEPTOR_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    r = {"status": job["status"], "step": job.get("step")}
    if job["status"] == "done":
        r["profile"] = job["profile"]
    if job["status"] == "error":
        r["error"] = job["error"]
    return r


class AdvancedDocking(BaseModel):
    """Every field is optional and defaults to Automatic — set only the ones
       an expert user actually overrode (see PROJECT_DOCUMENTATION.md's
       'Hide complexity, not evidence' principle: Automatic stays the
       default everywhere, this is opt-in per request, never persisted)."""
    exhaustiveness: Optional[int] = Field(default=None, ge=1, le=64)
    n_poses: Optional[int] = Field(default=None, ge=1, le=20)
    use_gnina: Optional[bool] = None          # None = automatic (use if installed); False = force off
    docking_mode: Optional[str] = None        # None/"site_specific" (default) | "blind"
    box_center: Optional[List[float]] = Field(default=None, min_length=3, max_length=3)
    box_size: Optional[List[float]] = Field(default=None, min_length=3, max_length=3)
    custom_profile: Optional[dict] = None     # from /api/docking/receptor/custom's job result


BLIND_CAVEAT = ("Blind docking: searching the ENTIRE protein surface, not a specific (or even assumed) "
                "pocket. Vina must cover a much larger volume than a site-specific box, which is both "
                "slower and substantially less reliable per-site — treat any hit as a candidate binding "
                "site to investigate, not a confirmed pose or affinity. Consider raising exhaustiveness "
                "well above the default 8 for a real blind search.")


def _resolve_docking_setup(target_id, advanced):
    """Turns a (possibly-empty) AdvancedDocking into (profile, engine,
       rescorer, n_poses, caveat) — shared by /api/docking/submit and the
       Screen pipeline so 'automatic unless overridden' behaves identically
       in both places.

       Raises HTTPException the same way the plain-automatic path already
       did when nothing usable is available."""
    from docking.engines import VinaEngine, GninaRescorer, NullRescorer
    adv = advanced or AdvancedDocking()
    blind = (adv.docking_mode == "blind")
    caveat = None

    if adv.custom_profile:
        profile = dict(adv.custom_profile)
        if profile.get("target_id") != target_id:
            raise HTTPException(400, "custom_profile's target_id doesn't match the request's target_id")
    else:
        try:
            profile = DOCK_PROFILE.load_profile(target_id)
        except Exception as e:
            raise HTTPException(404, str(e))

    if blind:
        from docking import receptor_prep as RP
        receptor_pdb = profile.get("receptor_pdb")
        if not receptor_pdb or not os.path.exists(receptor_pdb):
            raise HTTPException(404, f"no prepared receptor on disk for '{target_id}' — blind docking needs "
                                     "at least a stripped/repaired receptor.")
        try:
            center, size = RP.box_from_receptor(receptor_pdb)
        except Exception as e:
            raise HTTPException(500, f"could not compute a blind (whole-protein) box: {e}")
        profile = dict(profile, center=center, box_size=size, site_source="blind_whole_protein")
        caveat = (caveat + " " if caveat else "") + BLIND_CAVEAT
    elif adv.box_center or adv.box_size:
        profile = dict(profile)
        if adv.box_center: profile["center"] = adv.box_center
        if adv.box_size: profile["box_size"] = adv.box_size

    # Basic well-formedness check, not a scientific gate: a target with no
    # real co-crystallized ligand (site_source == "none_validated", see
    # scripts/audit_ligand_selection.py) has no center/box_size at all —
    # fail with a clear message instead of handing Vina None coordinates.
    if not blind and (not profile.get("center") or not profile.get("box_size")):
        raise HTTPException(400, f"'{target_id}' has no site-specific binding box on record — "
                                 "use Blind mode, or Advanced Settings to pick a structure with a real ligand.")

    engine = VinaEngine(exhaustiveness=adv.exhaustiveness or 8)
    rescorer = NullRescorer() if adv.use_gnina is False else GninaRescorer()
    n_poses = adv.n_poses or 9
    return profile, engine, rescorer, n_poses, caveat


class DockBody(BaseModel):
    target_id: str
    smiles: List[str]
    advanced: Optional[AdvancedDocking] = None
    # A1 — optional project-intake metadata (plant/organism this batch of
    # compounds came from), carried through to run_metadata/the export
    # package so every result stays traceable back to its natural source.
    plant_source: Optional[str] = None


def _raw_pdb_path_for_profile(target_id, profile):
    """Path to the RAW (un-stripped, still has every co-crystallized
       ligand) structure a profile was built from — needed for 'redock
       against a different ligand in the same structure' (see
       docking_alternate_ligands/docking_redock_alternate below), since
       receptor_pdb (the cleaned one) has had every heteroatom removed.
       A manually-picked Advanced Settings structure carries its own
       raw_pdb_path (see docking_receptor_custom); a registry-default
       target's raw file always lives at DOCKING_TARGETS_DIR/<target_id>/
       <pdb_source> (see scripts/validate_target.py, which fetches it
       there before receptor_prep.build_receptor ever runs)."""
    if profile.get("raw_pdb_path"):
        return profile["raw_pdb_path"]
    if profile.get("pdb_source"):
        return os.path.join(DOCK_PROFILE.DOCKING_TARGETS_DIR, target_id, profile["pdb_source"])
    return None


def _build_run_metadata(target_id, profile, engine, rescorer, n_poses, smiles, plant_source=None):
    """A6 — reproducibility snapshot for one docking/screen run: everything
       needed to redo it later (or write a Methods section from it) that
       isn't already implicit in the code itself. Captured at SUBMIT time
       (not read back from a job dict after completion, which nulls the
       profile out) — see docking_submit's own comment on the same issue
       for receptor_pdb_path."""
    import datetime
    from versions import snapshot as version_snapshot
    return {
        "target_id": target_id,
        "plant_source": plant_source,  # A1 — optional, traces this batch back to its natural source
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_compounds_submitted": len(smiles),
        "docking_mode": "blind" if profile.get("site_source") == "blind_whole_protein" else "site_specific",
        "box_center": profile.get("center"),
        "box_size": profile.get("box_size"),
        "pdb_source": profile.get("pdb_source"),
        "raw_pdb_path": _raw_pdb_path_for_profile(target_id, profile),
        "reference_ligand_resname": profile.get("reference_ligand_resname"),
        "engine": getattr(engine, "name", None),
        "exhaustiveness": getattr(engine, "exhaustiveness", None),
        "n_poses": n_poses,
        "rescorer": type(rescorer).__name__ if rescorer is not None else None,
        "ligand_embedding_seeds": [0xf00d, 1, 42, 7],  # fixed constant — see docking/ligand_prep.py's prepare_ligand()
        "software_versions": version_snapshot(),
    }


def _submit_docking(target_id, smiles_list, advanced, plant_source=None):
    """Shared by /api/docking/submit (a normal user submission) and
       /api/docking/job/{jid}/reproduce (A6 — resubmits with a past job's
       EXACT saved parameters instead of whatever's currently the UI/
       registry default) — identical job-creation logic either way."""
    import uuid
    if DOCK_AVAIL is None or not DOCK_AVAIL.status()["ready"]:
        raise HTTPException(503, "Docking is not available — install Vina and prep a receptor. See the Docking tab.")
    profile, engine, rescorer, n_poses, caveat = _resolve_docking_setup(target_id, advanced)
    smiles = [s.strip() for s in smiles_list if s and s.strip()]
    if not smiles:
        raise HTTPException(400, "No SMILES provided.")
    jid = uuid.uuid4().hex[:12]
    run_metadata = _build_run_metadata(target_id, profile, engine, rescorer, n_poses, smiles, plant_source=plant_source)
    # captured now (not read back off job["profile"], which _run_docking_job
    # nulls out once done) — the 3D pose viewer needs the receptor that was
    # ACTUALLY docked against, which for an Advanced Settings custom_profile
    # (any target, not just no-QSAR-model ones) is NOT the same file
    # /api/docking/receptor/{target_id} would serve (that's always the
    # registry's automatic default, wrong receptor for a manual override —
    # and a plain KeyError/500 for a GENE_ target_id, which has no registry
    # entry at all).
    receptor_pdb_path = profile.get("receptor_pdb")
    pdb_source = profile.get("pdb_source")
    _DOCK_JOBS[jid] = {"status": "queued", "total": len(smiles), "done": 0, "results": [], "caveat": caveat,
                       "profile": profile, "smiles": smiles, "engine": engine, "rescorer": rescorer, "n_poses": n_poses,
                       "receptor_pdb_path": receptor_pdb_path, "pdb_source": pdb_source, "run_metadata": run_metadata,
                       # Fresh Decoy Validation runs as its OWN job (see
                       # enrichment_fresh_submit) — keyed here by compound
                       # SMILES once it finishes, so the export package can
                       # include it instead of it only ever being shown
                       # inline and then discarded.
                       "fresh_decoy_results": {}}
    _run_docking_job(jid)      # background thread
    return {"job_id": jid, "total": len(smiles), "caveat": caveat, "pdb_source": pdb_source}


@app.post("/api/docking/submit")
def docking_submit(body: DockBody):
    return _submit_docking(body.target_id, body.smiles, body.advanced, plant_source=body.plant_source)


@app.post("/api/docking/job/{jid}/reproduce")
def docking_reproduce(jid: str):
    """A6 — 'Reproduce this analysis': resubmits a past job with its EXACT
       saved parameters (receptor FILE, box, exhaustiveness, GNINA on/off,
       compound list) rather than whatever the registry/UI currently
       defaults to. This matters because a target's default receptor/box
       CAN change later (re-validation, a newer structure) — pinning the
       actual file path recorded in run_metadata, via custom_profile,
       guarantees byte-identical inputs regardless of what's changed
       since, instead of silently drifting if we re-resolved 'automatic'
       from today's registry. Only works while the original job is still
       in memory (same lifetime limit /export_package already has)."""
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id")
    if not target_id:
        raise HTTPException(409, "this job has no reproducibility metadata to reproduce from")
    receptor_pdb_path = job.get("receptor_pdb_path")
    if not receptor_pdb_path or not os.path.exists(receptor_pdb_path):
        raise HTTPException(409, "the exact receptor file this job used is no longer on disk — cannot reproduce exactly")
    smiles = [r.get("smiles") for r in (job.get("results") or []) if r.get("smiles")]
    if not smiles:
        raise HTTPException(409, "no compound list recorded for this job")

    custom_profile = {
        "target_id": target_id,
        "receptor_pdb": receptor_pdb_path,
        "center": rm.get("box_center"),
        "box_size": rm.get("box_size"),
        "pdb_source": rm.get("pdb_source"),
    }
    advanced = AdvancedDocking(
        exhaustiveness=rm.get("exhaustiveness"),
        n_poses=rm.get("n_poses"),
        use_gnina=(False if rm.get("rescorer") == "NullRescorer" else None),
        docking_mode=rm.get("docking_mode"),
        custom_profile=custom_profile,
    )
    return _submit_docking(target_id, smiles, advanced, plant_source=rm.get("plant_source"))


def _alternate_ligands_response(raw_pdb_path, reference_ligand_resname, box_center):
    """Shared by docking_alternate_ligands/screen_alternate_ligands — the
       SAME raw PDB structure a job docked against often has more than one
       real co-crystallized ligand (different binding sites, or one copy
       per chain in a crystallographic dimer); the pipeline always
       silently picks just the single largest one as THE reference. This
       lists every real alternative (see docking/receptor_prep.py's
       list_ligands()) so the UI can offer 'dock again centered on a
       different one, same settings otherwise' beside 'Reproduce this
       analysis'."""
    if not raw_pdb_path or not os.path.exists(raw_pdb_path):
        return {"available": False, "note": "The raw structure this job used is no longer on disk.", "current": None, "ligands": []}
    from docking import receptor_prep as RP
    try:
        ligands = RP.list_ligands(raw_pdb_path)
    except Exception as e:
        return {"available": False, "note": f"Could not read ligands from the raw structure: {e}", "current": None, "ligands": []}
    current = None
    if reference_ligand_resname:
        same_name = [l for l in ligands if l["resname"] == reference_ligand_resname]
        if len(same_name) > 1 and box_center:
            bc = box_center
            current = min(same_name, key=lambda l: sum((a - b) ** 2 for a, b in zip(l["center"], bc)))
        elif same_name:
            current = same_name[0]
    return {"available": True, "current": current, "ligands": ligands}


@app.get("/api/docking/job/{jid}/alternate_ligands")
def docking_alternate_ligands(jid: str):
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    return _alternate_ligands_response(rm.get("raw_pdb_path"), rm.get("reference_ligand_resname"), rm.get("box_center"))


class AlternateLigandBody(BaseModel):
    resname: str
    chain: str
    resnum: int


@app.post("/api/docking/job/{jid}/alternate_ligand/build")
def docking_alternate_ligand_build(jid: str, body: AlternateLigandBody):
    """Step 1 of 'dock again against a different co-crystallized ligand in
       the same PDB structure' (the button beside 'Reproduce this
       analysis'): rebuilds a receptor centered on body's ligand INSTEAD
       of whichever one the original run used — reusing the SAME raw
       structure file already on disk (no re-fetch, unlike
       docking_receptor_custom's from-scratch flow, since this job's
       run_metadata already has it). Async job reusing _CUSTOM_RECEPTOR_JOBS — poll
       via the EXISTING /api/docking/receptor/custom/job/{jid}, so no new
       polling endpoint/frontend code is needed. Step 2
       (docking_alternate_ligand_submit) resubmits docking, with every
       OTHER setting pinned to the original run, once this profile is
       ready."""
    import threading, uuid
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id")
    raw_pdb_path = rm.get("raw_pdb_path")
    if not target_id or not raw_pdb_path or not os.path.exists(raw_pdb_path):
        raise HTTPException(409, "the raw structure this job used is no longer on disk")

    bjid = uuid.uuid4().hex[:12]
    _CUSTOM_RECEPTOR_JOBS[bjid] = {"status": "queued", "step": "Queued", "profile": None, "error": None}

    def work():
        j = _CUSTOM_RECEPTOR_JOBS[bjid]
        j["status"] = "running"
        def _step(label):
            j["step"] = label
        try:
            out_dir = os.path.join("docking_targets", "_custom")
            work_id = f"{target_id}__altlig_{uuid.uuid4().hex[:8]}"
            os.makedirs(os.path.join(out_dir, work_id), exist_ok=True)
            profile = _build_receptor_profile(raw_pdb_path, work_id, out_dir, target_id,
                                              f"{target_id} (alternate ligand {body.resname})",
                                              body.resname, body.chain, resnum=body.resnum, progress=_step)
            j["step"] = "Done"
            j["profile"] = profile
            j["status"] = "done"
        except Exception as e:
            j["status"] = "error"; j["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"job_id": bjid}


class AlternateLigandSubmitBody(BaseModel):
    profile: dict   # a ReceptorProfile from docking_alternate_ligand_build's job result


@app.post("/api/docking/job/{jid}/alternate_ligand/submit")
def docking_alternate_ligand_submit(jid: str, body: AlternateLigandSubmitBody):
    """Step 2: resubmits docking with the newly-built alternate-ligand
       profile as custom_profile, but every OTHER setting (exhaustiveness,
       poses, GNINA, docking mode, compound list, plant source) pinned to
       the ORIGINAL job — same 'keep settings identical' contract as
       docking_reproduce, just with a substituted receptor/box instead of
       the same one."""
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id")
    if not target_id:
        raise HTTPException(409, "this job has no reproducibility metadata to reproduce from")
    smiles = [r.get("smiles") for r in (job.get("results") or []) if r.get("smiles")]
    if not smiles:
        raise HTTPException(409, "no compound list recorded for this job")

    advanced = AdvancedDocking(
        exhaustiveness=rm.get("exhaustiveness"),
        n_poses=rm.get("n_poses"),
        use_gnina=(False if rm.get("rescorer") == "NullRescorer" else None),
        docking_mode=rm.get("docking_mode"),
        custom_profile=body.profile,
    )
    return _submit_docking(target_id, smiles, advanced, plant_source=rm.get("plant_source"))


def _run_docking_job(jid):
    import threading
    def work():
        job = _DOCK_JOBS[jid]
        job["status"] = "running"
        try:
            for s in job["smiles"]:
                _check_cancelled(job)
                result = DOCK_PIPE.dock_compound(
                    job["profile"], s, engine=job["engine"], rescorer=job["rescorer"],
                    n_poses=job["n_poses"], make_diagram=True)
                job["results"].append(result)
                job["done"] += 1
            # Redocking validation runs automatically here, using the SAME
            # structure/exhaustiveness/n_poses just used for the compounds
            # above — no separate action. Never fails the compound-docking
            # results if it errors (informational, not a gate).
            try:
                job["redocking_validation"] = DOCK_PIPE.redock_reference_for_profile(
                    job["profile"], exhaustiveness=job["engine"].exhaustiveness, n_poses=job["n_poses"])
            except Exception as e:
                job["redocking_validation"] = {"status": "error", "error": str(e)}
            job["status"] = "done"; job["smiles"] = None; job["profile"] = None
            job["engine"] = None; job["rescorer"] = None
        except JobCancelled:
            job["status"] = "cancelled"; job["smiles"] = None; job["profile"] = None
            job["engine"] = None; job["rescorer"] = None
        except Exception as e:
            job["status"] = "error"; job["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()


@app.post("/api/docking/cancel/{jid}")
def docking_cancel(jid: str):
    job = _DOCK_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    job["cancel_requested"] = True
    return {"ok": True}


@app.get("/api/docking/job/{jid}")
def docking_job(jid: str):
    job = _DOCK_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    r = {"status": job["status"], "done": job["done"], "total": job["total"], "caveat": job.get("caveat"),
         "pdb_source": job.get("pdb_source")}
    if job["status"] in ("done", "cancelled"):
        # Cancelled still returns whatever compounds finished docking
        # before the stop request landed, rather than throwing that work
        # away — the per-compound loop only checks cancel_requested
        # BETWEEN compounds, so job["results"] already has real ones.
        r["results"] = job["results"]
        r["receptor_pdb_path"] = job.get("receptor_pdb_path")
        r["redocking_validation"] = job.get("redocking_validation")
    if job["status"] == "error":
        r["error"] = job.get("error")
    return r


@app.get("/api/docking/job/{jid}/export_package")
def docking_export_package(jid: str):
    """B13's 'one-click experiment package' — metadata.json (A6's
       reproducibility snapshot), receptor.pdb, results.csv, and a
       poses/+interactions/ folder per compound, as one ZIP."""
    from fastapi.responses import StreamingResponse
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    from export_package import build_zip
    data = build_zip(job.get("run_metadata") or {}, job.get("receptor_pdb_path"), job.get("results") or [],
                     fresh_decoy_results=job.get("fresh_decoy_results"), redocking_validation=job.get("redocking_validation"))
    return StreamingResponse(iter([data]), media_type="application/zip",
                             headers={"Content-Disposition": f'attachment; filename="docking_{jid}.zip"'})


def _interaction_diagram_response(smiles, interactions, source, fmt):
    fmt = fmt.lower()
    if fmt not in _PLOT_MEDIA_TYPES:
        raise HTTPException(400, f"unsupported format '{fmt}' — use one of {sorted(_PLOT_MEDIA_TYPES)}")
    if not interactions:
        raise HTTPException(404, "no interaction data for this compound")
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(400, "could not parse this SMILES")
    from docking.interaction_diagram import diagram_bytes
    data = diagram_bytes(mol, interactions, fmt=fmt, title=smiles[:30], source=source or "")
    from fastapi.responses import StreamingResponse
    return StreamingResponse(iter([data]), media_type=_PLOT_MEDIA_TYPES[fmt],
                             headers={"Content-Disposition": f'attachment; filename="interaction_diagram.{fmt}"'})


@app.get("/api/docking/job/{jid}/interaction_diagram")
def docking_interaction_diagram(jid: str, smiles: str, fmt: str = "svg"):
    """B13 — vector/high-res export of one compound's 2D interaction
       diagram, regenerated on demand (from this job's already-stored
       smiles + detected interactions, re-parsed with RDKit — the same
       flat 2D depiction diagram_png already draws, just any matplotlib
       format instead of always base64 PNG) rather than kept per-format
       in job JSON."""
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    dock_row = next((r for r in (job.get("results") or []) if r.get("smiles") == smiles), None)
    if dock_row is None:
        raise HTTPException(404, "no docking result for this SMILES in this job")
    return _interaction_diagram_response(smiles, dock_row.get("interactions"), dock_row.get("interaction_source"), fmt)


def _failure_log_csv(jid, rows):
    import csv, io
    from fastapi.responses import StreamingResponse
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["smiles", "status", "category", "reason", "suggested_action"])
    n = 0
    for r in rows:
        if r.get("status") == "ok":
            continue
        n += 1
        w.writerow([r.get("smiles"), r.get("status"), r.get("category"), r.get("reason") or r.get("error"), r.get("suggested_action")])
    if n == 0:
        w.writerow(["(no failed compounds in this job)", "", "", "", ""])
    return StreamingResponse(iter([buf.getvalue()]), media_type="text/csv",
                             headers={"Content-Disposition": f'attachment; filename="failure_log_{jid}.csv"'})


@app.get("/api/docking/job/{jid}/failure_log")
def docking_failure_log(jid: str):
    """B13 — standalone failure log for a completed batch, without
       downloading the whole experiment package: every non-'ok' compound
       with its failure category/reason/suggested_action (docking/
       failure_diagnostics.py) as one small CSV."""
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    return _failure_log_csv(jid, job.get("results") or [])


class FreshDecoyBody(BaseModel):
    target_id: str
    smiles: str
    advanced: Optional[AdvancedDocking] = None
    n_decoys: int = Field(default=50, ge=5, le=200)
    # Optional linkage back to the docking/screen job this compound came
    # from — when given, the finished result is ALSO stored on that job
    # (keyed by this SMILES) so /export_package can include it instead of
    # it only ever being shown inline and then discarded once the button's
    # own component state resets.
    parent_kind: Optional[str] = None    # "docking" | "screen"
    parent_job_id: Optional[str] = None


_ENRICHMENT_JOBS = {}


@app.post("/api/docking/enrichment/fresh")
def enrichment_fresh_submit(body: FreshDecoyBody):
    """The EXPENSIVE, on-demand tier (see docking/enrichment.py's module
       docstring): ~(n_decoys + 1) fresh Vina runs against NEW decoys
       property-matched to THIS specific compound, using the exact same
       profile (receptor/box/exhaustiveness) the caller would submit for a
       normal dock — pass the same `advanced` block a prior /api/docking/
       submit call for this compound used, so the comparison is apples-to-
       apples with whatever structure/mode that used (including a manually-
       picked Advanced Settings structure, which the free per-compound
       annotation on /api/docking/submit deliberately skips)."""
    import uuid, threading
    if DOCK_AVAIL is None or not DOCK_AVAIL.status()["ready"]:
        raise HTTPException(503, "Docking is not available — install Vina and prep a receptor. See the Docking tab.")
    if not body.smiles.strip():
        raise HTTPException(400, "No SMILES provided.")
    profile, engine, rescorer, n_poses, caveat = _resolve_docking_setup(body.target_id, body.advanced)
    jid = uuid.uuid4().hex[:12]
    _ENRICHMENT_JOBS[jid] = {"status": "queued", "result": None, "error": None, "done": 0, "total": body.n_decoys + 1}

    def work():
        job = _ENRICHMENT_JOBS[jid]
        job["status"] = "running"
        try:
            from docking.enrichment import fresh_decoy_validation

            def on_progress(done, total):
                _check_cancelled(job)
                job["done"] = done
                job["total"] = total

            job["result"] = fresh_decoy_validation(body.target_id, body.smiles, profile,
                                                    engine=engine, n_decoys=body.n_decoys, progress_cb=on_progress)
            job["status"] = "done"
            if body.parent_kind and body.parent_job_id and not job["result"].get("error"):
                parent_jobs = _DOCK_JOBS if body.parent_kind == "docking" else _SCREEN_JOBS if body.parent_kind == "screen" else None
                parent = parent_jobs.get(body.parent_job_id) if parent_jobs is not None else None
                if parent is not None:
                    parent.setdefault("fresh_decoy_results", {})[body.smiles] = job["result"]
        except JobCancelled:
            job["status"] = "cancelled"
        except Exception as e:
            job["status"] = "error"; job["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"job_id": jid}


@app.post("/api/docking/enrichment/fresh/cancel/{jid}")
def enrichment_fresh_cancel(jid: str):
    job = _ENRICHMENT_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    job["cancel_requested"] = True
    return {"ok": True}


@app.get("/api/docking/enrichment/fresh/job/{jid}")
def enrichment_fresh_job(jid: str):
    job = _ENRICHMENT_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    r = {"status": job["status"], "done": job.get("done", 0), "total": job.get("total")}
    if job["status"] == "done":
        r["result"] = job["result"]
    if job["status"] == "error":
        r["error"] = job.get("error")
    return r


# ---------------- screen (the full STEP 1-8 pipeline, CLAUDE.md §7.1) ----------------
from serving import screen as SCREEN

_SCREEN_JOBS = {}


class ScreenBody(BaseModel):
    target_id: str
    smiles: List[str]
    advanced: Optional[AdvancedDocking] = None
    plant_source: Optional[str] = None


def _submit_screen(target_id, smiles_list, advanced, plant_source=None):
    """Shared by /api/screen/submit and /api/screen/job/{jid}/reproduce
       (A6) — see _submit_docking's docstring for why reproduce pins an
       explicit custom_profile rather than re-resolving 'automatic'."""
    import uuid
    smiles = [s.strip() for s in smiles_list if s and s.strip()]
    if not smiles:
        raise HTTPException(400, "No SMILES provided.")
    if target_id not in MA.list_target_ids():
        raise HTTPException(404, f"Unknown target '{target_id}'")
    jid = uuid.uuid4().hex[:12]
    _SCREEN_JOBS[jid] = {"status": "queued", "step": 0, "step_label": "Queued",
                         "done": None, "total": None, "result": None, "error": None,
                         "fresh_decoy_results": {}}
    _run_screen_job(jid, target_id, smiles, advanced.model_dump() if advanced else None, plant_source=plant_source)
    return {"job_id": jid}


@app.post("/api/screen/submit")
def screen_submit(body: ScreenBody):
    return _submit_screen(body.target_id, body.smiles, body.advanced, plant_source=body.plant_source)


@app.post("/api/screen/job/{jid}/reproduce")
def screen_reproduce(jid: str):
    """A6 — same 'reproduce with exact saved parameters' as
       docking_reproduce, adapted for the Screen pipeline's job/result
       shape (compounds come from result['shortlist'], the receptor path
       from result['receptor_pdb_path'])."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    result = job.get("result") or {}
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id") or result.get("target_id")
    if not target_id:
        raise HTTPException(409, "this job has no reproducibility metadata to reproduce from")
    smiles = [r.get("input_smiles") or r.get("smiles") for r in (result.get("shortlist") or []) if r.get("smiles")]
    if not smiles:
        raise HTTPException(409, "no compound list recorded for this job")

    advanced = None
    receptor_pdb_path = result.get("receptor_pdb_path")
    if receptor_pdb_path and os.path.exists(receptor_pdb_path):
        # Docking ran and we know exactly which receptor file it used —
        # pin it, same as docking_reproduce. If docking DIDN'T run in the
        # original job (QSAR-only), advanced stays None and reproduction
        # is QSAR-only again too — faithful either way.
        custom_profile = {
            "target_id": target_id,
            "receptor_pdb": receptor_pdb_path,
            "center": rm.get("box_center"),
            "box_size": rm.get("box_size"),
            "pdb_source": rm.get("pdb_source"),
        }
        advanced = AdvancedDocking(
            exhaustiveness=rm.get("exhaustiveness"),
            n_poses=rm.get("n_poses"),
            docking_mode=rm.get("docking_mode"),
            custom_profile=custom_profile,
        )
    return _submit_screen(target_id, smiles, advanced, plant_source=rm.get("plant_source"))


@app.get("/api/screen/job/{jid}/alternate_ligands")
def screen_alternate_ligands(jid: str):
    """Same as docking_alternate_ligands, adapted for the Screen
       pipeline's job/run_metadata shape."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    return _alternate_ligands_response(rm.get("raw_pdb_path"), rm.get("reference_ligand_resname"), rm.get("box_center"))


@app.post("/api/screen/job/{jid}/alternate_ligand/build")
def screen_alternate_ligand_build(jid: str, body: AlternateLigandBody):
    """Same as docking_alternate_ligand_build, adapted for the Screen
       pipeline's job/run_metadata shape. Still reuses
       _CUSTOM_RECEPTOR_JOBS/its existing poll endpoint — the build step
       itself doesn't care which tab originally ran the job."""
    import threading, uuid
    if DOCK_PROFILE is None:
        raise HTTPException(503, "Docking package not available")
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id") or (job.get("result") or {}).get("target_id")
    raw_pdb_path = rm.get("raw_pdb_path")
    if not target_id or not raw_pdb_path or not os.path.exists(raw_pdb_path):
        raise HTTPException(409, "the raw structure this job used is no longer on disk")

    bjid = uuid.uuid4().hex[:12]
    _CUSTOM_RECEPTOR_JOBS[bjid] = {"status": "queued", "step": "Queued", "profile": None, "error": None}

    def work():
        j = _CUSTOM_RECEPTOR_JOBS[bjid]
        j["status"] = "running"
        def _step(label):
            j["step"] = label
        try:
            out_dir = os.path.join("docking_targets", "_custom")
            work_id = f"{target_id}__altlig_{uuid.uuid4().hex[:8]}"
            os.makedirs(os.path.join(out_dir, work_id), exist_ok=True)
            profile = _build_receptor_profile(raw_pdb_path, work_id, out_dir, target_id,
                                              f"{target_id} (alternate ligand {body.resname})",
                                              body.resname, body.chain, resnum=body.resnum, progress=_step)
            j["step"] = "Done"
            j["profile"] = profile
            j["status"] = "done"
        except Exception as e:
            j["status"] = "error"; j["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"job_id": bjid}


@app.post("/api/screen/job/{jid}/alternate_ligand/submit")
def screen_alternate_ligand_submit(jid: str, body: AlternateLigandSubmitBody):
    """Step 2 for the Screen tab: resubmits with the newly-built alternate-
       ligand profile as custom_profile, every other setting pinned to the
       original run — same contract as screen_reproduce."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    result = job.get("result") or {}
    rm = job.get("run_metadata") or {}
    target_id = rm.get("target_id") or result.get("target_id")
    if not target_id:
        raise HTTPException(409, "this job has no reproducibility metadata to reproduce from")
    smiles = [r.get("input_smiles") or r.get("smiles") for r in (result.get("shortlist") or []) if r.get("smiles")]
    if not smiles:
        raise HTTPException(409, "no compound list recorded for this job")

    advanced = AdvancedDocking(
        exhaustiveness=rm.get("exhaustiveness"),
        n_poses=rm.get("n_poses"),
        docking_mode=rm.get("docking_mode"),
        custom_profile=body.profile,
    )
    return _submit_screen(target_id, smiles, advanced, plant_source=rm.get("plant_source"))


def _run_screen_job(jid, target_id, smiles, advanced=None, plant_source=None):
    import threading

    def on_progress(step, label, done=None, total=None):
        job = _SCREEN_JOBS[jid]
        _check_cancelled(job)   # step 6 (docking) calls this once per compound —
                                 # the only per-item checkpoint the pipeline has
        job["step"] = step; job["step_label"] = label; job["done"] = done; job["total"] = total

    def work():
        job = _SCREEN_JOBS[jid]
        job["status"] = "running"
        try:
            result = SCREEN.run(target_id, smiles, progress=on_progress, advanced=advanced)
            job["result"] = result
            import datetime
            from versions import snapshot as version_snapshot
            job["run_metadata"] = {
                "target_id": target_id,
                "plant_source": plant_source,
                "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "n_compounds_submitted": len(smiles),
                "docking_mode": result.get("docking_mode"),
                "box_center": result.get("box_center"),
                "box_size": result.get("box_size"),
                "pdb_source": result.get("pdb_source"),
                "raw_pdb_path": result.get("raw_pdb_path"),
                "reference_ligand_resname": result.get("reference_ligand_resname"),
                "exhaustiveness": result.get("exhaustiveness"),
                "n_poses": result.get("n_poses"),
                "ligand_embedding_seeds": [0xf00d, 1, 42, 7],
                "software_versions": version_snapshot(),
            }
            job["status"] = "done"
        except JobCancelled:
            job["status"] = "cancelled"
        except Exception as e:
            job["status"] = "error"; job["error"] = str(e)
    threading.Thread(target=work, daemon=True).start()


@app.post("/api/screen/cancel/{jid}")
def screen_cancel(jid: str):
    job = _SCREEN_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    job["cancel_requested"] = True
    return {"ok": True}


@app.get("/api/screen/job/{jid}")
def screen_job(jid: str):
    job = _SCREEN_JOBS.get(jid)
    if not job:
        raise HTTPException(404, "unknown job")
    r = {"status": job["status"], "step": job["step"], "step_label": job["step_label"],
        "total_steps": 8, "done": job["done"], "total": job["total"]}
    if job["status"] == "done":
        r["result"] = job["result"]
    if job["status"] == "error":
        r["error"] = job.get("error")
    return r


@app.get("/api/screen/job/{jid}/export_package")
def screen_export_package(jid: str):
    """Same 'one-click experiment package' as /api/docking/job/{jid}/
       export_package, adapted for the Screen pipeline's shortlist shape
       (each row's docking result nests under row['docking'] rather than
       being the row itself)."""
    from fastapi.responses import StreamingResponse
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    result = job["result"]
    flat_results = []
    for row in result.get("shortlist", []):
        d = dict(row.get("docking") or {})
        d["smiles"] = row.get("smiles")
        flat_results.append(d)
    from export_package import build_zip
    data = build_zip(job.get("run_metadata") or {}, result.get("receptor_pdb_path"), flat_results,
                     fresh_decoy_results=job.get("fresh_decoy_results"), redocking_validation=result.get("redocking_validation"))
    return StreamingResponse(iter([data]), media_type="application/zip",
                             headers={"Content-Disposition": f'attachment; filename="screen_{jid}.zip"'})


@app.get("/api/screen/job/{jid}/interaction_diagram")
def screen_interaction_diagram(jid: str, smiles: str, fmt: str = "svg"):
    """Same as docking_interaction_diagram, adapted for the Screen
       pipeline's shortlist shape."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    row = next((r for r in (job["result"].get("shortlist") or []) if r.get("smiles") == smiles), None)
    if row is None or not row.get("docking"):
        raise HTTPException(404, "no docking result for this SMILES in this job")
    d = row["docking"]
    return _interaction_diagram_response(smiles, d.get("interactions"), None, fmt)


@app.get("/api/screen/job/{jid}/failure_log")
def screen_failure_log(jid: str):
    """Same as docking_failure_log, adapted for the Screen pipeline's
       shortlist shape (each row's docking result nests under
       row['docking'])."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    flat = []
    for row in (job["result"].get("shortlist") or []):
        if row.get("docking") is None:
            continue   # docking never ran for this compound (QSAR-only) — not a failure to log
        d = dict(row["docking"])
        d["smiles"] = row.get("smiles")
        flat.append(d)
    return _failure_log_csv(jid, flat)


@app.get("/api/screen/job/{jid}/export.csv")
def screen_export_csv(jid: str):
    from fastapi.responses import StreamingResponse
    import io
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    rows = job["result"]["shortlist"]
    df = pd.DataFrame([{
        "rank": r["rank"], "input_smiles": r["input_smiles"], "smiles": r["smiles"],
        "predicted_pIC50": r["qsar"]["predicted_pIC50"], "in_domain": r["qsar"]["in_domain"],
        "qsar_confidence": r["qsar"]["confidence"],
        "vina_score": (r["docking"] or {}).get("vina_score"),
        "docking_confidence": (r["docking"] or {}).get("confidence"),
        "fused_score": r["fused_score"], "caveats": "; ".join(r["caveats"]),
    } for r in rows])
    buf = io.StringIO(); df.to_csv(buf, index=False); buf.seek(0)
    return StreamingResponse(iter([buf.getvalue()]), media_type="text/csv",
                             headers={"Content-Disposition": f'attachment; filename="screen_{jid}.csv"'})




# ============================================================
#  A3 — Natural-product similarity search (similarity.py)
# ============================================================
import similarity as SIM


class SimilaritySearchBody(BaseModel):
    smiles: str
    threshold: float = Field(default=0.4, ge=0.0, le=1.0)
    top_n: int = Field(default=50, ge=1, le=200)


@app.get("/api/similarity/status")
def similarity_status():
    """Whether the COCONUT-derived similarity index has been downloaded
       yet — this is a single shared resource (not per-target), so it's
       its own small download flow rather than downloads.py's per-target
       manifest-driven one."""
    return {"available": SIM.available(), "download_base_url": SIM.DOWNLOAD_BASE_URL}


@app.post("/api/similarity/download")
def similarity_download_start():
    if SIM.available():
        return {"job_id": None, "already_installed": True}
    return {"job_id": SIM.start_download()}


@app.get("/api/similarity/download/progress/{job_id}")
def similarity_download_progress(job_id: str):
    p = SIM.download_progress(job_id)
    if p is None:
        raise HTTPException(404, "unknown job")
    return p


@app.post("/api/similarity/download/cancel/{job_id}")
def similarity_download_cancel(job_id: str):
    if not SIM.cancel_download(job_id):
        raise HTTPException(404, "unknown job")
    return {"ok": True}


@app.post("/api/similarity/search")
def similarity_search(body: SimilaritySearchBody):
    """Morgan/ECFP4 Tanimoto similarity + Murcko scaffold matching + MCS
       (for the top hits) against the COCONUT natural-product index."""
    try:
        return SIM.search(body.smiles, threshold=body.threshold, top_n=body.top_n)
    except FileNotFoundError:
        raise HTTPException(503, "Similarity index not downloaded yet — see the Similarity Search tab.")
    except ValueError as e:
        raise HTTPException(400, str(e))


# ============================================================
#  A2 — Compound -> target prediction (target_fishing.py)
# ============================================================
import target_fishing as TF


class TargetFishingBody(BaseModel):
    smiles: str
    threshold: float = Field(default=0.4, ge=0.0, le=1.0)


@app.get("/api/target_fishing/status")
def target_fishing_status():
    """Whether the precomputed target-fishing fingerprint index is
       present — it's built offline (see scripts/build_target_fishing_index.py)
       and shipped as a static artifact, so this is normally always true;
       kept as an explicit check anyway so the frontend can show a clear
       error instead of a raw 503 if a build ever ships without it."""
    return {"available": TF.available()}


@app.post("/api/target_fishing/search")
def target_fishing_search(body: TargetFishingBody):
    """Ligand-based 'which targets might this compound hit' search:
       Tanimoto similarity to known bioactive compounds from a broad
       ChEMBL pull (see target_fishing.py's module docstring for the exact
       filter/provenance and the evidence_score formula) — the same
       'guilt by association' principle SwissTargetPrediction uses, not a
       calibrated probability."""
    try:
        return TF.search(body.smiles, threshold=body.threshold)
    except FileNotFoundError:
        raise HTTPException(503, "Target-fishing index not available in this build.")
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/target_fishing/suggest_compounds")
def target_fishing_suggest_compounds(q: str = "", limit: int = 8):
    """As-you-type suggestions from the target-fishing index (local, no
       network) — matches a SMILES fragment, a target's ChEMBL id, or a
       target's preferred name."""
    try:
        return {"results": TF.suggest_compounds(q, limit=limit)}
    except FileNotFoundError:
        raise HTTPException(503, "Target-fishing index not available in this build.")


# ============================================================
#  A2v2 — Target Prediction v2 (target_prediction_v2.py)
# ============================================================
# A separate, independently-validated method from A2 above, NOT a drop-in
# replacement wired into the same endpoint — see
# target_prediction_v2/METHODS_AND_VALIDATION.md for the full account.
# Density-adaptive retrieval (pool k=10 weighted votes when >=8 neighbours
# exist at Tanimoto>=0.5, else fall back to best-similarity), validated to
# significantly beat both v1's production ranking and a genuine SEA
# reimplementation on genuinely held-out data (G2 result). Kept as its own
# endpoint family, both methods available side by side, rather than
# silently swapping A2's production behaviour.
import target_prediction_v2 as TPV2


class TargetPredictionV2Body(BaseModel):
    smiles: str
    # Upper bound is above n_targets_indexed (4,882 as of the current v2
    # index) so the frontend can request "all" by passing a large value —
    # not truly unbounded, to keep a single pathological request bounded.
    top_k: int = Field(default=25, ge=1, le=10000)


@app.get("/api/target_prediction_v2/status")
def target_prediction_v2_status():
    """Whether the v2 index is present in this build."""
    return {"available": TPV2.available()}


@app.post("/api/target_prediction_v2/predict")
def target_prediction_v2_predict(body: TargetPredictionV2Body):
    """Density-adaptive ligand-based target prediction (v2) — ranked
       targets with per-prediction evidence (supporting neighbour
       compounds, potency weights, orthologue provenance where used) and
       an explicit, un-calibrated confidence label. See
       target_prediction_v2.py's module docstring for the validated
       method and target_prediction_v2/METHODS_AND_VALIDATION.md for the
       full validation account — this is a research tool's output, not a
       clinical or regulatory-grade probability."""
    try:
        return TPV2.predict(body.smiles, top_k=body.top_k)
    except FileNotFoundError:
        raise HTTPException(503, "Target Prediction v2 index not available in this build.")
    except ValueError as e:
        raise HTTPException(400, str(e))


# ============================================================
#  A4 — Literature intelligence (literature.py)
# ============================================================
import literature as LIT


class LiteratureSearchBody(BaseModel):
    query: str
    max_results: int = Field(default=10, ge=1, le=30)


@app.post("/api/literature/search")
def literature_search(body: LiteratureSearchBody):
    """Live PubMed search (NCBI E-utilities) for a target/compound/plant-
       source name — this app's first request-time external network call,
       so failures are surfaced as 502 (network/API) rather than crashing."""
    try:
        return LIT.search(body.query, max_results=body.max_results)
    except ValueError as e:
        raise HTTPException(400, str(e))
    except LIT.LiteratureError as e:
        raise HTTPException(502, str(e))


# ============================================================
#  A5 — Research story generator (research_report.py)
# ============================================================
class ResearchReportBody(BaseModel):
    smiles: str
    include_literature: bool = True


@app.post("/api/docking/job/{jid}/research_report")
def docking_research_report(jid: str, body: ResearchReportBody):
    """Assembles the A5 evidence-chain report for ONE compound from an
       already-completed Docking-tab job — reuses A1/A2/A4's data, never
       re-runs docking. See research_report.py's module docstring."""
    job = _DOCK_JOBS.get(jid)
    if not job or job["status"] not in ("done", "cancelled"):
        raise HTTPException(404, "job not found or not finished")
    dock_row = next((r for r in (job.get("results") or []) if r.get("smiles") == body.smiles), None)
    if dock_row is None:
        raise HTTPException(404, "no docking result for this SMILES in this job")
    rm = job.get("run_metadata") or {}
    import research_report as RR
    report = RR.build_report(rm.get("target_id"), body.smiles, dock_row, run_metadata=rm,
                              plant_source=rm.get("plant_source"), include_literature=body.include_literature)
    return {"report": report, "markdown": RR.to_markdown(report)}


@app.post("/api/screen/job/{jid}/research_report")
def screen_research_report(jid: str, body: ResearchReportBody):
    """Same as docking_research_report, adapted for the Screen pipeline's
       shortlist shape (each row's docking result nests under
       row['docking'])."""
    job = _SCREEN_JOBS.get(jid)
    if not job or job["status"] != "done":
        raise HTTPException(404, "job not found or not finished")
    result = job["result"]
    row = next((r for r in (result.get("shortlist") or []) if r.get("smiles") == body.smiles), None)
    if row is None:
        raise HTTPException(404, "no result for this SMILES in this job")
    dock_row = dict(row.get("docking") or {})
    dock_row["smiles"] = body.smiles
    rm = job.get("run_metadata") or {}
    import research_report as RR
    report = RR.build_report(rm.get("target_id") or result.get("target_id"), body.smiles, dock_row, run_metadata=rm,
                              plant_source=rm.get("plant_source"), include_literature=body.include_literature)
    return {"report": report, "markdown": RR.to_markdown(report)}
