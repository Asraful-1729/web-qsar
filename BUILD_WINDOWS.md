# PhytoScreen Desktop — Windows build (.exe)

A native window (PyWebview) wrapping the existing FastAPI JSON API
(`backend/app.py`, **unchanged — not one line edited by this conversion**)
plus the built frontend (`frontend/dist/`), served from one process on one
localhost port — see `desktop.py`. Adapted from the sibling `qsar-desktop`
project's own proven Windows build (same architecture, already shipping),
repointed at this project's current code and data layout.

**Everything in this document is new, additive infrastructure — `backend/`,
`frontend/src/`, and every existing route/component are byte-for-byte
unchanged.** `desktop.py` and `desktop_remote_gnina.py` are new top-level
files; `backend/remote_gnina_service.py` is a new file too, never imported
by `backend/app.py`.

## GNINA on Windows — read this first

**There is no official GNINA build for Windows.** Checked every GitHub
release (v1.0 through v1.3.3) — every asset is a Linux ELF binary.
`qsar-desktop`'s own prior build effort independently found and documented
the same thing ("fpocket and GNINA have no official Windows builds — both
already degrade gracefully"). So this build cannot bundle `gnina.exe` the
way it bundles `vina.exe` — that file doesn't exist.

Instead: `desktop.py` patches `docking.engines.GninaRescorer` (at the
class-object level, not a source edit — see that file's docstring) so
that when no local `gnina` binary is found, it optionally calls out to
`backend/remote_gnina_service.py` running on a real GPU machine, over
HTTP, instead of just being unavailable. **This is a genuinely separate
piece of infrastructure from the desktop .exe itself** — it runs
continuously on a server (not shipped inside the installer), and the
desktop app just needs its URL + a bearer token at build/run time (see
"Remote GNINA setup" below). If unset or unreachable, GNINA rescoring is
simply unavailable and Vina docking still works exactly as before — GNINA
has always been an optional "second opinion," never required.

fpocket **also has no official Windows build** (same finding as
`qsar-desktop`) — it degrades gracefully the same way (pocket detection
falls back to whatever the Docking tab's own fallback path already does
without it).

## 1. Install (on the Windows machine, in a venv)
    python -m venv .venv
    .venv\Scripts\activate
    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install -r backend\requirements.txt
    pip install -r requirements-desktop.txt

    cd frontend
    npm install
    npm run build
    cd ..

    # optional, for gold-standard interaction typing (falls back to a
    # built-in distance-based detector if absent):
    conda install -c conda-forge plip openbabel

AutoDock Vina **is** a hard requirement for docking to activate at all —
download the official Windows build and put `vina.exe` on `PATH`, or drop
it in `bin\vina.exe` next to `desktop.py` (`desktop.py` prepends `bin\` to
`PATH` at startup).

## 2. Prepare the bundled data (once, before packaging)
Everything below is small enough to ship directly in the installer
(unlike `models/`/`docking_targets/`, ~101GB/2.3GB — those stay
on-demand, see §4):

    docking_registry.json          (~2MB, metadata only)
    panel_results_v2.csv           (~37MB)
    models\curated\                (~15MB)
    target_prediction_v2_data\     (~241MB — v2_index, orthologue_index,
                                     target_name_map.json; already staged
                                     at the repo root, built from
                                     target_prediction_v2/phase2/v2_index,
                                     target_prediction_v2/phase3/orthologue_index,
                                     and target_prediction_v2/phase2/data/target_name_map.json —
                                     re-copy after any v2 index rebuild)

## 3. Remote GNINA setup (once, on a GPU server — not the Windows machine)
    cd backend
    export REMOTE_GNINA_TOKEN=$(python3 -c "import secrets; print(secrets.token_urlsafe(32))")
    uvicorn remote_gnina_service:app --host 0.0.0.0 --port 8300

Requires a working local `gnina` + matching `libcudnn` on THAT machine
(see the Docking tab's own status check, or `docking/availability.py`) —
this is infrastructure for the server side, unrelated to anything the
Windows build needs installed. Put this behind a real reverse proxy with
TLS before exposing it beyond a trusted network — the bearer token
travels in the `Authorization` header, in cleartext over plain HTTP, so a
plain `0.0.0.0:8300` exposure is only appropriate on a network you trust
(see "Security note" at the end of this document).

Then, when building/running the desktop app:

    set REMOTE_GNINA_URL=https://your-gpu-server:8300
    set REMOTE_GNINA_TOKEN=<the token generated above>

(PyInstaller doesn't bake env vars into the .exe — set these as real
Windows environment variables on the machine running the built app, or
bake sensible defaults into `desktop_remote_gnina.py`'s
`REMOTE_URL`/`REMOTE_TOKEN` before building, if every install should
point at the same server by default.)

## 4. Data: on-demand downloads (Cloudflare R2)
Unchanged from `qsar-desktop` — defaults to the project's public bucket,
nothing to set for a normal build:

    set DOWNLOAD_BASE_URL=https://<your-bucket-or-cdn>/   (only if forking)

`models/` and `docking_targets/` start out empty and are populated per-target
by the Downloads tab, exactly as documented in `backend/downloads.py`.

## 5. Run in dev first (confirm it works before packaging)
    python desktop.py
Validated already (this session, on Linux, headless): backend starts,
`/api/health` responds, the full frontend serves correctly at `/`
(pixel-identical to the plain web deployment — confirmed via screenshot),
`target_prediction_v2` and `docking` status both report correctly, and
pywebview falls back to printing a URL when no GTK/QT is available rather
than crashing. **Not yet validated**: actual native window rendering
(needs real Windows — this machine has no GUI), the PyInstaller `.exe`
itself (PyInstaller doesn't cross-compile — must be built ON Windows),
and a real Vina/GNINA run end-to-end inside a packaged build.

## 6. Build the .exe (PyInstaller)
From the project root:
    pip install pyinstaller
    pyinstaller --noconfirm --windowed --name PhytoScreen ^
      --paths backend ^
      --add-data "frontend/dist;frontend/dist" ^
      --add-data "docking_registry.json;." ^
      --add-data "panel_results_v2.csv;." ^
      --add-data "models/curated;models/curated" ^
      --add-data "target_prediction_v2_data;target_prediction_v2_data" ^
      --add-data "bin;bin" ^
      --collect-all rdkit ^
      --collect-all autogluon ^
      --collect-all chemprop ^
      --collect-all lightning ^
      --collect-all posebusters ^
      --collect-all meeko ^
      --collect-all gemmi ^
      --collect-all webview ^
      --collect-all admet_ai ^
      --collect-all cuik_molmaker ^
      --collect-all openmm ^
      --collect-all pdbfixer ^
      --collect-all lightgbm ^
      --collect-all catboost ^
      --collect-data xgboost ^
      --hidden-import xgboost ^
      --hidden-import uvicorn ^
      --hidden-import desktop_remote_gnina ^
      desktop.py

The only difference from `qsar-desktop`'s own command: the new
`--add-data "target_prediction_v2_data;..."` line, and
`--hidden-import desktop_remote_gnina` (`desktop.py` imports it by name
inside a function — `import desktop_remote_gnina as remote` — which
PyInstaller's static analysis should already catch since it's a
top-level `import` statement, but it's called out explicitly here since
missing it fails exactly the same silent way `--paths backend` does if
ever skipped: the app starts, logs "GNINA remote-rescoring patch FAILED",
and simply never offers GNINA rescoring, rather than crashing).

The .exe lands in `dist\PhytoScreen\PhytoScreen.exe` (onedir). Every
`--collect-all`/`--hidden-import` reasoning below is inherited unchanged
from `qsar-desktop`'s own `BUILD_WINDOWS.md` (same dependency stack) —
see that document for the full "why" on each one (`cuik_molmaker`,
`openmm`/`pdbfixer`, `lightgbm`/`catboost`/`xgboost`'s native `.dll`
loading, `webview`'s `WebView2Loader.dll`, etc.) if a specific module
turns up missing at runtime.

A packaged, double-clickable installer (Start Menu/Desktop shortcuts) is
built from this output via Inno Setup — see `installer/phytoscreen.iss`
(adapted from `qsar-desktop`'s own script — same structure, same file
list plus the new `target_prediction_v2_data\` directory).

## Security note — the remote GNINA token
`desktop_remote_gnina.py`'s token travels in an HTTP header on every
rescoring call. If `REMOTE_GNINA_URL` is baked into every installed copy
of the app (so it "just works" for a non-technical user with zero
configuration), the token is effectively embedded in a distributed
binary — extractable by anyone with the .exe, same trust model as an API
key shipped in a desktop app. This is an accepted, disclosed trade-off
for a low-stakes research tool with a small, known user base (unlike a
public product), not a hardened design — put the remote service behind
TLS at minimum, and rotate the token if the installer is ever shared
beyond its intended recipients. **Not yet decided**: whether to expose
`remote_gnina_service.py` on the open internet at all, and if so, at what
address — see the project conversation for that open decision.
