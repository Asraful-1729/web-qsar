"""
============================================================
  PhytoScreen DESKTOP  (desktop.py)
============================================================
  Wraps the existing FastAPI JSON API (backend/app.py, UNCHANGED -- not
  one line of it is edited by this conversion) in a native desktop
  window (no browser, no manual localhost, no separate frontend dev
  server) by mounting the built frontend (frontend/dist/) onto the same
  app instance, added last so it never shadows /api/*.

  Adapted from the sibling qsar-desktop project's own proven desktop.py
  (PyWebview + PyInstaller + Inno Setup, already used to ship a working
  Windows build) -- same architecture, repointed at THIS project's
  current code and data layout (which has moved on since qsar-desktop's
  snapshot: target_prediction_v2, the Docking GNINA fix, etc.), not
  rebuilt from scratch.

  GNINA on Windows: there is no official Windows GNINA build (checked
  every GitHub release, v1.0-v1.3.3 -- all Linux-only; qsar-desktop's
  own build effort found the same thing independently), so this cannot
  just bundle a gnina.exe the way it bundles vina.exe. Instead, BEFORE
  the existing backend is imported, this patches
  docking.engines.GninaRescorer's available()/rescore() methods (see
  desktop_remote_gnina.py) so that when no local gnina binary is found,
  it optionally calls out to a remote GNINA service running on a real
  GPU machine (backend/remote_gnina_service.py) instead of simply being
  unavailable. This is a class-level monkeypatch applied to the shared
  GninaRescorer class object itself -- not a source edit to
  docking/engines.py, docking/pipeline.py, or app.py, and it changes
  nothing for the plain web deployment (only desktop.py ever imports
  desktop_remote_gnina). If no REMOTE_GNINA_URL is configured (or it's
  unreachable), rescoring is simply unavailable, exactly like the
  original class's own behaviour with no local binary -- Vina docking
  always still works, GNINA has always been an optional second opinion.

  Run (dev):   python desktop.py        (needs `npm run build` in
                                          frontend/ first, or set
                                          PHYTO_SKIP_STATIC=1 to run
                                          API-only and open the Vite dev
                                          server separately)
  Build .exe:  see BUILD_WINDOWS.md

  --windowed (no console) means print()/stderr go nowhere a user can
  ever see -- a startup failure used to just look like "nothing
  happens." Everything below also writes to a persistent log file
  (%LOCALAPPDATA%\\PhytoScreen\\desktop.log on Windows) and, on a truly
  fatal failure, pops a native message box pointing at it.
============================================================
"""
import os, sys, threading, socket, time, contextlib, traceback, datetime

# A --windowed PyInstaller build has NO console at all on Windows --
# sys.stdout/sys.stderr are None, not just redirected/closed. print()
# specifically special-cases this and silently no-ops, but plenty of
# library code doesn't (e.g. uvicorn's colorized-logging formatter calls
# sys.stdout.isatty() during Config.__init__ and crashes with
# AttributeError on None) -- give every library a real, harmless
# writable stream instead, as early as possible.
if sys.stdout is None or sys.stderr is None:
    _null = open(os.devnull, "w")
    sys.stdout = sys.stdout or _null
    sys.stderr = sys.stderr or _null

HOST = "127.0.0.1"
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
FROZEN_ROOT = getattr(sys, "_MEIPASS", REPO_ROOT)
BACKEND_DIR = os.path.join(REPO_ROOT, "backend")
FRONTEND_DIST = os.path.join(FROZEN_ROOT, "frontend", "dist")
BIN_DIR = os.path.join(FROZEN_ROOT, "bin")


def _log_dir():
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
        return os.path.join(base, "PhytoScreen")
    return REPO_ROOT


LOG_PATH = os.path.join(_log_dir(), "desktop.log")


def _log(msg):
    """Never let logging itself crash the app -- best-effort only."""
    try:
        os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(f"{datetime.datetime.now().isoformat(timespec='seconds')} {msg}\n")
    except Exception:
        pass
    print(msg)


def _fatal(msg):
    _log(f"[desktop] FATAL: {msg}")
    if os.name == "nt":
        try:
            import ctypes
            ctypes.windll.user32.MessageBoxW(
                0,
                f"PhytoScreen failed to start:\n\n{msg}\n\nDetails were written to:\n{LOG_PATH}",
                "PhytoScreen — startup error",
                0x10,  # MB_ICONERROR
            )
        except Exception:
            pass
    sys.exit(1)


def _free_port():
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind((HOST, 0))
        return s.getsockname()[1]


def _wait_up(port, timeout=30):
    import urllib.request
    end = time.time() + timeout
    while time.time() < end:
        try:
            urllib.request.urlopen(f"http://{HOST}:{port}/api/health", timeout=1)
            return True
        except Exception:
            time.sleep(0.25)
    return False


def _put_bundled_binaries_on_path():
    """Docking engines (backend/docking/engines.py, availability.py) find
       vina/fpocket via shutil.which() -- bundling them in bin/ and
       prepending it to PATH here means that code needs zero changes.
       GNINA is deliberately NOT expected here on Windows (see module
       docstring) -- its availability comes from the remote-rescoring
       patch below instead, or from a genuinely local gnina if someone
       runs this build on Linux/macOS with one on PATH already."""
    if os.path.isdir(BIN_DIR):
        os.environ["PATH"] = BIN_DIR + os.pathsep + os.environ.get("PATH", "")


def _patch_gnina_rescorer():
    """Applied BEFORE `import app` (and therefore before anything else
       imports docking.engines/docking.pipeline) so every subsequent
       `from docking.engines import GninaRescorer` -- wherever it
       happens, module-level or lazy -- resolves to the SAME class
       object this function mutates. See desktop_remote_gnina.py's own
       docstring for the full reasoning."""
    if BACKEND_DIR not in sys.path:
        sys.path.insert(0, BACKEND_DIR)
    try:
        import docking.engines as engines
        import desktop_remote_gnina as remote

        original_rescore = engines.GninaRescorer.rescore
        engines.GninaRescorer.available = remote.patched_available
        engines.GninaRescorer.rescore = remote.make_patched_rescore(original_rescore)
        _log(f"[desktop] GNINA remote-rescoring patch applied "
             f"(REMOTE_GNINA_URL={'set' if remote._remote_configured() else 'not set'})")
    except Exception:
        _log(f"[desktop] GNINA remote-rescoring patch FAILED (non-fatal -- "
             f"GNINA rescoring simply won't be offered):\n{traceback.format_exc()}")


def build_app():
    """Import the existing serving app (backend/app.py, unchanged) and,
       unless skipped, mount the built frontend onto it."""
    if BACKEND_DIR not in sys.path:
        sys.path.insert(0, BACKEND_DIR)
    import app as serving          # backend/app.py (unchanged)

    if not os.environ.get("PHYTO_SKIP_STATIC"):
        if os.path.isdir(FRONTEND_DIST):
            from fastapi.staticfiles import StaticFiles
            serving.app.mount(
                "/", StaticFiles(directory=FRONTEND_DIST, html=True), name="frontend"
            )
        else:
            _log(f"[desktop] no built frontend at {FRONTEND_DIST} "
                 f"(run `npm run build` in frontend/, or set PHYTO_SKIP_STATIC=1 "
                 f"and open the Vite dev server separately) -- API-only for now.")
    return serving.app


_server_error = []  # [] = still starting/ok, [exc_text] = the thread died


def start_server(port):
    try:
        import uvicorn
        app = build_app()
        _log(f"[desktop] backend built OK, starting uvicorn on port {port}")
        uvicorn.run(app, host=HOST, port=port, log_level="warning", log_config=None)
    except Exception:
        _server_error.append(traceback.format_exc())
        _log(f"[desktop] backend thread crashed:\n{_server_error[-1]}")


def start_admet_worker(port):
    """The learned ADMET-AI layer (backend/admet_service.py) is normally
       a SEPARATE process a web deployment's operator starts by hand --
       nobody does that for a double-clicked desktop app, so this runs
       it in its own background thread instead, giving the same
       automatic parity with zero setup. admet.py's own health-check +
       graceful-degrade design already covers admet-ai missing/broken,
       so failure here is caught broadly and never allowed to affect
       the main backend thread (unlike start_server(), whose failure
       IS fatal)."""
    try:
        if BACKEND_DIR not in sys.path:
            sys.path.insert(0, BACKEND_DIR)
        import uvicorn
        import admet_service
        _log(f"[desktop] starting ADMET-AI worker on port {port}")
        uvicorn.run(admet_service.app, host=HOST, port=port, log_level="warning", log_config=None)
    except Exception:
        _log(f"[desktop] ADMET-AI worker failed to start (non-fatal -- "
             f"the deterministic ADMET layer still works):\n{traceback.format_exc()}")


def main():
    if getattr(sys, "frozen", False):
        os.chdir(os.path.dirname(sys.executable))

        # See qsar-desktop's BUILD_WINDOWS.md for the full "why" -- a
        # PyInstaller 6.x onedir build's --add-data files land in
        # _internal/ (== FROZEN_ROOT/sys._MEIPASS), a SIBLING of the exe,
        # not the same folder as cwd (deliberately set to the exe's own
        # folder above, so models/docking_targets/ land somewhere a user
        # would actually find them, not buried in _internal/). Every one
        # of these env vars must be set before build_app() imports
        # anything that reads them.
        os.environ.setdefault("DOCKING_REGISTRY", os.path.join(FROZEN_ROOT, "docking_registry.json"))
        os.environ.setdefault("PANEL_RESULTS_CSV", os.path.join(FROZEN_ROOT, "panel_results_v2.csv"))
        os.environ.setdefault("CURATED_DATA_DIR", os.path.join(FROZEN_ROOT, "models", "curated"))
        # target_prediction_v2's serving data (~241MB total: the v2
        # index, the orthologue index, the target-name lookup) -- small
        # enough to bundle directly (see BUILD_WINDOWS.md), same pattern
        # as the three lines above. backend/target_prediction_v2.py
        # already reads these exact env vars (added when it was built,
        # specifically so a deployment could relocate the data without
        # editing that file) -- nothing new needed there.
        os.environ.setdefault(
            "TARGET_PREDICTION_V2_INDEX_DIR",
            os.path.join(FROZEN_ROOT, "target_prediction_v2_data", "v2_index"),
        )
        os.environ.setdefault(
            "TARGET_PREDICTION_V2_ORTHO_DIR",
            os.path.join(FROZEN_ROOT, "target_prediction_v2_data", "orthologue_index"),
        )
        os.environ.setdefault(
            "TARGET_PREDICTION_V2_NAME_MAP",
            os.path.join(FROZEN_ROOT, "target_prediction_v2_data", "target_name_map.json"),
        )

    _log(f"[desktop] starting (frozen={getattr(sys, 'frozen', False)}, cwd={os.getcwd()})")
    _put_bundled_binaries_on_path()
    _patch_gnina_rescorer()

    # Must happen BEFORE build_app() ever imports backend/app.py -> admet.py,
    # since admet.py reads ADMET_SERVICE_URL once, at import time.
    if "ADMET_SERVICE_URL" not in os.environ:
        admet_port = _free_port()
        os.environ["ADMET_SERVICE_URL"] = f"http://{HOST}:{admet_port}"
        threading.Thread(target=start_admet_worker, args=(admet_port,), daemon=True).start()

    port = _free_port()
    threading.Thread(target=start_server, args=(port,), daemon=True).start()
    if not _wait_up(port):
        if _server_error:
            _fatal(f"the backend crashed on startup:\n\n{_server_error[-1]}")
        else:
            _fatal(f"the backend did not respond on port {port} within 30s "
                    f"(no crash was caught -- it may be hanging on a slow import).")
    _log("[desktop] backend is up, serving")
    url = f"http://{HOST}:{port}/"
    try:
        import webview
        # Off by default in pywebview -- the embedded browser (WebView2 on
        # Windows) otherwise silently cancels every download, including
        # the frontend's Blob/<a download> CSV and PDB export buttons.
        webview.settings["ALLOW_DOWNLOADS"] = True
        window = webview.create_window("PhytoScreen", url, width=1280, height=860, min_size=(1000, 700))
        # App/window icon -- on Windows the .exe's own embedded icon
        # (PyInstaller's --icon flag, see BUILD_WINDOWS.md) is what
        # actually shows in the taskbar/title bar/Explorer; this covers
        # the platforms where pywebview reads a separate icon file
        # itself (primarily GTK on Linux). assets/icon.ico is bundled
        # the same way as docking_registry.json etc (see BUILD_WINDOWS.md/
        # the CI workflow's --add-data list) -- FROZEN_ROOT-relative when
        # frozen, repo-relative in dev.
        icon_path = os.path.join(FROZEN_ROOT, "assets", "icon.ico")
        webview.start(icon=icon_path if os.path.exists(icon_path) else None)          # blocks until the window closes
        return
    except ImportError:
        _log(f"[desktop] pywebview not installed. Open {url} in a browser, or `pip install pywebview`.")
    except Exception as e:
        _log(f"[desktop] could not open a native window ({e}). Open {url} in a browser instead.")
    while True:
        time.sleep(1)


if __name__ == "__main__":
    main()
