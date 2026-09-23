@echo off
rem Launcher installed next to PhytoScreen.exe (see phytoscreen.iss) — sets
rem DOWNLOAD_BASE_URL (models/docking_targets on-demand downloads, see
rem backend/downloads.py) and REMOTE_GNINA_URL/REMOTE_GNINA_TOKEN (optional
rem GNINA CNN rescoring via a remote GPU server, see
rem backend/remote_gnina_service.py and BUILD_WINDOWS.md — GNINA has no
rem official Windows build, so this is the only way it's ever offered on
rem Windows), then starts the app. The __...__ placeholders are substituted
rem at CI build time (see .github/workflows/build-windows.yml) — edit this
rem file after install to point at different values without a rebuild.
set DOWNLOAD_BASE_URL=__DOWNLOAD_BASE_URL__
set REMOTE_GNINA_URL=__REMOTE_GNINA_URL__
set REMOTE_GNINA_TOKEN=__REMOTE_GNINA_TOKEN__
start "" "%~dp0PhytoScreen.exe"
