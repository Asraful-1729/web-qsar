; PhytoScreen Windows installer (Inno Setup 6).
; Built by .github/workflows/build-windows.yml on a windows-latest runner
; from the PyInstaller onedir output (dist\PhytoScreen\); not meant to be
; run by hand except to sanity-check the packaging (see BUILD_WINDOWS.md).
;
; Ships without models/ or docking_targets/ (~101GB / 2.3GB) — those are
; pulled on demand by the Downloads tab (backend/downloads.py) from
; whatever PhytoScreen.bat's DOWNLOAD_BASE_URL points at. GNINA CNN
; rescoring (no official Windows build — see BUILD_WINDOWS.md) is
; likewise not bundled; PhytoScreen.bat's REMOTE_GNINA_URL points the
; installed app at a remote GPU server for that one optional feature.
;
; Adapted from the sibling qsar-desktop project's own proven installer
; script — same structure, plus target_prediction_v2_data\ (this
; project's own new dependency) in the file list, and a distinct AppId
; (a real separate product, not a variant of qsar-desktop's own install).

#define MyAppName "PhytoScreen"
#ifndef MyAppVersion
  #define MyAppVersion "0.0.0-dev"
#endif
#define MyAppExeName "PhytoScreen.bat"

[Setup]
AppId={{326B1FD0-0189-40E6-87C9-3F90714315F3}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
; Per-user, always-writable location — NOT {autopf} (Program Files).
; models/ and docking_targets/ are written next to the exe at runtime
; (see desktop.py's os.chdir + downloads.py), and Program Files is
; admin/UAC-protected: a normal (non-elevated) launch of the installed
; app can't write there at all, so every on-demand download fails with
; a WinError (PermissionError/WinError 5) for every user, every time.
; {localappdata} needs no elevation for either the installer or the app.
DefaultDirName={localappdata}\{#MyAppName}
PrivilegesRequired=lowest
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
OutputBaseFilename=PhytoScreenSetup
OutputDir=dist_installer
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
; Installer wizard/uninstaller icon — separate from the [Icons] entries
; below (those are the app's OWN shortcuts, shown after install).
SetupIconFile=..\assets\icon.ico

[Files]
Source: "..\dist\PhytoScreen\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion
Source: "PhytoScreen.bat"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
; IconFilename explicitly points at PhytoScreen.exe (which HAS the real
; icon embedded via PyInstaller's --icon flag, see BUILD_WINDOWS.md) --
; without this, Inno Setup defaults to the shortcut's own Filename
; target's icon, which is PhytoScreen.bat (a plain batch file launcher,
; see PhytoScreen.bat's own docstring for why the app launches through it
; rather than the exe directly) -- Windows shows a generic batch-file
; icon for that, regardless of what icon the real exe has. This was the
; actual cause of "no logo when desktop application show in desktop".
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\PhytoScreen.exe"; WorkingDir: "{app}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\PhytoScreen.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppName}"; Flags: nowait postinstall skipifsilent
