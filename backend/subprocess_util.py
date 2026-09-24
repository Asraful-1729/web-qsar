"""
Windows-only fix for the desktop build: a subprocess.run()/Popen() call
to a CONSOLE program (vina.exe, obabel, prepare_receptor) pops a visible
black console window for that child process, EVEN THOUGH the parent app
itself has none (it's a PyInstaller --windowed build) -- Windows'
CreateProcess creates a fresh console for a console-subsystem child
unless explicitly told not to, regardless of the parent's own console
state. During docking/screening, this fires once per compound (each a
separate vina invocation), which is exactly why it was seen flashing
open/closed repeatedly until the run finished.

No effect on Linux/macOS (STARTF_USESHOWWINDOW/CREATE_NO_WINDOW don't
exist there) -- os.name == "nt" gates it, so the identical code path
runs unchanged for the plain web/server deployment.

Usage: subprocess.run([...], **hidden_subprocess_kwargs())
"""
import os
import subprocess


def hidden_subprocess_kwargs():
    if os.name != "nt":
        return {}
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = subprocess.SW_HIDE
    return {"startupinfo": startupinfo, "creationflags": subprocess.CREATE_NO_WINDOW}
