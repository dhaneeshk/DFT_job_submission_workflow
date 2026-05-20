@echo off
setlocal
cd /d "%~dp0"

if not exist ".venv" (
    python -m venv .venv
)

".venv\Scripts\python.exe" -m pip install --upgrade pip
".venv\Scripts\python.exe" -m pip install -e .

echo Environment ready.
echo Run the GUI with: run_gui.bat
