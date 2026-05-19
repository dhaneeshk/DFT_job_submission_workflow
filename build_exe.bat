@echo off
setlocal

if not exist ".venv\Scripts\python.exe" (
    echo Virtual environment not found. Run setup_env.bat first.
    exit /b 1
)

if exist ".venv\Lib\site-packages\openbabel\bin\data" (
    set "BABEL_DATADIR=%CD%\.venv\Lib\site-packages\openbabel\bin\data"
)

echo Building an experimental executable. The supported workflow is running from the virtual environment.
".venv\Scripts\python.exe" -m pip install -e ".[dev]"
".venv\Scripts\pyinstaller.exe" packaging\dft-workflow.spec --clean --noconfirm

echo Build complete. See dist\DFTWorkflow\DFTWorkflow.exe
