@echo off
setlocal
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
    echo Virtual environment not found. Run setup_env.bat first.
    exit /b 1
)

if exist ".venv\Lib\site-packages\openbabel\bin\data" (
    set "BABEL_DATADIR=%CD%\.venv\Lib\site-packages\openbabel\bin\data"
)

".venv\Scripts\dft-workflow.exe"
