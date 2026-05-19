$ErrorActionPreference = "Stop"

if (-not (Test-Path -LiteralPath ".venv")) {
    python -m venv .venv
}

& ".\.venv\Scripts\python.exe" -m pip install --upgrade pip
& ".\.venv\Scripts\python.exe" -m pip install -e .

"Environment ready. Activate with: .\.venv\Scripts\Activate.ps1"
"Run the GUI with: dft-workflow"
