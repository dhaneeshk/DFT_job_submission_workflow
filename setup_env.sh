#!/usr/bin/env bash
set -euo pipefail

if [ ! -d ".venv" ]; then
  python3 -m venv .venv
fi

.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e .

echo "Environment ready. Activate with: source .venv/bin/activate"
echo "Run the GUI with: dft-workflow"
