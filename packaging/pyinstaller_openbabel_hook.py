import os
import sys
from pathlib import Path

base = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
candidate = base / "openbabel" / "bin" / "data"
if candidate.exists():
    os.environ["BABEL_DATADIR"] = str(candidate)
