from __future__ import annotations

import configparser
from pathlib import Path


CONFIG_FILE_NAME = "dft_workflow_config.ini"


def load_xtb_command(config_path: str | Path | None = None) -> str:
    path = Path(config_path) if config_path is not None else Path.cwd() / CONFIG_FILE_NAME
    if not path.exists():
        return "xtb"

    parser = configparser.ConfigParser()
    parser.read(path, encoding="utf-8")
    command = parser.get("executables", "xtb", fallback="").strip()
    return command or "xtb"
