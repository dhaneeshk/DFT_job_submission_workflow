# Contributing

## Development Setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
```

On Windows, use `setup_env.bat` or `setup_env.ps1` instead.

## Before Opening A Pull Request

- Run `python -m pytest`.
- Do not commit `.venv/`, generated job folders, or VASP `POTCAR` files.
- Keep VASP templates free of private credentials before publishing.
