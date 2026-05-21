from __future__ import annotations

from pathlib import Path

from ase import Atoms
from ase.io import read


MOLECULE_FORMATS = {".xyz": "xyz", ".sdf": "sdf"}
SUPPORTED_MOLECULE_FORMATS = set(MOLECULE_FORMATS)


def load_molecule(path: str | Path) -> Atoms:
    molecule_path = Path(path)
    suffix = molecule_path.suffix.lower()
    if suffix not in SUPPORTED_MOLECULE_FORMATS:
        supported = ", ".join(sorted(SUPPORTED_MOLECULE_FORMATS))
        raise ValueError(f"Unsupported molecule format '{suffix}'. Use one of: {supported}")

    molecule = read(molecule_path, format=MOLECULE_FORMATS[suffix])
    if not isinstance(molecule, Atoms):
        raise ValueError(f"Could not read a single molecule from {molecule_path}")
    if len(molecule) == 0:
        raise ValueError(f"Molecule file contains no atoms: {molecule_path}")
    return molecule
