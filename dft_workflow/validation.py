from __future__ import annotations

import numpy as np
from ase import Atoms


def validate_export(
    molecule: Atoms | None,
    slab: Atoms | None,
    assembled: Atoms | None,
    flags: list[tuple[bool, bool, bool]],
    vacuum: float,
) -> list[str]:
    warnings: list[str] = []
    if molecule is None:
        warnings.append("No molecule loaded.")
    if slab is None:
        warnings.append("No slab built.")
    if vacuum < 10:
        warnings.append("Vacuum is below 10 A; check whether this is enough for your system.")
    if assembled is not None and _minimum_distance(assembled) < 0.7:
        warnings.append("At least two atoms are closer than 0.7 A.")
    if flags and all(all(axis for axis in flag) for flag in flags):
        warnings.append("No atoms are frozen; check the frozen-layer setting.")
    return warnings


def _minimum_distance(atoms: Atoms) -> float:
    positions = atoms.get_positions()
    if len(positions) < 2:
        return float("inf")
    minimum = float("inf")
    for index, position in enumerate(positions[:-1]):
        distances = np.linalg.norm(positions[index + 1 :] - position, axis=1)
        minimum = min(minimum, float(distances.min()))
    return minimum
