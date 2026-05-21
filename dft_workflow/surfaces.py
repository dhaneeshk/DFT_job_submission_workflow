from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import fcc100, fcc110, fcc111


SUPPORTED_FACETS = ("111", "100", "110")


def build_fcc_slab(
    metal: str,
    facet: str,
    size_x: int,
    size_y: int,
    layers: int,
    vacuum: float,
    lattice_constant: float | None = None,
    bottom_vacuum: float = 0.0,
    orthogonal_111: bool = False,
) -> Atoms:
    if facet not in SUPPORTED_FACETS:
        raise ValueError(f"Unsupported facet '{facet}'. Choose from {', '.join(SUPPORTED_FACETS)}")
    if size_x < 1 or size_y < 1 or layers < 1:
        raise ValueError("Slab dimensions and layer count must be positive")
    if vacuum <= 0:
        raise ValueError("Vacuum must be positive")
    if bottom_vacuum < 0:
        raise ValueError("Bottom vacuum cannot be negative")
    if orthogonal_111 and facet != "111":
        raise ValueError("Rectangular cell shape is only available for the 111 facet")
    if orthogonal_111 and size_y % 2 != 0:
        raise ValueError("Rectangular fcc(111) slabs require an even Y size. Try Y = 2, 4, 6, ...")

    kwargs = {
        "symbol": metal,
        "size": (size_x, size_y, layers),
        "vacuum": None,
        "periodic": True,
    }
    if lattice_constant is not None and lattice_constant > 0:
        kwargs["a"] = lattice_constant

    if facet == "111":
        kwargs["orthogonal"] = orthogonal_111
        slab = fcc111(**kwargs)
    elif facet == "100":
        slab = fcc100(**kwargs)
    else:
        slab = fcc110(**kwargs)

    _position_slab_with_vacuum(slab, vacuum, bottom_vacuum)
    slab.set_pbc((True, True, False))
    return slab


def _position_slab_with_vacuum(slab: Atoms, vacuum: float, bottom_vacuum: float) -> None:
    positions = slab.get_positions()
    positions[:, 2] -= positions[:, 2].min()
    positions[:, 2] += bottom_vacuum
    slab.set_positions(positions)

    slab_top = positions[:, 2].max()
    cell = slab.cell.array.copy()
    cell[2] = np.array([0.0, 0.0, slab_top + vacuum])
    slab.set_cell(cell, scale_atoms=False)
