from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class Assembly:
    atoms: Atoms
    mobile_flags: list[tuple[bool, bool, bool]]
    element_order: list[str]
    adsorbate_surface_distances: list[float]


@dataclass(frozen=True)
class AdsorbatePlacement:
    atoms: Atoms
    translation: tuple[float, float, float]
    rotation_deg: tuple[float, float, float]


def assemble_adsorbate_on_slab(
    molecule: Atoms,
    slab: Atoms,
    metal: str,
    translation: tuple[float, float, float],
    rotation_deg: tuple[float, float, float],
    frozen_layers: int,
    layer_tolerance: float = 0.35,
) -> Assembly:
    return assemble_adsorbates_on_slab(
        [AdsorbatePlacement(molecule, translation, rotation_deg)],
        slab,
        metal,
        frozen_layers,
        layer_tolerance,
    )


def assemble_adsorbates_on_slab(
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms,
    metal: str,
    frozen_layers: int,
    layer_tolerance: float = 0.35,
) -> Assembly:
    if not adsorbates:
        raise ValueError("At least one adsorbate is required")

    moved_slab = slab.copy()
    moved_adsorbates = []
    adsorbate_surface_distances = []
    for adsorbate in adsorbates:
        moved_adsorbate = adsorbate.atoms.copy()
        _rotate_about_center(moved_adsorbate, adsorbate.rotation_deg)
        _place_molecule_above_slab(moved_adsorbate, moved_slab, adsorbate.translation)
        adsorbate_surface_distances.append(_minimum_atom_distance(moved_adsorbate, moved_slab))
        moved_adsorbates.append(moved_adsorbate)

    combined = moved_adsorbates[0].copy()
    for moved_adsorbate in moved_adsorbates[1:]:
        combined += moved_adsorbate
    combined += moved_slab
    combined.set_cell(moved_slab.cell)
    combined.set_pbc(moved_slab.pbc)

    adsorbate_count = sum(len(adsorbate) for adsorbate in moved_adsorbates)
    frozen_surface_indices = _frozen_slab_indices(moved_slab, frozen_layers, layer_tolerance)

    flags: list[tuple[bool, bool, bool]] = []
    for index in range(len(combined)):
        if index < adsorbate_count:
            flags.append((True, True, True))
            continue
        slab_index = index - adsorbate_count
        flags.append((False, False, False) if slab_index in frozen_surface_indices else (True, True, True))

    ordered_atoms, ordered_flags, element_order = order_adsorbate_first_metal_last(
        combined, flags, adsorbate_count, metal
    )
    return Assembly(ordered_atoms, ordered_flags, element_order, adsorbate_surface_distances)


def order_adsorbate_first_metal_last(
    atoms: Atoms,
    flags: list[tuple[bool, bool, bool]],
    molecule_count: int,
    metal: str,
) -> tuple[Atoms, list[tuple[bool, bool, bool]], list[str]]:
    molecule_symbols = atoms.get_chemical_symbols()[:molecule_count]
    element_order = []
    for symbol in molecule_symbols:
        if symbol != metal and symbol not in element_order:
            element_order.append(symbol)
    element_order.append(metal)

    ordered_indices = []
    symbols = atoms.get_chemical_symbols()
    for symbol in element_order:
        ordered_indices.extend(index for index, value in enumerate(symbols) if value == symbol)

    return atoms[ordered_indices], [flags[index] for index in ordered_indices], element_order


def _rotate_about_center(atoms: Atoms, rotation_deg: tuple[float, float, float]) -> None:
    center = atoms.get_positions().mean(axis=0)
    for angle, axis in zip(rotation_deg, ("x", "y", "z")):
        if angle:
            atoms.rotate(angle, axis, center=center, rotate_cell=False)


def _place_molecule_above_slab(
    molecule: Atoms, slab: Atoms, translation: tuple[float, float, float]) -> None:
    tx, ty, height = translation
    slab_positions = slab.get_positions()
    molecule_positions = molecule.get_positions()
    slab_center_xy = slab_positions[:, :2].mean(axis=0)
    molecule_center_xy = molecule_positions[:, :2].mean(axis=0)
    slab_top = slab_positions[:, 2].max()
    molecule_bottom = molecule_positions[:, 2].min()

    delta = np.array([
        slab_center_xy[0] - molecule_center_xy[0] + tx,
        slab_center_xy[1] - molecule_center_xy[1] + ty,
        slab_top - molecule_bottom + height,
    ])
    molecule.translate(delta)


def _minimum_atom_distance(adsorbate: Atoms, slab: Atoms) -> float:
    adsorbate_positions = adsorbate.get_positions()
    slab_positions = slab.get_positions()
    minimum = float("inf")
    for position in adsorbate_positions:
        distances = np.linalg.norm(slab_positions - position, axis=1)
        minimum = min(minimum, float(distances.min()))
    return minimum


def _frozen_slab_indices(slab: Atoms, frozen_layers: int, tolerance: float) -> set[int]:
    if frozen_layers <= 0:
        return set()

    z_values = slab.get_positions()[:, 2]
    layer_values: list[float] = []
    for z_value in sorted(z_values):
        if not layer_values or abs(z_value - layer_values[-1]) > tolerance:
            layer_values.append(float(z_value))

    frozen_cutoff_layers = layer_values[:frozen_layers]
    frozen_indices = set()
    for index, z_value in enumerate(z_values):
        if any(abs(z_value - layer_z) <= tolerance for layer_z in frozen_cutoff_layers):
            frozen_indices.add(index)
    return frozen_indices
