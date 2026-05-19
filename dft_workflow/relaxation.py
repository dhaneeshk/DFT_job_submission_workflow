from __future__ import annotations

import shutil
import subprocess
import tempfile
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from dft_workflow.assembly import AdsorbatePlacement, _place_molecule_above_slab, _rotate_about_center


TOP_LAYER_TOLERANCE = 0.35
DEFAULT_MAX_RELAXATION_ATOMS = 600
DEFAULT_RELAXATION_TIMEOUT_SECONDS = 3600
DEFAULT_GFNFF_STAGES = 30
DEFAULT_GFNFF_STEPS_PER_STAGE = 5


@dataclass(frozen=True)
class RelaxedAdsorbate:
    atoms: Atoms
    translation: tuple[float, float, float]
    rotation_deg: tuple[float, float, float]


@dataclass(frozen=True)
class RelaxationResult:
    success: bool
    message: str
    log: list[str]
    adsorbates: list[RelaxedAdsorbate]
    trajectory: list[Atoms]


def xtb_is_available(command: str = "xtb") -> bool:
    return shutil.which(command) is not None


def top_layer_indices(slab: Atoms, tolerance: float = TOP_LAYER_TOLERANCE) -> list[int]:
    z_positions = slab.get_positions()[:, 2]
    top_z = float(z_positions.max())
    return [index for index, z_value in enumerate(z_positions) if abs(float(z_value) - top_z) <= tolerance]


def build_gfnff_relaxation_system(
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms,
    tolerance: float = TOP_LAYER_TOLERANCE,
) -> tuple[Atoms, list[tuple[int, int]], list[int]]:
    placed_adsorbates = _placed_adsorbates(adsorbates, slab)
    relaxation_system = placed_adsorbates[0].copy()
    adsorbate_slices = []
    start = 0
    for placed_adsorbate in placed_adsorbates:
        if start == 0:
            relaxation_system = placed_adsorbate.copy()
        else:
            relaxation_system += placed_adsorbate
        stop = start + len(placed_adsorbate)
        adsorbate_slices.append((start, stop))
        start = stop

    top_indices = top_layer_indices(slab, tolerance)
    top_layer = slab[top_indices]
    relaxation_system += top_layer
    fixed_indices = list(range(start, start + len(top_layer)))
    relaxation_system.set_cell(slab.cell)
    relaxation_system.set_pbc((False, False, False))
    return relaxation_system, adsorbate_slices, fixed_indices


def run_gfnff_relaxation(
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms,
    xtb_command: str = "xtb",
    max_atoms: int = DEFAULT_MAX_RELAXATION_ATOMS,
    timeout_seconds: int = DEFAULT_RELAXATION_TIMEOUT_SECONDS,
    stages: int = DEFAULT_GFNFF_STAGES,
    steps_per_stage: int = DEFAULT_GFNFF_STEPS_PER_STAGE,
) -> RelaxationResult:
    log = ["Method: GFN-FF via xTB"]
    if not xtb_is_available(xtb_command):
        return RelaxationResult(False, "xTB executable not found", log + ["xTB executable not found in PATH."], [], [])

    relaxation_system, adsorbate_slices, fixed_indices = build_gfnff_relaxation_system(adsorbates, slab)
    log.append(f"Adsorbates: {sum(stop - start for start, stop in adsorbate_slices)} atoms")
    log.append(f"Fixed top-layer surface atoms: {len(fixed_indices)}")
    log.append(f"Total relaxation atoms: {len(relaxation_system)}")
    if len(relaxation_system) > max_atoms:
        return RelaxationResult(
            False,
            "Relaxation system too large",
            log + [f"Relaxation skipped: {len(relaxation_system)} atoms exceeds limit of {max_atoms}."],
            [],
            [],
        )

    with tempfile.TemporaryDirectory(prefix="dft_gfnff_") as tmpdir:
        workdir = Path(tmpdir)
        xcontrol = workdir / "constraints.inp"
        _write_xtb_control_file(xcontrol, fixed_indices, steps_per_stage)
        optimized = relaxation_system.copy()
        staged_trajectory = [optimized.copy()]
        log.append(f"Running staged xTB optimization: {stages} stages, {steps_per_stage} optimizer steps per stage.")

        for stage in range(1, stages + 1):
            input_xyz = workdir / f"stage_{stage:03d}.xyz"
            write(input_xyz, optimized)
            command = [xtb_command, input_xyz.name, "--gfnff", "--opt", "--verbose", "--input", xcontrol.name, "--no-restart"]
            log.append(f"Stage {stage}/{stages}: " + " ".join(command))
            try:
                completed = subprocess.run(
                    command,
                    cwd=workdir,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    capture_output=True,
                    timeout=timeout_seconds,
                    check=False,
                )
            except subprocess.TimeoutExpired:
                return RelaxationResult(False, "xTB relaxation timed out", log + [f"xTB timed out during stage {stage}."], [], [])

            if completed.returncode != 0:
                if completed.stdout:
                    log.extend(_tail_lines(completed.stdout, f"xTB stdout stage {stage}"))
                if completed.stderr:
                    log.extend(_tail_lines(completed.stderr, f"xTB stderr stage {stage}"))
                return RelaxationResult(False, f"xTB failed with return code {completed.returncode} in stage {stage}", log, [], [])

            output_xyz = workdir / "xtbopt.xyz"
            if not output_xyz.exists():
                return RelaxationResult(False, f"xTB did not write xtbopt.xyz in stage {stage}", log, [], [])
            next_optimized = read(output_xyz)
            displacement = float(np.linalg.norm(next_optimized.get_positions() - optimized.get_positions(), axis=1).max())
            optimized = next_optimized
            staged_trajectory.append(optimized.copy())
            log.append(f"Stage {stage}: max displacement {displacement:.6f} A")
            if displacement < 1e-4:
                log.append(f"Stopping staged optimization after stage {stage}; geometry change is below threshold.")
                break

        trajectory = _read_xtb_trajectory(workdir, log)
        if not trajectory and len(staged_trajectory) > 1:
            trajectory = staged_trajectory
            log.append(f"Using {len(trajectory)} staged xTB geometries for playback.")

    relaxed = _extract_relaxed_adsorbates(optimized, adsorbate_slices, slab)
    return RelaxationResult(True, "GFN-FF relaxation completed", log + ["Relaxation completed."], relaxed, trajectory)


def run_mmff_relaxation(
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms | None = None,
    max_atoms: int = DEFAULT_MAX_RELAXATION_ATOMS,
    steps: int = 500,
    batch_size: int = 10,
) -> RelaxationResult:
    return run_openbabel_relaxation("MMFF94", "MMFF94", adsorbates, slab, max_atoms, steps, batch_size)


def run_uff_relaxation(
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms | None = None,
    max_atoms: int = DEFAULT_MAX_RELAXATION_ATOMS,
    steps: int = 500,
    batch_size: int = 10,
) -> RelaxationResult:
    return run_openbabel_relaxation("UFF", "UFF", adsorbates, slab, max_atoms, steps, batch_size)


def run_openbabel_relaxation(
    forcefield_name: str,
    label: str,
    adsorbates: list[AdsorbatePlacement],
    slab: Atoms | None = None,
    max_atoms: int = DEFAULT_MAX_RELAXATION_ATOMS,
    steps: int = 500,
    batch_size: int = 10,
) -> RelaxationResult:
    log = [f"Method: Open Babel {label}"]
    if slab is None:
        log.append(f"{label} relaxes adsorbates only; the surface is not included.")
        total_atoms = sum(len(adsorbate.atoms) for adsorbate in adsorbates)
        relaxation_system = None
        adsorbate_slices = []
        fixed_indices = []
    else:
        log.append(f"{label} will attempt all adsorbates plus the fixed top surface layer.")
        relaxation_system, adsorbate_slices, fixed_indices = build_gfnff_relaxation_system(adsorbates, slab)
        total_atoms = len(relaxation_system)
        log.append(f"Fixed top-layer surface atoms: {len(fixed_indices)}")
    log.append(f"Relaxation atoms: {total_atoms}")
    if total_atoms > max_atoms:
        return RelaxationResult(
            False,
            f"{label} relaxation system too large",
            log + [f"Relaxation skipped: {total_atoms} atoms exceeds limit of {max_atoms}."],
            [],
            [],
        )

    try:
        import openbabel
        from openbabel import openbabel as ob
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        return RelaxationResult(False, "Open Babel is not available", log + [f"Open Babel import failed: {exc}"], [], [])
    _configure_openbabel_data_path(openbabel, log)

    if slab is not None:
        assert relaxation_system is not None
        try:
            optimized, trajectory = _run_openbabel_forcefield(
                relaxation_system, forcefield_name, fixed_indices, ob, log, steps, batch_size
            )
            relaxed_adsorbates = _extract_relaxed_adsorbates(optimized, adsorbate_slices, slab)
            return RelaxationResult(True, f"{label} relaxation completed", log + ["Relaxation completed."], relaxed_adsorbates, trajectory)
        except Exception as exc:  # pragma: no cover - depends on optional dependency/input chemistry
            return RelaxationResult(
                False,
                f"{label} failed for adsorbates plus surface layer",
                log + [str(exc), f"Try {'UFF' if label == 'MMFF94' else 'GFN-FF'} or disable the top-layer surface for Open Babel."],
                [],
                [],
            )

    relaxed_adsorbates = []
    for index, adsorbate in enumerate(adsorbates, start=1):
        try:
            relaxed_atoms, _trajectory = _run_openbabel_forcefield(adsorbate.atoms, forcefield_name, [], ob, log, steps, batch_size)
            log.append(f"Adsorbate {index}: Open Babel {label} completed.")
            relaxed_adsorbates.append(RelaxedAdsorbate(relaxed_atoms, adsorbate.translation, adsorbate.rotation_deg))
        except Exception as exc:  # pragma: no cover - depends on optional dependency/input chemistry
            return RelaxationResult(False, f"{label} failed for adsorbate {index}", log + [str(exc)], [], [])

    return RelaxationResult(True, f"{label} relaxation completed", log + ["Relaxation completed."], relaxed_adsorbates, [])


def _placed_adsorbates(adsorbates: list[AdsorbatePlacement], slab: Atoms) -> list[Atoms]:
    placed_adsorbates = []
    for adsorbate in adsorbates:
        moved_adsorbate = adsorbate.atoms.copy()
        _rotate_about_center(moved_adsorbate, adsorbate.rotation_deg)
        _place_molecule_above_slab(moved_adsorbate, slab, adsorbate.translation)
        placed_adsorbates.append(moved_adsorbate)
    return placed_adsorbates


def _extract_relaxed_adsorbates(
    optimized: Atoms,
    adsorbate_slices: list[tuple[int, int]],
    slab: Atoms,
) -> list[RelaxedAdsorbate]:
    slab_positions = slab.get_positions()
    slab_center_xy = slab_positions[:, :2].mean(axis=0)
    slab_top = float(slab_positions[:, 2].max())
    relaxed_adsorbates = []
    for start, stop in adsorbate_slices:
        relaxed_atoms = optimized[start:stop]
        positions = relaxed_atoms.get_positions()
        center_xy = positions[:, :2].mean(axis=0)
        bottom_z = float(positions[:, 2].min())
        translation = (
            float(center_xy[0] - slab_center_xy[0]),
            float(center_xy[1] - slab_center_xy[1]),
            float(bottom_z - slab_top),
        )
        relaxed_adsorbates.append(RelaxedAdsorbate(relaxed_atoms, translation, (0.0, 0.0, 0.0)))
    return relaxed_adsorbates


def _write_xtb_control_file(path: Path, fixed_indices: list[int], maxcycle: int) -> None:
    # xTB uses one-based atom indices in constraint files.
    fixed = ", ".join(str(index + 1) for index in fixed_indices)
    path.write_text(
        f"$fix\n  atoms: {fixed}\n$end\n"
        f"$opt\n  maxcycle={maxcycle}\n$end\n",
        encoding="utf-8",
    )


def _tail_lines(text: str, label: str, limit: int = 20) -> list[str]:
    lines = text.splitlines()
    tail = lines[-limit:]
    return [f"{label}:"] + tail


def _read_xtb_trajectory(workdir: Path, log: list[str]) -> list[Atoms]:
    produced_files = sorted(path.name for path in workdir.iterdir() if path.is_file())
    log.append("xTB produced files: " + (", ".join(produced_files) if produced_files else "none"))

    candidate_names = (
        "xtbopt.log",
        "xtb.trj",
        "xtbopt.trj",
        "trajectory.xyz",
        "xtbopt.xyz",
    )
    for name in candidate_names:
        path = workdir / name
        if not path.exists():
            continue
        try:
            frames = read(path, index=":")
        except Exception as exc:
            log.append(f"Could not read {name} as a trajectory: {exc}")
            continue
        if isinstance(frames, Atoms):
            frames = [frames]
        if len(frames) > 1:
            log.append(f"Loaded {len(frames)} relaxation frames from {name}.")
            return frames
        log.append(f"{name} contains only one frame; not using it for playback.")

    log.append("No multi-frame xTB trajectory found; final optimized geometry was applied.")
    return []


def _run_openbabel_forcefield(
    atoms: Atoms,
    forcefield_name: str,
    fixed_indices: list[int],
    ob,
    log: list[str],
    steps: int = 500,
    batch_size: int = 10,
) -> tuple[Atoms, list[Atoms]]:
    mol = _ase_atoms_to_openbabel_mol(atoms, ob)
    forcefield = ob.OBForceField.FindForceField(forcefield_name)
    if forcefield is None:
        raise ValueError(f"Open Babel force field not found: {forcefield_name}")

    constraints = ob.OBFFConstraints()
    for atom_index in fixed_indices:
        constraints.AddAtomConstraint(atom_index + 1)
    if fixed_indices:
        forcefield.SetConstraints(constraints)

    if not forcefield.Setup(mol):
        raise ValueError(f"Open Babel {forcefield_name} could not set up this system")

    trajectory = [atoms.copy()]
    forcefield.ConjugateGradientsInitialize(steps)
    previous_positions = atoms.get_positions()
    actual_steps = 0
    for completed_steps in range(batch_size, steps + batch_size, batch_size):
        forcefield.ConjugateGradientsTakeNSteps(batch_size)
        actual_steps = min(completed_steps, steps)
        forcefield.GetCoordinates(mol)
        frame = _openbabel_mol_to_ase_atoms(mol, atoms)
        positions = frame.get_positions()
        displacement = float(np.linalg.norm(positions - previous_positions, axis=1).max())
        trajectory.append(frame.copy())
        previous_positions = positions
        if displacement < 1e-4:
            log.append(
                f"Open Babel {forcefield_name}: stopped after {min(completed_steps, steps)} steps; "
                "geometry change is below threshold."
            )
            break

    forcefield.GetCoordinates(mol)
    log.append(f"Open Babel {forcefield_name}: completed {actual_steps} conjugate-gradient steps.")
    optimized = _openbabel_mol_to_ase_atoms(mol, atoms)
    if len(trajectory) > 1:
        log.append(f"Loaded {len(trajectory)} Open Babel relaxation frames for playback.")
    return optimized, trajectory


def _ase_atoms_to_openbabel_mol(atoms: Atoms, ob):
    xyz_lines = [str(len(atoms)), "generated by DFT workflow"]
    for symbol, position in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        xyz_lines.append(f"{symbol} {position[0]:.12f} {position[1]:.12f} {position[2]:.12f}")
    conversion = ob.OBConversion()
    conversion.SetInFormat("xyz")
    mol = ob.OBMol()
    if not conversion.ReadString(mol, "\n".join(xyz_lines) + "\n"):
        raise ValueError("Open Babel could not parse XYZ geometry")
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()
    return mol


def _openbabel_mol_to_ase_atoms(mol, template: Atoms) -> Atoms:
    positions = []
    for atom_index in range(1, mol.NumAtoms() + 1):
        atom = mol.GetAtom(atom_index)
        positions.append([atom.GetX(), atom.GetY(), atom.GetZ()])
    optimized = template.copy()
    optimized.set_positions(np.array(positions, dtype=float))
    return optimized


def _configure_openbabel_data_path(openbabel_module, log: list[str]) -> None:
    package_dir = Path(openbabel_module.__file__).resolve().parent
    candidates = [
        package_dir / "bin" / "data",
        package_dir.parent / "openbabel" / "bin" / "data",
    ]
    for candidate in candidates:
        if (candidate / "UFF.prm").exists() and (candidate / "mmff94.ff").exists():
            os.environ["BABEL_DATADIR"] = str(candidate)
            log.append(f"Open Babel data directory: {candidate}")
            return
    log.append("Open Babel data directory not found; force-field setup may fail.")
