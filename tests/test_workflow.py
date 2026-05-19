from pathlib import Path

import pytest

from dft_workflow.assembly import AdsorbatePlacement, assemble_adsorbate_on_slab, assemble_adsorbates_on_slab
from dft_workflow.molecule_io import load_molecule
from dft_workflow.relaxation import (
    _ase_atoms_to_openbabel_mol,
    build_gfnff_relaxation_system,
    run_gfnff_relaxation,
    run_uff_relaxation,
    top_layer_indices,
)
from dft_workflow.surfaces import build_fcc_slab
from dft_workflow.vasp_export import write_job_folder


ROOT = Path(__file__).resolve().parents[1]


def test_xyz_loads_example_molecule():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")

    assert len(molecule) == 92
    assert molecule.get_chemical_symbols()[0] == "H"


def test_assembly_keeps_slab_cell_and_metal_last():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    slab = build_fcc_slab("Cu", "111", 2, 2, 2, 12.0)

    assembly = assemble_adsorbate_on_slab(
        molecule,
        slab,
        "Cu",
        translation=(0.0, 0.0, 2.5),
        rotation_deg=(0.0, 0.0, 0.0),
        frozen_layers=1,
    )

    assert len(assembly.atoms) == len(molecule) + len(slab)
    assert assembly.atoms.cell.volume > 0
    assert assembly.element_order[-1] == "Cu"
    assert any(flag == (False, False, False) for flag in assembly.mobile_flags)
    assert assembly.adsorbate_surface_distances[0] >= 2.5


def test_slab_bottom_is_at_cell_bottom_with_vacuum_above():
    vacuum = 14.0
    slab = build_fcc_slab("Cu", "111", 3, 3, 3, vacuum)
    z_positions = slab.get_positions()[:, 2]

    assert z_positions.min() == pytest.approx(0.0)
    assert slab.cell.lengths()[2] - z_positions.max() == pytest.approx(vacuum)


def test_slab_can_include_bottom_vacuum():
    top_vacuum = 14.0
    bottom_vacuum = 3.0
    slab = build_fcc_slab("Cu", "111", 3, 3, 3, top_vacuum, bottom_vacuum=bottom_vacuum)
    z_positions = slab.get_positions()[:, 2]

    assert z_positions.min() == pytest.approx(bottom_vacuum)
    assert slab.cell.lengths()[2] - z_positions.max() == pytest.approx(top_vacuum)


def test_export_writes_expected_job_files(tmp_path):
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    slab = build_fcc_slab("Ag", "100", 2, 2, 2, 12.0)
    assembly = assemble_adsorbate_on_slab(
        molecule,
        slab,
        "Ag",
        translation=(0.0, 0.0, 2.5),
        rotation_deg=(0.0, 0.0, 0.0),
        frozen_layers=1,
    )

    write_job_folder(
        tmp_path,
        assembly.atoms,
        assembly.mobile_flags,
        assembly.element_order,
        ROOT / "Templates",
        "python potcar_generator.py {elements} > POTCAR",
    )

    assert (tmp_path / "POSCAR").exists()
    assert (tmp_path / "INCAR").exists()
    assert (tmp_path / "KPOINTS").exists()
    assert (tmp_path / "sub").exists()
    assert (tmp_path / "POTCAR_COMMAND.txt").read_text(encoding="utf-8").splitlines()[1] == (
        "Element order: H C N O Ag"
    )

    poscar = (tmp_path / "POSCAR").read_text(encoding="utf-8")
    assert "Selective dynamics" in poscar
    assert "Direct" in poscar
    assert "F   F   F" in poscar


def test_export_can_write_cartesian_poscar(tmp_path):
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    slab = build_fcc_slab("Cu", "111", 2, 2, 2, 12.0)
    assembly = assemble_adsorbate_on_slab(
        molecule,
        slab,
        "Cu",
        translation=(0.0, 0.0, 2.5),
        rotation_deg=(0.0, 0.0, 0.0),
        frozen_layers=1,
    )

    write_job_folder(
        tmp_path,
        assembly.atoms,
        assembly.mobile_flags,
        assembly.element_order,
        ROOT / "Templates",
        "python potcar_generator.py {elements} > POTCAR",
        coordinate_mode="Cartesian",
    )

    lines = (tmp_path / "POSCAR").read_text(encoding="utf-8").splitlines()
    assert lines[8] == "Cartesian"
    first_position = [float(value) for value in lines[9].split()[:3]]
    assert first_position == pytest.approx(assembly.atoms.get_positions()[0])


def test_multiple_adsorbates_export_as_one_system(tmp_path):
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    second_molecule = molecule[:3]
    slab = build_fcc_slab("Au", "110", 2, 2, 2, 12.0)

    assembly = assemble_adsorbates_on_slab(
        [
            AdsorbatePlacement(molecule, (0.0, 0.0, 2.5), (0.0, 0.0, 0.0)),
            AdsorbatePlacement(second_molecule, (3.0, 0.0, 2.5), (0.0, 0.0, 45.0)),
        ],
        slab,
        "Au",
        frozen_layers=1,
    )

    assert len(assembly.atoms) == len(molecule) + len(second_molecule) + len(slab)
    assert assembly.element_order[-1] == "Au"
    assert len(assembly.adsorbate_surface_distances) == 2
    assert all(distance >= 2.5 for distance in assembly.adsorbate_surface_distances)

    write_job_folder(
        tmp_path,
        assembly.atoms,
        assembly.mobile_flags,
        assembly.element_order,
        ROOT / "Templates",
        "python potcar_generator.py {elements} > POTCAR",
    )

    command = (tmp_path / "POTCAR_COMMAND.txt").read_text(encoding="utf-8")
    assert "Element order: H C N O Au" in command


def test_relaxation_system_uses_adsorbates_and_top_layer_only():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")[:5]
    slab = build_fcc_slab("Cu", "111", 3, 3, 3, 12.0)
    top_indices = top_layer_indices(slab)

    relaxation_system, adsorbate_slices, fixed_indices = build_gfnff_relaxation_system(
        [AdsorbatePlacement(molecule, (0.0, 0.0, 2.5), (0.0, 0.0, 0.0))],
        slab,
    )

    assert len(relaxation_system) == len(molecule) + len(top_indices)
    assert adsorbate_slices == [(0, len(molecule))]
    assert fixed_indices == list(range(len(molecule), len(molecule) + len(top_indices)))


def test_gfnff_reports_missing_xtb_before_running():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")[:5]
    slab = build_fcc_slab("Cu", "111", 2, 2, 2, 12.0)

    result = run_gfnff_relaxation(
        [AdsorbatePlacement(molecule, (0.0, 0.0, 2.5), (0.0, 0.0, 0.0))],
        slab,
        xtb_command="definitely_missing_xtb_binary",
    )

    assert not result.success
    assert "not found" in result.message


def test_gfnff_blocks_large_relaxation_system():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    slab = build_fcc_slab("Cu", "111", 3, 3, 3, 12.0)

    result = run_gfnff_relaxation(
        [AdsorbatePlacement(molecule, (0.0, 0.0, 2.5), (0.0, 0.0, 0.0))],
        slab,
        max_atoms=1,
    )

    if result.message == "xTB executable not found":
        pytest.skip("xTB is not available in this test environment")
    assert not result.success
    assert "too large" in result.message


def test_openbabel_xyz_conversion_perceives_connectivity():
    pytest.importorskip("openbabel")
    from openbabel import openbabel as ob

    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    mol = _ase_atoms_to_openbabel_mol(molecule, ob)

    assert mol.NumAtoms() == len(molecule)
    assert mol.NumBonds() > 0


def test_uff_blocks_large_relaxation_system():
    molecule = load_molecule(ROOT / "Model_compound1_TPP_oOBA_freebase.xyz")
    slab = build_fcc_slab("Cu", "111", 3, 3, 3, 12.0)

    result = run_uff_relaxation(
        [AdsorbatePlacement(molecule, (0.0, 0.0, 2.5), (0.0, 0.0, 0.0))],
        slab,
        max_atoms=1,
    )

    assert not result.success
    assert "too large" in result.message
