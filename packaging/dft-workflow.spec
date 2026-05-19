# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path
import openbabel

project_root = Path.cwd()
openbabel_data = Path(openbabel.__file__).resolve().parent / "bin" / "data"

a = Analysis(
    [str(project_root / "dft_workflow" / "app.py")],
    pathex=[str(project_root)],
    binaries=[],
    datas=[
        (str(project_root / "dft_workflow" / "Templates"), "dft_workflow/Templates"),
        (str(openbabel_data), "openbabel/bin/data"),
    ],
    hiddenimports=[
        "openbabel",
        "openbabel.openbabel",
        "ase.io.extxyz",
        "ase.io.xyz",
        "ase.io.vasp",
        "ase.io.sdf",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[str(project_root / "packaging" / "pyinstaller_openbabel_hook.py")],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="DFTWorkflow",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="DFTWorkflow",
)
