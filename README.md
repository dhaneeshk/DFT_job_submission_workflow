# DFT Job Submission Workflow

Python GUI for assembling adsorbate-on-metal-surface VASP job folders using ASE.

## Features

- Load one or more 3D adsorbates from `.xyz` or `.sdf`.
- Build fcc metal slabs for `111`, `100`, and `110` facets.
- Move and rotate each adsorbate independently.
- Freeze selected bottom slab layers for VASP selective dynamics.
- Export `POSCAR`, `INCAR`, `KPOINTS`, `sub`, and `POTCAR_COMMAND.txt`.
- Choose `Direct` or `Cartesian` POSCAR coordinates.
- Relax adsorbates using `GFN-FF via xTB`, `MMFF94 (Open Babel)`, or `UFF (Open Babel)`.
- Replay relaxation trajectories when frames are available.

## Recommended Install From GitHub

Use a Python virtual environment. This is the supported install path and is more reliable than a bundled executable for scientific Python packages such as ASE, Open Babel, and xTB.

If you have Git, clone the repository:

```bash
git clone https://github.com/dhaneeshk/DFT_job_submission_workflow.git
cd DFT_job_submission_workflow
```

If you do not have Git, download the files as a ZIP:

1. Open `https://github.com/dhaneeshk/DFT_job_submission_workflow` in a browser.
2. Click the green `Code` button.
3. Click `Download ZIP`.
4. Extract the ZIP somewhere convenient.
5. Open a terminal in the extracted folder.

Windows without PowerShell activation:

```bat
setup_env.bat
run_gui.bat
```

Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_env.ps1
.\.venv\Scripts\dft-workflow.exe
```

Linux/macOS:

```bash
./setup_env.sh
source .venv/bin/activate
dft-workflow
```

Manual setup:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
dft-workflow
```

On Windows, replace `source .venv/bin/activate` with `.\.venv\Scripts\Activate.ps1`, or run `.\.venv\Scripts\dft-workflow.exe` directly.

## External Tools

### xTB

`GFN-FF via xTB` requires the external `xtb` executable. xTB is not installed by this app.

Download/install options:

- xTB releases: `https://github.com/grimme-lab/xtb/releases`
- xTB setup documentation: `https://xtb-docs.readthedocs.io/en/latest/setup.html`
- Conda install:

```bash
conda install -c conda-forge xtb
```

If xTB is available on `PATH`, leave `dft_workflow_config.ini` unchanged.

Check with:

```bash
xtb --version
```

If `xtb --version` does not work, edit `dft_workflow_config.ini` in the project folder and set the full executable path:

```ini
[executables]
xtb = C:\Users\your-name\software\xtb\bin\xtb.exe
```

Linux/macOS example:

```ini
[executables]
xtb = /home/your-name/software/xtb/bin/xtb
```

If xTB is unavailable, use `MMFF94 (Open Babel)` or `UFF (Open Babel)`.

### VASP POTCAR Files

This app does not bundle or generate licensed VASP pseudopotentials. It writes `POTCAR_COMMAND.txt` using the POSCAR element order so the POTCAR can be generated on the appropriate server.

## Workflow

1. Add one or more 3D adsorbates from `.xyz` or `.sdf`.
2. Build an fcc metal slab for `111`, `100`, or `110`.
3. Select an adsorbate from the dropdown and move/rotate it using numeric controls.
4. Optionally relax adsorbates.
5. Freeze the bottom slab layers.
6. Export the VASP job folder.

The slab builder uses separate top and bottom vacuum values. `Top vacuum` is the empty space above the slab; `Bottom vacuum` lifts the bottom slab layer above `z=0`.

The relaxation panel can relax adsorbates with `GFN-FF via xTB`, `MMFF94 (Open Babel)`, or `UFF (Open Babel)`. Relaxation runs in a background thread so the GUI remains responsive. GFN-FF uses all adsorbates plus the top slab layer with the top layer fixed. MMFF94 and UFF use Open Babel; by default they relax adsorbates only, but `Include top surface layer` can be enabled to attempt the fixed top-layer subsystem.

## Templates

The exporter copies these files from `Templates/`:

- `INCAR`
- `KPOINTS`
- `sub`

`POSCAR` is generated from the assembled ASE structure. `POTCAR` is not generated locally.

## POSCAR Rules

- `Selective dynamics` is always written.
- Coordinates can be written as `Direct` or `Cartesian`.
- Adsorbate elements are written first.
- The slab metal is always written last.
- Frozen slab atoms receive `F F F`; mobile atoms receive `T T T`.

## Development

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

Build source/wheel distributions:

```bash
python -m build
```

## Before Publishing To GitHub

- Review `Templates/sub` for private usernames, emails, account names, and cluster paths.
- Do not commit generated VASP job folders.
- Do not commit VASP `POTCAR` files.
- Confirm `.venv/` is not committed.

## License

MIT License. See `LICENSE`.
