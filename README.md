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

## Installation

This app is meant to be run from a Python virtual environment. A virtual environment is a local Python folder that keeps this app's packages separate from the rest of your computer.

### 1. Install Python

Install Python 3.10, 3.11, or 3.12 from `https://www.python.org/downloads/`.

On Windows, enable `Add python.exe to PATH` during installation if the installer offers that option.

### 2. Download This App

If you do not use Git:

1. Open `https://github.com/dhaneeshk/DFT_job_submission_workflow` in a browser.
2. Click the green `Code` button.
3. Click `Download ZIP`.
4. Extract the ZIP somewhere convenient, such as your Desktop or Documents folder.
5. Open the extracted folder.

On Windows 11, right-click inside the folder and choose `Open in Terminal`. On older Windows versions, open PowerShell and use `cd` to move into the extracted folder before running the setup command.

If you use Git:

```bash
git clone https://github.com/dhaneeshk/DFT_job_submission_workflow.git
cd DFT_job_submission_workflow
```

### 3. Set Up The Python Environment

Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_env.ps1
```

macOS/Linux terminal:

```bash
./setup_env.sh
```

This creates a `.venv` folder and installs the required Python packages.

### 4. Start The App

Windows PowerShell:

```powershell
.\.venv\Scripts\dft-workflow.exe
```

macOS/Linux terminal:

```bash
source .venv/bin/activate
dft-workflow
```

If you close the terminal and want to run the app again later, return to the project folder and run the command from this step again. You do not need to repeat the setup step unless you download a new version.

### Troubleshooting Installation

If `python` is not recognized, install Python again and make sure it is added to `PATH`, or use the Python Launcher on Windows:

```powershell
py -m venv .venv
```

If package installation fails, try updating pip manually:

```powershell
.\.venv\Scripts\python.exe -m pip install --upgrade pip
```

## External Tools

### xTB

The `GFN-FF via xTB` relaxation method requires the external `xtb` program. xTB is not installed by this app.

You can still use the app without xTB. If xTB is not installed, use `MMFF94 (Open Babel)` or `UFF (Open Babel)` in the relaxation menu.

Where to get xTB:

- xTB releases: `https://github.com/grimme-lab/xtb/releases`
- xTB setup documentation: `https://xtb-docs.readthedocs.io/en/latest/setup.html`

If you use conda, this is usually the easiest install command:

```bash
conda install -c conda-forge xtb
```

After installing xTB, check whether your terminal can find it:

```powershell
xtb --version
```

If that command works, leave `dft_workflow_config.ini` unchanged.

If that command does not work, you can tell this app exactly where xTB is installed:

1. Open `dft_workflow_config.ini` in the project folder.
2. Find the `[executables]` section.
3. Set `xtb` to the full path of the xTB executable.

Windows example:

```ini
[executables]
xtb = C:\Users\your-name\software\xtb\bin\xtb.exe
```

Linux/macOS example:

```ini
[executables]
xtb = /home/your-name/software/xtb/bin/xtb
```

Save the file, then restart the GUI.

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

For the `111` facet, the `111 cell shape` control selects either the primitive rhombic surface cell or a rectangular/orthogonal cell. Rectangular `111` slabs require an even `Size Y` value, such as `2`, `4`, or `6`.

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
