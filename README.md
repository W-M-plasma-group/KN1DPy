# KN1D

KN1D is a 1D-space, 2D-velocity kinetic neutrals code developed by B. LaBombard (MIT PSFC). KN1DPy is a Python translation of the original IDL code (the original IDL version is also included in the `IDL/` directory of this repository).

This translation was written by Nicholas Brown at William & Mary (njbrown@wm.edu), and is now maintained by Jamie Dunsmore at MIT (jduns@mit.edu).

A comprehensive description of the algorithm can be found in this [PSFC report (LaBombard, 2001)](https://library.psfc.mit.edu/catalog/reports/2000/01rr/01rr003/01rr003_full.pdf).

---

## Getting Started

### Installation

KN1DPy requires Python ≥ 3.7. Three install options are available:

**Option 1 — pip (recommended for most users)**

Works with your existing Python environment and integrates easily with other packages.

```bash
git clone https://github.com/your-org/KN1DPy
cd KN1DPy
pip install -e .
```

**Option 2 — pixi (recommended if you use conda-based environments)**

[pixi](https://pixi.sh) manages a self-contained conda environment. Useful if you need conda-forge packages alongside KN1DPy.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
cd KN1DPy
pixi install
pixi shell
```

**Option 3 — pixi with locked environment (for exact reproduction of reference results)**

A `pixi.lock` file pinning the exact environment used to generate the reference outputs is provided in `lock_files/`. Follow the same steps as Option 2, but after `cd KN1DPy` and before `pixi install`, copy the lock file to the repo root:

```bash
cp lock_files/pixi.lock .
```

### Running your own case

The example run scripts (e.g. `examples/C-Mod/run_cmod_test_python.py`) can be used as templates. The inputs are loaded from a `.sav` file, but any format (`.pkl`, `.npz`, etc.) can be used — just read the data into numpy arrays and pass them directly to `kn1d()`. See [Inputs](#inputs) for a full description of the required arguments.

---

## Examples

Three example cases are provided in the `examples/` directory: **C-Mod**, **MAST-U**, and **DIII-D**. Each includes input files, run scripts for both Python and IDL, pre-computed output files, and scripts for comparing and plotting the results.

> **Note:** The Python version will only reproduce the IDL results exactly when run using the `config.json` that is saved alongside each set of outputs in `examples/<device>/python_output/`. 

All commands must be run from the **KN1DPy root directory**.

Each example directory has the following structure:

```
examples/<case>/
  input/                          input .sav file
  python_output/                  KN1DPy output files and the config used
  IDL_output/                     IDL output files
  run_<case>_test_python.py       Python run script
  run_<case>_test_idl.pro         IDL run script
  plot_<case>_test_comparison.py  side-by-side comparison plots
  check_<case>_outputs.py         numerical comparison table
```

---

### Python

The python cases can be run straightforwardly from the root directory using the following commands:

```bash
python examples/C-Mod/run_cmod_test_python.py
python examples/MAST-U/run_mastu_test_python.py
python examples/DIII-D/run_diiid_test_python.py
```

---

### IDL

> **Note:** An IDL licence is required. If you do not have one, you can still process the pre-computed IDL outputs that are already included in each `examples/<device>/IDL_output/` directory.

IDL requires a Fortran shared library (`fast_b2val.so`) located in the `IDL/` directory. The `LD_LIBRARY_PATH` environment variable must be set when launching IDL so that it can find this library. Launch IDL from the terminal:

```bash
LD_LIBRARY_PATH=$(pwd)/IDL:$LD_LIBRARY_PATH idl
```

Then at the IDL prompt, replacing `<device>` and `<proc>` with the values from the table below:

```idl
.compile IDL/kn1d.pro
.run examples/<device>/run_<device>_test_idl.pro
<proc>
```

| Case   | `<device>` | `<proc>`             |
|--------|------------|----------------------|
| C-Mod  | `C-Mod`    | `run_cmod_test_idl`  |
| MAST-U | `MAST-U`   | `run_mastu_test_idl` |
| DIII-D | `DIII-D`   | `run_diiid_test_idl` |

> **Warning:** Do not place any `.so` files in the KN1DPy root directory. IDL resolves shared library paths relative to the working directory first, so a stray `fast_b2val.so` in the root will be picked up before the one in `IDL/` and may cause the code to crash.

---

## Inputs

The main entry point is `kn1d()` in `KN1DPy/kn1d.py`. All array inputs should be 1D numpy arrays of length `nx`, defined on the radial coordinate grid `x`.

| Parameter  | Type            | Units         | Description |
|------------|-----------------|---------------|-------------|
| `x`        | ndarray (nx)    | m             | x-coordinate |
| `xlimiter` | float           | m             | Limiter position |
| `xsep`     | float           | m             | Separatrix position |
| `GaugeH2`  | float           | mTorr         | Molecular neutral pressure at the wall |
| `mu`       | float           | —             | Ion mass: 1 = hydrogen, 2 = deuterium |
| `Ti`       | ndarray (nx)    | eV            | Ion temperature profile |
| `Te`       | ndarray (nx)    | eV            | Electron temperature profile |
| `n`        | ndarray (nx)    | m⁻³           | Electron density profile |
| `vxi`      | ndarray (nx)    | m s⁻¹         | Plasma flow velocity (negative = towards wall). Generally set this to 0 |
| `LC`       | ndarray (nx)    | m             | Connection length to nearest limiter along field lines (0 = infinity) |
| `PipeDia`  | ndarray (nx)    | m             | Effective diameter of the pressure gauge pipe for side-wall collisions (0 = disabled, which is the most common setting) |

---

## Outputs

When `File` is specified, `kn1d()` writes four `.npz` files and a copy of the active `config.json` to the output directory. The `.npz` files can be loaded with `numpy.load()`.

| File | Description |
|------|-------------|
| `KN1D_H.npz` | Atomic hydrogen results. The core output is `fH`, the 2D velocity distribution function (vr × vx) on the atomic spatial grid `xH`. All atomic quantities — density `nH`, particle flux `GammaxH`, temperature `TH`, ionization source `Sion`, emissivities, and higher-order moments — are derived from this distribution function. |
| `KN1D_H2.npz` | Molecular hydrogen results. The core output is `fH2`, the 2D velocity distribution function (vr × vx) on the molecular spatial grid `xH2`. Derived quantities include `nH2`, `GammaxH2`, `TH2`, and the atomic and ion source terms `SH` and `SP`. |
| `KN1D_input.npz` | The input profiles (`Ti`, `Te`, `n`, etc.) interpolated onto both the atomic and molecular spatial grids, together with the velocity grid arrays (`vrA`, `vxA`, `vrM`, `vxM`). Useful for plotting inputs and outputs on the same axes. |
| `KN1D_mesh.npz` | The raw velocity and spatial grid parameters used internally. |
| `config.json` | A copy of the configuration used for this run, so the outputs are fully self-documenting. |

---

## Configuration

As well as the inputs, the settings for each run (e.g choice of atomic rate coefficients and mesh sizes) can be controlled via `config.json` in the root directory (or a custom path passed via `config_path`). The settings that can be changed in the `config.json` file are as follows...

### `kinetic_h` and `kinetic_h2`

| Key              | Description |
|------------------|-------------|
| `mesh_size`      | Number of velocity grid points. Likely needs to be increased from default for cases with > 500 eV pedestals. |
| `grid_fctr`      | Scales the physics-based maximum spatial grid spacing. Smaller values give a finer mesh. Default `0.3` (matching the IDL default). Lower values mean finer resolution. |
| `extra_energy_bins_eV` | Additional velocity grid points at the specified energies (eV). `kinetic_h` has no hardcoded energy bins, but `kinetic_h2` already has hardcoded energy bins at `0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0`, so any `kinetic_h2` inputs here will be in addition to these. Default `[]`. |
| `ion_rate`       | Ionization rate method: `"adas"` (recommended), `"jh"` (Johnson–Hinnov), `"collrad"`, or `"janev"`. (`kinetic_h` only) |

### `collisions`

Each flag enables or disables a specific collision channel:

| Key         | Description |
|-------------|-------------|
| `H2_H2_EL`  | H₂ → H₂ elastic self-collisions |
| `H2_P_EL`   | H₂ → H⁺ elastic collisions |
| `H2_H_EL`   | H₂ ↔ H elastic collisions |
| `H2_P_CX`   | H₂ → H₂⁺ charge exchange |
| `H_H_EL`    | H → H elastic self-collisions |
| `H_P_CX`    | H → H⁺ charge exchange |
| `H_P_EL`    | H → H⁺ elastic collisions |
| `SIMPLE_CX` | Use simplified CX collisions (neutrals born with ion distribution) |

---

## kn1d_lite

`kn1d_lite` is a simplified version of `kn1d` that is designed to be run on closed field lines (inside the separatrix).
It ignores molecules, and simply calculates the atomic neutral density for specified plasma profiles on the closed field lines. See [docs/kn1d_lite.md](docs/kn1d_lite.md) for a full description and instructions on how to run.
