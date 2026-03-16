# KN1D
KN1D is a 1D-space, 2-D velocity neutral kinetic code developed by B. LaBombard (MIT).
This repo contains the updated python version, KN1DPy, of the original KN1D code.
Contact: njbrown@wm.edu

## NOTE
This translation is still in development and not fully tested. Certain functionality and accuracy may still be missing.

## Requirements
This translation was developed in python 3.12.3. It uses the pixi package manager. In order to set up correctly, clone the repository and run the following commands in the terminal:

```
curl -fsSL https://pixi.sh/install.sh | bash
cd KN1DPy
pixi install
pixi shell
```

## Examples

Example scripts are in `examples/C-Mod/`. All commands should be run from the **KN1DPy root directory**.

### Python (KN1DPy)

From the KN1DPy root directory, type the following commands one at a time:

```
pixi shell
python examples/C-Mod/run_cmod_python.py
```

### IDL

Start IDL from the **KN1DPy root directory** with the following command to ensure the correct shared library is loaded:

```bash
LD_LIBRARY_PATH=/path/to/KN1DPy/IDL:$LD_LIBRARY_PATH idl
```

Then type the following commands one at a time in the IDL prompt:

```
.compile IDL/bs2dr.pro
.compile IDL/kn1d.pro
.run examples/C-Mod/run_cmod_idl.pro
run_cmod_idl
```

Outputs are written to `examples/C-Mod/cmod_example_idl/`.

### Comparison plots

After running both versions, generate side-by-side comparison plots with:

```bash
pixi run python examples/C-Mod/plot_cmod_comparison.py
```

This saves `examples/C-Mod/cmod_comparison.png`.

## Limitations
Currently, anything using the Johnson-Hinov Tables are not working.
This includes Lyman_Alpha and Balmer Alpha, which will return 0 for the moment.
As such, the default choice for ionization coefficients has been set to Collrad Ionization.

There are also various other features that are currently not implemented.
These may be added later once the core program is completed.

## Configuration File
The file config.json is used to handle several settings

### Kinetic_H

- mesh_size - sets the size of the mesh generated for the kinetic_h calculations
- ion_rate - sets the method with which kinetic_h will perform ionization rate calculation
    - 'collrad' to use collrad ionization
    - 'jh' to use johnson-hinov ionization
    - 'janev' to use janev coefficients
    - KN1DPy will throw an exception if this value is not set to one of these three


### Kinetic_H2

- mesh_size - sets the size of the mesh generated for the kinetic_h2 calculations


### Collisions

- H2_H2_EL	- if set, then include H2 -> H2 elastic self collisions
- H2_P_EL	- if set, then include H2 -> H(+) elastic collisions
- H2_H_EL	- if set, then include H2 <-> H elastic collisions
- H2_HP_CX	- if set, then include H2 -> H2(+) charge exchange collisions
- H_H_EL	- if set, then include H -> H elastic self collisions
- H_P_CX	- if set, then include H -> H(+) charge exchange collisions
- H_P_EL	- if set, then include H -> H(+) elastic collisions
- SIMPLE_CX	- if set, then use CX source option (B): Neutrals are born
              in velocity with a distribution proportional to the local
              ion distribution function.
