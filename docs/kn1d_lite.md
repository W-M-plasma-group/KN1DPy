# Using kn1d_lite

`kn1d_lite` is a simplified version of `kn1d` that is designed to be run on closed field lines (inside the separatrix).
It ignores molecules, and simply calculates the atomic neutral density for specified plasma profiles on the closed field lines.

There are two ways to use it depending on your setup.

---

## Option 1: Direct import

If KN1DPy runs in your local or installed python environment, you can import and call `kn1d_lite` directly:

```python
import numpy as np
from KN1DPy.kn1d_lite import kn1d_lite

# Toy profiles for illustrative purposes
x   = np.linspace(0, 0.05, 50)   # m, runs from separatrix inwards
n   = np.linspace(1e19, 5e19, 50) # m^-3
Te  = np.linspace(5.0, 100.0, 50) # eV
Ti  = Te.copy()
vxi = np.zeros_like(x)
mu = 2 # deutetrium
incident_n0 = 1e15 # m^-3. Density of inwardly-propagating neutrals.
energy_eV = 3.0 # assume all neutrals at separatrix have 3eV of energy.

result = kn1d_lite(
    x=x, mu=mu, Ti=Ti, Te=Te, n=n, vxi=vxi,
    incident_n0=incident_n0,
    energies_eV=[energy_eV],
    config_path='./config.json',
)

print(result.nH)   # neutral density profile on the atomic grid
```

---

## Option 2: Subprocess via file I/O

If you want to keep the KN1DPy pixi environment separate from your own working environment, you can call `run_kn1d_lite.py` as a subprocess. This script acts as a thin command-line wrapper around `kn1d_lite`: it reads inputs from a `.npz` file, runs the code, and writes the outputs to another `.npz` file. This allows `kn1d_lite` to interface with other codes that require different environments.

### Step 1: find your pixi Python path

After running `pixi install` in the KN1DPy directory, the Python executable will be at:

```
/path/to/KN1DPy/.pixi/envs/default/bin/python
```

You can find the exact path by running from the KN1DPy directory:

```bash
pixi run which python
```

### Step 2: write a wrapper function

In your own script, write a small wrapper that saves inputs to a temporary file, calls the `kn1d_lite` subprocess, and reads back the results:

```python
import subprocess
import tempfile
import numpy as np

KN1DPY_PYTHON = '/path/to/KN1DPy/.pixi/envs/default/bin/python'
KN1DPY_SCRIPT = '/path/to/KN1DPy/run_kn1d_lite.py'
KN1DPY_CONFIG = '/path/to/KN1DPy/config.json'

def call_kn1d_lite(x, mu, Ti, Te, n, vxi, incident_n0, velocities_ms, fractions,
                   config_path=KN1DPY_CONFIG):

    with tempfile.TemporaryDirectory() as tmp:
        input_path  = tmp + '/input.npz'
        output_path = tmp + '/output.npz'

        np.savez(input_path,
                 x=x, mu=mu, Ti=Ti, Te=Te, n=n, vxi=vxi,
                 incident_n0=incident_n0,
                 velocities_ms=np.atleast_1d(velocities_ms),
                 fractions=np.atleast_1d(fractions))

        subprocess.run([
            KN1DPY_PYTHON, KN1DPY_SCRIPT,
            '--input',  input_path,
            '--output', output_path,
            '--config', config_path,
        ], check=True)

        return np.load(output_path)
```

### Example usage

```python
import numpy as np

# Let's specify the velocity of the incoming neutrals on this run, rather than their energy. Stick with a single energy.

v_in = 15000 #m/s

result = call_kn1d_lite(x, mu=mu, Ti=Ti, Te=Te, n=n, vxi=vxi,
                        incident_n0=1e15,
                        velocities_ms=[v_in],
                        fractions=[1.0],
                        config_path='/path/to/my_config.json')
xH = result['xH']
fH = result['fH']
nH = result['nH']
```

### What `run_kn1d_lite.py` does

The script accepts three command-line arguments:

| Flag       | Description                          |
|------------|--------------------------------------|
| `--input`  | Path to a `.npz` file containing the input profiles |
| `--output` | Path where the output `.npz` will be written |
| `--config` | Path to the KN1DPy `config.json`     |

It loads the inputs, calls `kn1d_lite`, and saves all output arrays to the output file. The calling script is responsible for writing the input file and reading back the output — the paths can be anywhere on disk.

---

## Output arrays

When using direct import, outputs are attributes of the result object. When using the subprocess, they are keys in the output `.npz`. Note that `kn1d_lite` returns a subset of the full `kinetic_h` outputs — the higher-order moments and error diagnostics are omitted.

| Name          | Description |
|---------------|-------------|
| `xH`          | Spatial grid for atomic neutrals (m) |
| `fH`          | Full 2D velocity distribution function (nvr × nvx × nx) |
| `nH`          | Neutral density profile (m⁻³) |
| `GammaxH`     | Neutral particle flux (m⁻² s⁻¹) |
| `VxH`         | Neutral flow velocity (m s⁻¹) |
| `TH`          | Neutral temperature (eV) |
| `qxH_total`   | Total neutral energy flux (W m⁻²) |
| `Sion`        | Ionisation source (m⁻³ s⁻¹) |
| `fHBC`        | Boundary distribution function |
| `GammaxHBC`   | Incident neutral flux at the boundary |
| `vr`, `vx`    | Velocity grid arrays |
| `Tnorm`       | Normalisation temperature for the velocity grid (eV) |
