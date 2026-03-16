"""
C-Mod KN1D example - Python (KN1DPy)
=====================================
Run this script from the KN1DPy root directory so that the package is found:

    python examples/C-Mod/run_cmod_python.py

Inputs are loaded from the shared .sav file (same source as the IDL script).

Config settings (mesh size, collision flags, etc.) are read from CONFIG_PATH
below. A copy of the active config is automatically saved alongside the outputs
as config.json inside the output directory.
"""

from KN1DPy.kn1d import kn1d
from scipy.io import readsav
import numpy as np
import time

np.set_printoptions(linewidth=225)

# ------------------------------------------------------------------ #
#  Paths (edit these as needed)
# ------------------------------------------------------------------ #
data_file   = 'examples/C-Mod/cmod_test_in.sav'
config_path = './config.json'   # root config; swap for a per-example json if needed
output_file = 'examples/C-Mod/cmod_example'  # outputs saved in examples/C-Mod/cmod_example/

# ------------------------------------------------------------------ #
#  Load inputs
# ------------------------------------------------------------------ #
print(f'Loading inputs: {data_file}')
sav = readsav(data_file)

# Unit conversions (sav file stores keV and 1e20 m^-3)
Ti = sav['Ti'] * 1e3     # keV  -> eV
Te = sav['Te'] * 1e3     # keV  -> eV
n  = sav['n']  * 1e20    # 1e20 m^-3 -> m^-3

# ------------------------------------------------------------------ #
#  Run KN1D
# ------------------------------------------------------------------ #
print('Starting KN1D ...')
t0 = time.time()

results = kn1d(
    x        = sav['x'],
    xlimiter = sav['xlimiter'],
    xsep     = sav['xsep'],
    GaugeH2  = sav['GaugeH2'],
    mu       = sav['mu'],
    Ti       = Ti,
    Te       = Te,
    n        = n,
    vxi      = sav['vxi'],
    LC       = sav['LC'],
    PipeDia  = sav['PipeDia'],
    max_gen        = 100,
    compute_errors = 1,
    debrief        = 1,
    Hdebrief       = 1,
    H2debrief      = 1,
    Hdebug         = 0,
    H2debug        = 0,
    File           = output_file,
    config_path    = config_path,
)

print(f'Elapsed time: {time.time() - t0:.1f} s')

# ------------------------------------------------------------------ #
#  Quick summary of key outputs
# ------------------------------------------------------------------ #
print('\n--- Results summary ---')
print(f'  xH shape  : {results.xH.shape}')
print(f'  nH peak   : {results.nH.max():.3e} m^-3')
print(f'  Sion peak : {results.Sion.max():.3e} m^-3 s^-1')
print(f'  Results written to {output_file}/ (KN1D_input.npz, KN1D_mesh.npz, KN1D_H2.npz, KN1D_H.npz, config.json)')
