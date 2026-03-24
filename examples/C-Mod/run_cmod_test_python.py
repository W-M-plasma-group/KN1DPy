"""
C-Mod KN1D example - Python (KN1DPy)
=====================================
Usage:

    python examples/C-Mod/run_cmod_test_python.py

Inputs are loaded from the shared .sav file (same source as the IDL script).

Config settings (mesh size, collision flags, etc.) are read from CONFIG_PATH
below. A copy of the active config is automatically saved alongside the outputs
as config.json inside the output directory.
"""

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from KN1DPy.kn1d import kn1d
from scipy.io import readsav
import numpy as np
import time

np.set_printoptions(linewidth=225)

# ------------------------------------------------------------------ #
#  Paths (edit these as needed)
# ------------------------------------------------------------------ #
data_file   = 'examples/C-Mod/input/cmod_test_in.sav'
config_path = 'examples/C-Mod/python_output/config.json'
output_file = 'examples/C-Mod/python_output'  # outputs saved in examples/C-Mod/python_output/

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
