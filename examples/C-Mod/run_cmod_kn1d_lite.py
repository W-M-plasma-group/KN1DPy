"""
C-Mod kn1d_lite example
========================
Runs kn1d_lite on the closed-field lines part of the C-Mod example
profile, using a mono-energetic incident neutral boundary condition.
"""

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
import numpy as np
from scipy.io import readsav
from KN1DPy.kn1d_lite import kn1d_lite

#  Settings
data_file   = 'examples/C-Mod/input/cmod_test_in.sav'
config_path = 'config.toml'
incident_n0 = 1e15      # m^-3 — incident neutral density at separatrix
energy_eV   = 3.0       # eV   — incident neutral energy

#  Load and prepare profile
sav = readsav(data_file)

# Unit conversions
Ti  = sav['ti']  * 1e3    # keV -> eV
Te  = sav['te']  * 1e3    # keV -> eV
n   = sav['n']   * 1e20   # 1e20 m^-3 -> m^-3
x   = sav['x']            # m, x=0 at wall, increases into plasma
vxi = sav['vxi']          # m/s. Should all be zeros anyway
mu  = int(sav['mu'])      # ion mass. Just 2 for deuterium
xsep = float(sav['xsep']) # separatrix position (m)

# Adapt inputs for kn1d_lite:
# 1) slice to only include points inside the separatrix (x >= xsep)
# 2) remap x so that x=0 at the separatrix and x increases into the plasma
# 3) add in another point at exactly x=0, because this is where the incident neutrals BC is applied.
mask  = x >= xsep
x_in   = np.concatenate([[0.0],  x[mask] - xsep])
Ti_in  = np.concatenate([[np.interp(xsep, x, Ti)],  Ti[mask]])
Te_in  = np.concatenate([[np.interp(xsep, x, Te)],  Te[mask]])
n_in   = np.concatenate([[np.interp(xsep, x, n)],   n[mask]])
vxi_in = np.concatenate([[np.interp(xsep, x, vxi)], vxi[mask]])

# ------------------------------------------------------------------ #
#  Run kn1d_lite
# ------------------------------------------------------------------ #
print(f'\nRunning kn1d_lite  (n0={incident_n0:.1e} m^-3, E_in={energy_eV} eV) ...')

result = kn1d_lite(
    x=x_in, mu=mu, Ti=Ti_in, Te=Te_in, n=n_in,
    vxi=vxi_in,
    incident_n0=incident_n0,
    energies_eV=[energy_eV],
    config_path=config_path,
)

print('Done.\n')