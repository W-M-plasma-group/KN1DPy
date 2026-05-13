'''
New for the python version. Option to use the ADAS ionisation and recombination rates.
'''

import numpy as np
from scipy.interpolate import RectBivariateSpline
import os
import urllib.request

ADAS_DIR = os.path.dirname(os.path.abspath(__file__))

_ADAS_FILES = {
    'scd12_h.dat': 'https://open.adas.ac.uk/download/adf11/scd12/scd12_h.dat',
    'acd12_h.dat': 'https://open.adas.ac.uk/download/adf11/acd12/acd12_h.dat',
}

def _adas_path(filename):
    """Return the full path to an ADAS data file stored next to this module."""
    return os.path.join(ADAS_DIR, filename)

def _ensure_adas_data():
    for filename, url in _ADAS_FILES.items():
        path = _adas_path(filename)
        if not os.path.exists(path):
            print(f"Downloading ADAS data file: {filename}")
            urllib.request.urlretrieve(url, path)

def read_adf11(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    header_tokens = lines[0].split('/')[0].split()
    iz0    = int(header_tokens[0])
    n_ne   = int(header_tokens[1])
    n_Te   = int(header_tokens[2])
    is1min = int(header_tokens[3])
    is1max = int(header_tokens[4])

    def is_separator(line):
        return '/' in line

    def parse_numbers(line):
        return [float(x) for x in line.split()]

    sections = []
    current  = []
    for line in lines[1:]:
        if is_separator(line):
            sections.append(current)
            current = []
        else:
            stripped = line.strip()
            if stripped and not stripped.startswith('-' * 10):
                try:
                    current.extend(parse_numbers(line))
                except ValueError:
                    pass
    sections.append(current)

    grid_values = np.array(sections[0])
    ne_log10    = grid_values[:n_ne]
    Te_log10    = grid_values[n_ne:n_ne + n_Te]
    ne = 10.0 ** ne_log10
    Te = 10.0 ** Te_log10

    blocks   = []
    z1_vals  = []
    block_idx = 0
    for line in lines[1:]:
        if is_separator(line) and 'Z1=' in line.upper():
            z1_str = line.upper().split('Z1=')[1].split('/')[0].strip()
            z1_vals.append(int(z1_str))
            raw        = np.array(sections[block_idx + 1])
            rate_log10 = raw[:n_ne * n_Te].reshape((n_Te, n_ne)).T
            blocks.append(10.0 ** rate_log10)
            block_idx += 1

    return {'Te': Te, 'ne': ne, 'data': blocks, 'z1': z1_vals,
            'iz0': iz0, 'n_ne': n_ne, 'n_Te': n_Te}


def make_adf11_interpolator(filename, block=0):
    d    = read_adf11(filename)
    Te   = d['Te']
    ne   = d['ne']
    rate = d['data'][block]

    log_Te   = np.log10(Te)
    log_ne   = np.log10(ne)
    log_rate = np.log10(rate)
    spl      = RectBivariateSpline(log_ne, log_Te, log_rate, kx=3, ky=3)

    ne_lo, ne_hi = log_ne[0], log_ne[-1]
    Te_lo, Te_hi = log_Te[0], log_Te[-1]

    def interpolator(Te_eV, ne_cm3):
        Te_eV  = np.asarray(Te_eV,  dtype=float)
        ne_cm3 = np.asarray(ne_cm3, dtype=float)
        scalar = (Te_eV.ndim == 0 and ne_cm3.ndim == 0)
        lTe    = np.clip(np.log10(np.atleast_1d(Te_eV).ravel()),  Te_lo, Te_hi)
        lne    = np.clip(np.log10(np.atleast_1d(ne_cm3).ravel()), ne_lo, ne_hi)
        result = np.array([float(spl(lne[i], lTe[i]).item(0)) for i in range(lTe.size)])
        result = 10.0 ** result
        return float(result[0]) if scalar else result

    return interpolator


_ensure_adas_data()
_scd_interp = make_adf11_interpolator(_adas_path('scd12_h.dat'), block=0)
_acd_interp = make_adf11_interpolator(_adas_path('acd12_h.dat'), block=0)


def scd_adas(ne_m3, Te_eV):
    """
    Effective ionisation rate coefficient from ADAS SCD data.

    Parameters
    ----------
    ne_m3 : array-like, electron density in m^-3
    Te_eV : array-like, electron temperature in eV

    Returns
    -------
    ndarray, ionisation rate coefficient in m^3 s^-1
    """
    return _scd_interp(Te_eV, np.asarray(ne_m3) * 1e-6) * 1e-6


def acd_adas(ne_m3, Te_eV):
    """
    Effective recombination rate coefficient from ADAS ACD data.

    Parameters
    ----------
    ne_m3 : array-like, electron density in m^-3
    Te_eV : array-like, electron temperature in eV

    Returns
    -------
    ndarray, recombination rate coefficient in m^3 s^-1
    """
    return _acd_interp(Te_eV, np.asarray(ne_m3) * 1e-6) * 1e-6
