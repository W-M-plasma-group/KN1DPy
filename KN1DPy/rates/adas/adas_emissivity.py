'''ADAS PEC-based Lyman-alpha and Balmer-alpha emissivity for atomic hydrogen.'''

import re
import urllib.request
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

ADAS_DIR = Path(__file__).resolve().parent

_PEC_FILENAME = 'pec12_h_pju_h0.dat'
_PEC_URL = 'https://open.adas.ac.uk/download/adf15/pec12][h/pec12][h_pju][h0.dat'

# ISEL indices within pec12][h_pju][h0.dat
_ISEL_LYA_EXC = 1    # Lyman-alpha 1215.2 A, excitation
_ISEL_LYA_REC = 67   # Lyman-alpha 1215.2 A, recombination
_ISEL_BAL_EXC = 3    # Balmer-alpha 6561.9 A, excitation
_ISEL_BAL_REC = 69   # Balmer-alpha 6561.9 A, recombination

# Photon energies (J): E = hc/lambda
_E_LYA = 6.626e-34 * 2.998e8 / 1215.2e-10  # 1.635e-18 J  (10.2 eV)
_E_BAL = 6.626e-34 * 2.998e8 / 6561.9e-10  # 3.027e-19 J  (1.89 eV)


def _ensure_pec_data():
    path = ADAS_DIR / _PEC_FILENAME
    if not path.exists():
        print(f'Downloading ADAS PEC data: {_PEC_FILENAME}')
        urllib.request.urlretrieve(_PEC_URL, path)
    return path


def _read_adf15_block(lines, target_isel):
    '''Extract ne, Te grid and PEC values for a given ISEL from ADF15 file lines.

    Returns ne (m⁻³), Te (eV), pec (m³/s) with shape (n_ne, n_Te).
    '''
    header_idx = None
    n_ne = n_Te = None
    for i, line in enumerate(lines):
        m = re.search(r'/ISEL\s*=\s*(\d+)', line, re.IGNORECASE)
        if m and int(m.group(1)) == target_isel:
            parts = line.split()
            n_ne = int(parts[1])
            n_Te = int(parts[2])
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f'ISEL={target_isel} not found in ADF15 file')

    needed = n_ne + n_Te + n_ne * n_Te
    vals = []
    for line in lines[header_idx + 1:]:
        if 'ISEL' in line.upper():
            break  # reached the next block header
        try:
            vals.extend(float(x) for x in line.split())
        except ValueError:
            continue
        if len(vals) >= needed:
            break

    ne_m3  = np.array(vals[:n_ne]) * 1e6                           # cm⁻³ → m⁻³
    Te_eV  = np.array(vals[n_ne : n_ne + n_Te])
    pec_m3 = np.array(vals[n_ne + n_Te : needed]).reshape(n_ne, n_Te) * 1e-6  # cm³/s → m³/s

    return ne_m3, Te_eV, pec_m3


def _make_pec_interp(ne, Te, pec):
    '''Log-log RectBivariateSpline for PEC(ne, Te), clipped to grid bounds.

    ne in m⁻³, Te in eV, pec in m³/s.  Returns a callable (ne_m3, Te_eV) → pec_m3.
    '''
    log_ne  = np.log10(ne)
    log_Te  = np.log10(Te)
    log_pec = np.log10(np.maximum(pec, 1e-99))
    spl = RectBivariateSpline(log_ne, log_Te, log_pec, kx=3, ky=3)
    ne_lo, ne_hi = log_ne[0], log_ne[-1]
    Te_lo, Te_hi = log_Te[0], log_Te[-1]

    def interp(ne_m3, Te_eV):
        lne = np.clip(np.log10(np.asarray(ne_m3, dtype=float).ravel()), ne_lo, ne_hi)
        lTe = np.clip(np.log10(np.asarray(Te_eV, dtype=float).ravel()), Te_lo, Te_hi)
        return 10.0 ** np.array([float(spl(lne[i], lTe[i]).item(0))
                                  for i in range(lne.size)])

    return interp


def _load_interpolators():
    path = _ensure_pec_data()
    with path.open() as f:
        lines = f.readlines()

    ne_e, Te_e, pec_lya_exc = _read_adf15_block(lines, _ISEL_LYA_EXC)
    ne_r, Te_r, pec_lya_rec = _read_adf15_block(lines, _ISEL_LYA_REC)
    _,    _,    pec_bal_exc = _read_adf15_block(lines, _ISEL_BAL_EXC)
    _,    _,    pec_bal_rec = _read_adf15_block(lines, _ISEL_BAL_REC)

    return (
        _make_pec_interp(ne_e, Te_e, pec_lya_exc),
        _make_pec_interp(ne_r, Te_r, pec_lya_rec),
        _make_pec_interp(ne_e, Te_e, pec_bal_exc),
        _make_pec_interp(ne_r, Te_r, pec_bal_rec),
    )


try:
    _lya_exc, _lya_rec, _bal_exc, _bal_rec = _load_interpolators()
    available = True
except Exception:
    available = False


def lyman_alpha(ne, Te, n0, ni=None):
    '''
    Lyman-alpha emissivity (W/m³) using ADAS PEC coefficients.

    emissivity = (ne * n0 * PEC_exc + ne * ni * PEC_rec) * E_photon

    Parameters
    ----------
    ne : electron density (m⁻³)
    Te : electron temperature (eV)
    n0 : neutral hydrogen density (m⁻³)
    ni : ion density (m⁻³); defaults to ne (quasi-neutral pure hydrogen)
    '''
    ne = np.asarray(ne, dtype=float).ravel()
    Te = np.asarray(Te, dtype=float).ravel()
    n0 = np.asarray(n0, dtype=float).ravel()
    ni = ne if ni is None else np.asarray(ni, dtype=float).ravel()
    photons = ne * (n0 * _lya_exc(ne, Te) + ni * _lya_rec(ne, Te))
    return photons * _E_LYA


def balmer_alpha(ne, Te, n0, ni=None):
    '''
    Balmer-alpha emissivity (W/m³) using ADAS PEC coefficients.

    emissivity = (ne * n0 * PEC_exc + ne * ni * PEC_rec) * E_photon

    Parameters
    ----------
    ne : electron density (m⁻³)
    Te : electron temperature (eV)
    n0 : neutral hydrogen density (m⁻³)
    ni : ion density (m⁻³); defaults to ne (quasi-neutral pure hydrogen)
    '''
    ne = np.asarray(ne, dtype=float).ravel()
    Te = np.asarray(Te, dtype=float).ravel()
    n0 = np.asarray(n0, dtype=float).ravel()
    ni = ne if ni is None else np.asarray(ni, dtype=float).ravel()
    photons = ne * (n0 * _bal_exc(ne, Te) + ni * _bal_rec(ne, Te))
    return photons * _E_BAL
