# Utility Functions for KN1DPy
import json
import tomllib
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import tomli_w
from numpy.typing import NDArray
from scipy import interpolate

# --- Configuration dataclasses ---

@dataclass
class KHConfig:
    '''Configuration settings for KN1DPy'''
    mesh_size: int = 10
    grid_fctr: float = 0.3
    ion_rate: str = 'adas'
    extra_energy_bins_eV: NDArray = field(default_factory=lambda: np.array([]))
    ci_test: bool = False
    alpha_cx_test: bool = False

@dataclass
class KHCollisions:
    '''Collision settings for Kinetic H procedure'''
    H2_H_EL: bool = False
    H_H_EL: bool = False
    H_P_EL: bool = False
    H_P_CX: bool = False
    SIMPLE_CX: bool = False

@dataclass
class KH2Collisions:
    '''Collision settings for Kinetic H2 procedure'''
    H2_H_EL: bool = False
    H2_H2_EL: bool = False
    H2_P_EL: bool = False
    H2_P_CX: bool = False
    SIMPLE_CX: bool = False


# --- File Paths ---

def get_local_directory(file_const) -> Path:
    return Path(file_const).resolve().parent


# --- Config Files ---

def get_json(file_path: Path | str) -> dict[str, Any]:
    with Path(file_path).open('r') as f:
        return json.load(f)

def _warn_json_deprecated(json_path: Path | str):
    abs_path = Path(json_path).resolve()
    filename  = abs_path.name
    toml_name = abs_path.with_suffix('.toml')
    warnings.warn(
        f"\n\nKN1DPy now uses TOML config files rather than JSON config files.\n"
        f"To convert '{filename}' to TOML, run the following Python commands from any directory:\n\n"
        f"    from KN1DPy.utils import convert_config_json_to_toml\n"
        f"    convert_config_json_to_toml(r'{abs_path}')\n\n"
        f"This will create '{toml_name}' in the same folder as '{filename}'.\n",
        DeprecationWarning,
        stacklevel=3,
    )

def get_config(config_path: str = './config.toml') -> dict[str, Any]:
    '''Load config file (.toml or .json). JSON files trigger a deprecation warning.'''
    config_path = Path(config_path)
    ext = config_path.suffix.lower()

    if ext == '.toml':
        if config_path.exists():
            with config_path.open('rb') as f:
                return tomllib.load(f)
        # Fall back to .json if .toml is missing
        json_path = config_path.with_suffix('.json')
        if json_path.exists():
            _warn_json_deprecated(json_path)
            return get_json(json_path)
        raise FileNotFoundError(f"Config file not found: {config_path}")

    elif ext == '.json':
        _warn_json_deprecated(config_path)
        return get_json(config_path)

    else:
        raise ValueError(f"Unsupported config format '{ext}'. Use .toml (recommended) or .json.")

def convert_config_json_to_toml(json_path: str) -> str:
    '''
    Convert a KN1DPy JSON config file to TOML format.
    The .toml file is saved in the same folder as the .json file.

    Parameters
    ----------
    json_path : str
        Path to the existing config .json file (absolute or relative).

    Returns
    -------
    str
        Absolute path of the newly created .toml file.

    Example
    -------
    Run from anywhere in Python:

        from KN1DPy.utils import convert_config_json_to_toml
        convert_config_json_to_toml('/path/to/your/config.json')
    '''
    abs_json = Path(json_path).resolve()
    abs_toml = abs_json.with_suffix('.toml')
    config = get_json(abs_json)
    with Path(abs_toml).open('wb') as f:
        tomli_w.dump(config, f)
    print(f"Converted '{Path(abs_json).name}' -> '{Path(abs_toml).name}'")
    print(f"Saved to:  {abs_toml}")
    return abs_toml


def convert_config_dict_to_dataclasses(config: dict[str, Any] | None = None) -> tuple[KHConfig, KHCollisions, KHConfig, KH2Collisions]:
    if not isinstance(config, dict):
        config = {}
    kh_cfg_dict = config.get('kinetic_h', {})
    kh2_cfg_dict = config.get('kinetic_h2', {})
    coll_dict = config.get('collisions', {})
    kh_cfg = KHConfig(
        kh_cfg_dict.get('mesh_size', 10),
        kh_cfg_dict.get('grid_fctr', 0.3),
        kh_cfg_dict.get('ion_rate', 'adas'),
        np.atleast_1d(kh_cfg_dict.get('extra_energy_bins_eV', [])),
        kh_cfg_dict.get('ci_test', False),
        kh_cfg_dict.get('alpha_cx_test', False),
    )
    kh_coll = KHCollisions(
        coll_dict.get('H2_H_EL', False),
        coll_dict.get('H_H_EL', False),
        coll_dict.get('H_P_EL', False),
        coll_dict.get('H_P_CX', False),
        coll_dict.get('SIMPLE_CX', False),
    )
    kh2_cfg = KHConfig(
        kh2_cfg_dict.get('mesh_size', 10),
        kh2_cfg_dict.get('grid_fctr', 0.3),
        kh2_cfg_dict.get('ion_rate', 'adas'),
        np.atleast_1d(kh2_cfg_dict.get('extra_energy_bins_eV', [])),
        kh2_cfg_dict.get('ci_test', False),
        kh2_cfg_dict.get('alpha_cx_test', False),
    )
    kh2_coll = KH2Collisions(
        coll_dict.get('H2_H_EL', False),
        coll_dict.get('H2_H2_EL', False),
        coll_dict.get('H2_P_EL', False),
        coll_dict.get('H2_P_CX', False),
        coll_dict.get('SIMPLE_CX', False),
    )
    return kh_cfg, kh_coll, kh2_cfg, kh2_coll

def convert_config_file_to_dataclasses(config_path: str = './config.toml') -> tuple[KHConfig, KHCollisions, KHConfig, KH2Collisions]:
    if not isinstance(config_path, str):
        config_path = str(config_path)
    config = get_config(config_path)
    return convert_config_dict_to_dataclasses(config)


# --- Printing ---

def debrief(statement: str, condition: bool):
    ''' Print statement if condition is true '''

    if condition:
        print(statement)

def sval(s,length=None):
  ''' removes leading / trailing spaces and truncates string to a specified length '''

  return str(s).strip()[:length]


# --- Bounding  ---

class Bound:
    '''
    Defines boundaries of some array, stores first and last index of bounding

    Attributes
    ----------
        start: int
            first index
        end: int
            last index
    '''

    def __init__(self, first : int, last : int):
        self.start = first
        self.end = last

    def range(self):
        '''Return inclusive range of values between start and end of bound'''
        return range(self.start, self.end+1)

    def slice(self, start_offset=0, end_offset=0):
        '''Returns slice object for bounds, with offset'''
        return slice(self.start+start_offset, self.end+end_offset)


# --- Polynomials ---

def poly(x, c):
    '''
    Evaluate a polynomial at one or more points

    Parameters
    ----------
        x : float or ndarray
            Variable/s to evaluate the polynomial at
        c : ndarray
            array of polynomial coefficients

    Returns
    -------
        y : float or ndarray
            Value of the polynomial evaluated at x, array of values if x is an array
    '''

    x = np.asarray(x)
    n = len(c)-1
    y = c[n]
    for i in range(n-1, -1, -1):
        y = y*x + c[i]
    return y


# --- Interpolation ---

def interp_1d(funx: NDArray, funy: NDArray, x: NDArray, kind: str = 'linear', axis: int = -1,
        copy: bool = True, bounds_error: Any | None = None, fill_value: float = np.nan, assume_sorted: bool = False):
    ''' Wrapper function for creating a scipy 1d interpolation function and run it on an array '''

    interpfunc = interpolate.interp1d(funx, funy, kind=kind, axis=axis, copy=copy, bounds_error=bounds_error, fill_value=fill_value, assume_sorted=assume_sorted)
    return interpfunc(x)

def path_interp_2d(p, px, py, x, y):
    interp = interpolate.RegularGridInterpolator((px, py), p, method='linear')
    points = np.column_stack([x, y])
    return interp(points)

def bs2dr(x, y, kx_ord, ky_ord, xknot, yknot, bscoef):
    '''
    IDL bs2dr translation equivalent
    '''
    try:
        import scipy.interpolate._dfitpack as _dfitpack
        return _dfitpack.bispeu(yknot, xknot, bscoef, kx_ord-1, ky_ord-1, y, x)[0]
    except ImportError:
        # Fallback for scipy versions without _dfitpack: evaluate point-by-point
        # using bisplev (grid evaluator), available in all scipy versions.
        tck = (yknot, xknot, bscoef, kx_ord-1, ky_ord-1)
        return np.array([interpolate.bisplev(np.array([y[i]]), np.array([x[i]]), tck).ravel()[0]
                         for i in range(len(x))])


# --- Table Searching ---

def locate(table, value):
    '''
    Finds the index of a value (or values) in a sorted table using np.searchsorted.

    Parameters
    ----------
        table : ndarray
            Sorted list or array of numbers (ascending or descending).
        value : float, ndarray
            Value(s) to search for

    Returns
    -------
        ndarray
            Array of indices (integers) corresponding to the positions where the values meet the conditions.
    '''

    # Convert inputs to NumPy arrays if they are scalars or lists
    table = np.asarray(table)
    value = np.atleast_1d(value)  # Ensure `value` is an array

    # Determine if the table is in ascending or descending order
    asc = table[0] <= table[-1]

    if not asc:
        # If the table is in descending order, temporarily reverse it
        table = table[::-1]

    # Use np.searchsorted to find the indices
    indices = np.searchsorted(table, value, side='right' if asc else 'left') - 1

    # Adjust indices for descending tables
    if not asc:
        indices = len(table) - indices - 1

    # Handle special cases: out-of-range values
    indices[value < table[0]] = -1  # Values less than the first element
    indices[value >= table[-1]] = len(table) - 1  # Values greater than or equal to the last element

    if(len(indices) == 1): #Convert to scalar if only one value
        indices = indices[0]

    return indices


# --- Reverse Function from reverse.pro ---

def reverse(a, subscript=1):
    '''
        reverses the order of a list at the given dimension (subscript)
        initially assume at least 1 dimension
    '''

    ndims = 1
    b = a

    #if the 1st variable is also a list then a dimension is added, recurring until no longer true
    while type(b[0]) is list:
        ndims += 1
        if len(b) == 0:
            break
        b = b[0]
    if subscript > ndims:
        raise Exception('Subscript_index must be less than or equal to number of dimensions.')
    if subscript == 1: #unique case where it is reversing the 1st dim
        a = a[::-1]
        return a
    return rev_rec(a, subscript, 1)

def rev_rec(a, subscript, dim_tracker):
    ''' Recursive function that iterates over everything in a, and reverses everything in the specified dim '''

    i = 0
    while i < len(a):
        if dim_tracker == subscript-1:
            a[i] = a[i][::-1]
        else:
            a[i] = rev_rec(a[i], subscript, dim_tracker+1)
        i += 1
    return a
