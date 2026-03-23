# Utility Functions for KN1DPy
from __future__ import annotations

import json
from typing import Any
import os

from numpy.typing import NDArray
import numpy as np
from scipy import interpolate
from scipy.io import readsav

# --- File Paths ---

def get_local_directory(file_const):
    dir_path = os.path.dirname(os.path.realpath(file_const))
    return dir_path


# --- Json Files ---

def get_json(file_path:str) -> dict[str, Any]:
    ''' Load json file '''

    with open(file_path, 'r') as config:
        return json.load(config)
    
def get_config(config_path: str = './config.json') -> dict[str, Any]:
    ''' Load config file from the given path (defaults to ./config.json) '''

    return get_json(config_path)


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
    while type(b[0]) == list:
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


# --- Read Functions ---

def sav_read(sav_path, nc_path):
    '''
    Used to read and save .sav files

    Parameters
    ----------
        sav_path : str
            Path to the .sav input file
        nc_path : str
            Path to the .nc file being created
            
    Returns
    -------
        input_dict : dict
            Dictionary of all inputs from the input file
    '''

    import netCDF4 as nc
    sav_data = readsav(sav_path)
    fn = nc_path
    ds = nc.Dataset(fn, 'w', format = 'NETCDF4') 
    for k,v in sav_data.items(): 
        setattr(ds, k, v)
    input_dict = ds.__dict__
    return input_dict

def nc_read(nc_path):
    '''
    Used to read and save .nc files (netCDF)

    Parameters
    ----------
        nc_path : str
            Path to the .nc file being created
            
    Returns
    -------
        input_dict : dict
            Dictionary of all inputs from the input file
    '''
    
    import netCDF4 as nc
    fn = nc_path
    ds = nc.Dataset(fn)
    input_dict = ds.__dict__
    return input_dict