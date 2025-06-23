import numpy as np
from numpy.typing import NDArray
from scipy import interpolate
import os.path

from .generate_jh_coefficients import generate_jh_coeffs
from KN1DPy.common.JH_Coef import JH_Coef

# Evaluates the ionization rate coefficients, S (m^-3 s^-1), from  Johnson-Hinnov table 2 (MKS units).
#; Input:
#  	Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
#  	Te	- fltarr, electron temperature (eV)
#
# Keywords:
#	create	- if set, then create bi-cubic spline coefficients for
#		  interpolation of S and save them in the default save set. 
#	No_Null	- if set, then rather than generate a NULL value when Density and Te
#                 are outside the data range, compute the rate based on the min or max
#		  data range values.
#________________________________________________________________________________
# History:
#    Coding by B. LaBombard  6/29/99
#    Coefficients from J. Terry's idl code JH_RATES.PRO
def jhs_coef(Density : NDArray, Te : NDArray, jh_coeffs : JH_Coef, create = 0, no_null = 0) -> NDArray:

    generate_jh_coeffs(jh_coeffs, create)

    # Evaluate S coefficients 
    if np.size(Density) != np.size(Te):
        raise Exception('Number of elements of Density and Te are different!')
    result = Density ; result[:] = 1.0e32
    LDensity = np.log(Density)
    LTe = np.log(Te)
    if no_null:
        LDensity = np.maximum(LDensity, min(jh_coeffs.DKnot))
        LDensity = np.minimum(LDensity, max(jh_coeffs.DKnot))
        LTe = np.maximum(LTe, min(jh_coeffs.TKnot))
        LTe = np.minimum(LTe, max(jh_coeffs.TKnot))
        count = np.size(LDensity)
        ok = np.arange(count)
    else:
        for i in range(0, len(Density)):
            if min(jh_coeffs.DKnot) < LDensity[i] < max(jh_coeffs.DKnot) and min(jh_coeffs.TKnot) < LTe[i] < max(jh_coeffs.TKnot):
                ok = np.append(ok, i)
    count = np.size(ok)
    ok = ok.astype(int) 
    if count > 0: 
        for i in ok: # fixed how result is defined not completely confident in this - GG
            #result[i] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LDensity[i], LTe[i], 
                                                 (jh_coeffs.DKnot, jh_coeffs.TKnot, jh_coeffs.LogS_BSCoef, 3, 3), 0, 0)) # updated 
    return result
