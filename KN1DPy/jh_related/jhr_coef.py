import numpy as np
from numpy.typing import NDArray
from scipy import interpolate

from .generate_jh_coefficients import generate_jh_coeffs
from KN1DPy.common.JH_Coef import JH_Coef

# Evaluates the r0(p) and r1(p) coefficients from Johnson-Hinnov tables 1 or 2.
# Gwendolyn Galleher 


def jhr_coef(Density : NDArray, Te : NDArray, Ion : int, p : int, jh_coeffs : JH_Coef, create = 0, no_null = 0):


# Input:
#  	Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
#  	Te	- fltarr, electron temperature (eV
#  	Ion	- integer, =0: return "recombination" coeffcient, r0(p)
#		   =1: return "ionization" coeffcient, r1(p)
#  	p	- integer, hydrogen energy level, p=1 is ground state
# Key Words: 
#	create	- if set, then create bi-cubic spline coefficients for
#		  interpolation of r0(p) r1(p) and save them in the
#		  default save set. 
#	No_Null	- if set, then rather than generate a NULL value when Density and Te
#                 are outside the data range, compute the rate based on the min or max
#		  data range values.

    generate_jh_coeffs(jh_coeffs, create)
    
    # Evaluates R coefficients 
    if np.size(Density) != np.size(Te):
        raise Exception('Number of elements of Density and Te are different!')
    if np.size(Ion) != 1:
        raise Exception('"Ion" must be a scalar')
    if np.size(p) != 1:
        raise Exception('"p" must be a scalar')
    if p < 2 or p > 6:
        raise Exception('"p" must in range 1 < p < 7')
    if Ion < 0 or Ion > 1:
        raise Exception('"Ion" must 0 or 1')
    result = Density ; result[:] = 1.0e32
    LDensity = np.log(Density)
    LTe = np.log(Te)
    if no_null:
        LDensity = np.maximum( LDensity, min(jh_coeffs.DKnot))
        LDensity = np.minimum( LDensity, max(jh_coeffs.DKnot))
        LTe = np.maximum( LTe, min(jh_coeffs.TKnot))
        LTe = np.minimum( LTe, max(jh_coeffs.TKnot))
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
            # result[i] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LTe[i], LDensity[i],
                                                 (jh_coeffs.TKnot, jh_coeffs.DKnot, jh_coeffs.LogR_BSCoef[:,Ion,p-2], jh_coeffs.order, jh_coeffs.order),0,0)) # updated
            
    print("jhr_coef", result)
    input()
    return result 
