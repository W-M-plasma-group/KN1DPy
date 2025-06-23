import numpy as np
from numpy.typing import NDArray

from .nh_saha import nh_saha
from .jhr_coef import jhr_coef # fixed capitalization - GG
from .generate_jh_coefficients import generate_jh_coeffs
from KN1DPy.common.JH_Coef import JH_Coef

#   Computes Balmer-alpha emissivity (watts m^-3) given the local
# electron density, electron temperature, and ground-state
# neutral density.

# Method :
#   (1) Compute the local n=3 population density using the Johnson-Hinnov
#       rate equations and coefficients [L.C.Johnson and E. Hinnov, J. Quant. 
#       Spectrosc. Radiat. Transfer. vol. 13 pp.333-358]
#   (2) Multiply by the n=3->2 spontaneous emission coefficient
#   (3) Convert to watts/m^3

def balmer_alpha(Density : NDArray, Te : NDArray, N0 : NDArray, jh_coeffs : JH_Coef, photons = 0, create = 0, no_null = 0):
    #________________________________________________________________________________
    # Input:
    #  	Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
    # 	Te	- fltarr, electron temperature (eV
    # 	N0	- fltarr, ground state neutral density (m^-3)
    #
    # Keywords:
    #	photons - returns emissivity in number of photons m^-3 s^-1
    #	create	- if set, then create bi-cubic spline coefficients for
    #		  interpolation of r0(p) r1(p) and save them in the
    #		  default save set. 
    #	No_Null	- if set, then rather than generate a NULL value when Density and Te
    #                 are not null but still outside the data range, compute the rate based on the min or max
    #		  data range values.
    #________________________________________________________________________________
    # History:
    #    Coding by B. LaBombard  6/29/99
    #    Coefficients from J. Terry's idl code JH_RATES.PRO
    
    generate_jh_coeffs(jh_coeffs, create)

    # From Johnson-Hinnov, eq (11):
    # n(3) = ( r(0) + r1(3) * n(1) / NHsaha(1) ) * NHsaha(3)

    if np.size(Density) != np.size(Te):
        raise Exception('Number of elements of Density and Te are different!')
    if np.size(Density) != np.size(N0):
        raise Exception(' Number of elements of Density and N0 are different! ')
    result = np.full(Density.shape,1.0e32)
    photons = np.full(Density.shape,1.0e32)
    r03 = jhr_coef(Density, Te, 0, 3, jh_coeffs, no_null = no_null)
    r13 = jhr_coef(Density, Te, 1, 3, jh_coeffs, no_null = no_null)
    NHSaha1 = nh_saha(Density, Te, 1)
    NHSaha3 = nh_saha(Density, Te, 3)
    ok=np.array([])
    for i in range(0, np.size(Density)):
        if 0 < N0[i] < 1e32 and r03[i] < 1.0e32 and r13[i] < 1.0e32 and NHSaha1[i] < 1.0e32 and NHSaha3[i] < 1.0e32:
            ok = np.append(ok, i)
    count = np.size(ok)
    if count > 0:
        for i in range(0, np.size(ok)):
            photons[i] = jh_coeffs.A_Balmer[0] * ( r03[i] + r13[i] * N0[i] / NHSaha1[i] ) * NHSaha3[i]
            result[i] = 13.6057 * ( 0.25 - 1.0 / 9.0 ) * photons[i] * 1.6e-19
    return result
