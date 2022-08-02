import numpy as np
import scipy 
from create_jh_bscoef import create_jh_bscoef

# Evaluates the r0(p) and r1(p) coefficients from Johnson-Hinnov tables 1 or 2.
# Gwendolyn Galleher 


def JHR_Coef(Density, Te, Ion, p, create = 0, no_null = 0):


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

# variables in common block 
    Dknot = None
    Tknot = None
    LogR_BSCoef=None
    LogS_BSCoef=None
    LogAlpha_BSCoef=None
    A_Lyman=None
    A_Balmer=None
    if create:
        create_jh_bscoef()
    if LogR_BSCoef is None:
        # this is where old data is restored I don't entirely know how we want to do this yet or if we are doing this at all 
        pass
    
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
        LDensity = np.maximum( LDensity, min(Dknot))
        LDensity = np.minimum( LDensity, max(Dknot))
        LTe = np.maximum( LTe, min(Tknot))
        LTe = np.minimum( LTe, Tknot)
        count = np.size(LDensity)
        ok = np.arange(count)
    else:
         for i in range(0, len(Density)):
            if min(Dknot) < LDensity[i] < max(Dknot) and min(Tknot) < LTe[i] < max(Tknot):
                ok = np.append(ok, i)
    count = np.size(ok)
    ok = ok.astype(int) 
    if count > 0: 
        for i in ok: # fixed how result is defined not completely confident in this - GG
            result[i] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
    return result 