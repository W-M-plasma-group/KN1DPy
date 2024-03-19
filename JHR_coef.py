import numpy as np
from scipy import interpolate
from create_jh_bscoef import Create_JH_BSCoef
import os.path
# Evaluates the r0(p) and r1(p) coefficients from Johnson-Hinnov tables 1 or 2.
# Gwendolyn Galleher 


def JHR_Coef(Density, Te, Ion, p, create = 0, no_null = 0, g=None):


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
    Dknot = g.JH_Coef_DKnot
    Tknot = g.JH_Coef_TKnot
    LogR_BSCoef=g.JH_Coef_LogR_BSCoef
    LogS_BSCoef=g.JH_Coef_LogS_BSCoef
    LogAlpha_BSCoef=g.JH_Coef_LogAlpha_BSCoef
    A_Lyman=g.JH_Coef_A_Lyman
    A_Balmer=g.JH_Coef_A_Balmer
    if create or not os.path.exists('jh_bscoef.npz'):
        Create_JH_BSCoef()
    if LogR_BSCoef is None:
        # this is where old data is restored 
        s=np.load('jh_bscoef.npz')
        Dknot=s['DKnot']
        Tknot=s['TKnot']
        order=s['order']
        LogR_BSCoef=s['LogR_BSCoef']
        LogS_BSCoef=s['LogS_BSCoef']
        LogAlpha_BSCoef=s['LogAlpha_BSCoef']
        A_Lyman=s['A_Lyman']
        A_Balmer=s['A_Balmer']
        # update global vars JH_coef common block
        g.JH_Coef_DKnot=Dknot
        g.JH_Coef_TKnot=Tknot
        g.JH_Coef_order=order
        g.JH_Coef_LogR_BSCoef=LogR_BSCoef
        g.JH_Coef_LogS_BSCoef=LogS_BSCoef
        g.JH_Coef_LogAlpha_BSCoef=LogAlpha_BSCoef
        g.JH_Coef_A_Lyman=A_Lyman
        g.JH_Coef_A_Balmer=A_Balmer
    
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
            # result[i] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LDensity[i],LTe[i],(Dknot,Tknot,LogS_BSCoef,3,3),0,0)) # updated
    return result 
