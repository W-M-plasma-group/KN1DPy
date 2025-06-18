import numpy as np
from scipy import interpolate
import os.path

from .create_jh_bscoef import create_jh_bscoef

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
def jhs_coef(Density, Te, create = 0, no_null = 0, g=None):
    # these are the variables that will be called from the classes for now I am defining them at the top. 
    # They are not defined as there actual values because the actual values used are defined in other files.
    Dknot = g.JH_Coef_DKnot
    Tknot = g.JH_Coef_TKnot
    LogR_BSCoef=g.JH_Coef_LogR_BSCoef
    LogS_BSCoef=g.JH_Coef_LogS_BSCoef
    LogAlpha_BSCoef=g.JH_Coef_LogAlpha_BSCoef
    A_Lyman=g.JH_Coef_A_Lyman
    A_Balmer=g.JH_Coef_A_Balmer
    if create or not os.path.exists('jh_bscoef.npz'):
        create_jh_bscoef()
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

    # Evaluate S coefficients 
    if np.size(Density) != np.size(Te):
        raise Exception('Number of elements of Density and Te are different!')
    result = Density ; result[:] = 1.0e32
    LDensity = np.log(Density)
    LTe = np.log(Te)
    if no_null:
        LDensity = np.maximum(LDensity, min(Dknot))
        LDensity = np.minimum(LDensity, max(Dknot))
        LTe = np.maximum(LTe, min(Tknot))
        LTe = np.minimum(LTe, max(Tknot))
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
            #result[i] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LDensity[i],LTe[i],(Dknot,Tknot,LogS_BSCoef,3,3),0,0)) # updated 
    return result
