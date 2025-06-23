import numpy as np
from numpy.typing import NDArray
from scipy import interpolate
import os.path

from .create_jh_bscoef import create_jh_bscoef
from .common.JH_Coef import JH_Coef

# Evaluates the alpha coefficients 
# - the comment on the original file seemed to be wrong so I will probaly have to expand on this later
def jhalpha_coef(Density : NDArray, Te : NDArray, jh_coeffs : JH_Coef, create = 0, no_null = 0):
    # these are the variables that will be called from the classes for now I am defining them at the top. 
    # They are not defined as there actual values because the actual values used are defined in other files.
    Dknot = jh_coeffs.DKnot
    Tknot = jh_coeffs.TKnot
    LogR_BSCoef = jh_coeffs.LogR_BSCoef
    LogS_BSCoef = jh_coeffs.LogS_BSCoef
    LogAlpha_BSCoef = jh_coeffs.LogAlpha_BSCoef
    A_Lyman = jh_coeffs.A_Lyman
    A_Balmer = jh_coeffs.A_Balmer
    if create or not os.path.exists('jh_bscoef.npz'):
        create_jh_bscoef()
    print("running")
    if LogR_BSCoef is None:
        print("running 2")
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
        jh_coeffs.DKnot=Dknot
        jh_coeffs.TKnot=Tknot
        jh_coeffs.order=order
        jh_coeffs.LogR_BSCoef=LogR_BSCoef
        jh_coeffs.LogS_BSCoef=LogS_BSCoef
        jh_coeffs.LogAlpha_BSCoef=LogAlpha_BSCoef
        jh_coeffs.A_Lyman=A_Lyman
        jh_coeffs.A_Balmer=A_Balmer

    # Evaluate Alpha coefficients 
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
        for i in ok:
            #result[ok] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LDensity[i],LTe[i],(Dknot,Tknot,LogS_BSCoef,3,3),0,0)) # updated
    return result
