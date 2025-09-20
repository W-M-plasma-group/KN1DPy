import numpy as np
from numpy.typing import NDArray
from scipy import interpolate

from .generate_jh_coefficients import generate_jh_coeffs
from KN1DPy.common.JH_Coef import JH_Coef

# Evaluates the alpha coefficients 
# - the comment on the original file seemed to be wrong so I will probaly have to expand on this later
def jhalpha_coef(Density : NDArray, Te : NDArray, jh_coeffs : JH_Coef, create = 0, no_null = 0):
    # these are the variables that will be called from the classes for now I am defining them at the top. 
    # They are not defined as there actual values because the actual values used are defined in other files.
    
    generate_jh_coeffs(jh_coeffs, create)

    # Evaluate Alpha coefficients 
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
        for i in ok:
            #result[ok] = np.exp( ) # currently missing the python equivalent to bs2dr will come back to this later 
            result[i]=np.exp(interpolate.bisplev(LTe[i], LDensity[i], 
                                                 (jh_coeffs.TKnot, jh_coeffs.DKnot, jh_coeffs.LogAlpha_BSCoef, jh_coeffs.order, jh_coeffs.order), 0, 0)) # updated
            
    print("jhalpha_coef", result)
    input()

    return result
