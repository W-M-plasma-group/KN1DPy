import numpy as np
import os.path

from .create_jh_bscoef import create_jh_bscoef
from KN1DPy.common.JH_Coef import JH_Coef

# This function generates JH_Coefficients and stores them in the JH_Coef object passed in.
# It also generates jh_bscoef.npz if it hasn't already been generated
# If both the coefficients and the file have been generated, this method is ignored
def generate_jh_coeffs(jh_coeffs : JH_Coef, create = 0) -> None:

    if create or not os.path.exists('jh_bscoef.npz'):
        create_jh_bscoef()

    if jh_coeffs.LogR_BSCoef is None:
        # this is where old data is restored 
        s = np.load('jh_bscoef.npz')
        # update global vars JH_coef common block
        jh_coeffs.DKnot = s['DKnot']
        jh_coeffs.TKnot = s['TKnot']
        jh_coeffs.order = s['order']
        jh_coeffs.LogR_BSCoef = s['LogR_BSCoef']
        jh_coeffs.LogS_BSCoef = s['LogS_BSCoef']
        jh_coeffs.LogAlpha_BSCoef = s['LogAlpha_BSCoef']
        jh_coeffs.A_Lyman = s['A_Lyman']
        jh_coeffs.A_Balmer = s['A_Balmer']
        
    return