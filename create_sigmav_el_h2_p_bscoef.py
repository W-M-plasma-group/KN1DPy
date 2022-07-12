import numpy as np
import scipy.interpolate
from sval import sval
from make_sigmav import make_sigmav
#   Creates save set storing Bi-cubic spine interpoaltion 
# coefficients for parameters in Johnson-Hinov rate equations.
# Gwendolyn Galleher
def Create_Sigmav_EL_H2_P_BSCOEF():
    Sigma_Function = 'sigma_el_P_HH'

    # particle is HH
    # target is P 

    mE = 5
    Emin = 0.1 ; Emax = 1.03e3
    nT=5
    Tmin=0.1 ; Tmax=1.0e3
    E_Particle = 10 ** (np.log10(Emin) + (np.log10(Emax) - np.log10(Emin)) * np.arange(mE)/(mE-1))
    Mu_Particle = 2.0 
    T_Target = 10 ** (np.log10(Tmin) + (np.log10(Tmax) - np.log10(Tmin)) * np.arange(nT)/(nT-1))
    Mu_Target = 1.0 
    SigmaV = np.zeros(nT, mE)
    for i in range(0, np.size(T_Target)-1):
        print('Processing T=' + sval(T_Target[i]))
        SigmaV[(0,i)] = make_sigmav(E_Particle,Mu_Particle,T_Target[i],Mu_Target,Sigma_Function)


    print('Computing B-Spline coefficients')
    LogE_Particle=np.log(E_Particle)
    LogT_Target=np.log(T_Target)
    LogSigmaV_EL_H2_P = np.log(SigmaV)
    LogSigmav_Interp = scipy.interpolate.RectBivariateSpline(LogE_Particle, LogT_Target, LogSigmaV_EL_H2_P) 
    LogSigmaV_EL_H2_P_BSCoef = LogSigmav_Interp.get_coeffs()
    print('Saving results in files: ')
    # We havent discussed how we are saving things yet so I left the last lines out 
    return