import numpy as np
import scipy.interpolate

from .sval import sval
from .make_sigmav import make_sigmav
from .sigma.sigma_el_p_h import sigma_el_p_h # added import for sigma function

#Creates a save set storing Bi-cubic spline interpolation
def create_sigmav_el_h_p_bscoef():

    #particle is h, target is p

    sigma_function = sigma_el_p_h # changed string variable to function

    mE=6
    Emin=0.1
    Emax=2.0e4
    nT=6
    Tmin=0.1
    Tmax=2.0e4
    findgen1 = [None] * mE
    for i in range(mE):
        findgen1[i] = float(i)
    findgen2 = [None] * nT
    for j in range(nT):
        findgen2[j] = float(j)
    E_Particle=10^(np.log10(Emin)+(np.log10(Emax)-np.log10(Emin))*findgen1/(mE-1))
    mu_particle=1.0
    T_target=10^(np.log10(Tmin)+(np.log10(Tmax)-np.log10(Tmin))*findgen2/(nT-1))
    mu_target=1.0
    SigmaV=np.zeros((nT,mE)) # changed SigmaV creation method
    for iT in range(len(T_target)): # fixed range
        print('Processing T='+sval(T_target(iT)))
        SigmaV[iT,:]=make_sigmav(E_Particle,mu_particle,T_target[iT],mu_target,sigma_function) # fixed SigmaV indexing

    print('Computing B-Spline coefficients')
    order_EL_H_P=4
    order=order_EL_H_P
    LogE_Particle=np.alog(E_Particle)
    LogT_Target=np.alog(T_target)
    LogSigmaV_EL_H_P=np.alog(SigmaV)
    LogSigmav_Interp = scipy.interpolate.RectBivariateSpline(LogE_Particle, LogT_Target, LogSigmaV_EL_H_P) 
    LogSigmaV_EL_H_P_BSCoef = LogSigmav_Interp.get_coeffs()
    EKnot_EL_H_P,TKnot_EL_H_P=LogSigmav_Interp.get_knots() 

    np.savez('sigmav_el_h_p_bscoef',
             EKnot_EL_H_P=EKnot_EL_H_P,
             TKnot_EL_H_P=TKnot_EL_H_P,
             order_EL_H_P=order_EL_H_P,
             LogSigmaV_EL_H_P_BSCoef=LogSigmaV_EL_H_P_BSCoef) # now saves to current directory

    return
