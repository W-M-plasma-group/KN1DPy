import numpy as np 

#   Evaluates the elastic scattering (momentum transfer) <sigma v> for
#  H2  molecules with energy E impacting on protons with temperature T.

#  Data generated from cross sections tablulated in:
# Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 -
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and
# Helium Ions, Neutrals, and Isotopes  by D.R. Schultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 305.

def sigmav_el_hh_p(T, E, use_Bspline = 0):
    #  Input:
    #	T	- fltarr(*) or float, proton temperature (eV)
    #	E	- fltarr(*) or float, atom mono-energy (eV)
    #
    #  Output:
    #	returns <sigma V> for 0.1 < T < 1e3 and 0.1 < E < 1e5 
    #	units: m^3/s
    #       if T and/or E is outside this range, the value on the boundary is returned
    if np.size(E) != np.size(T):
        raise Exception('Number of elements of E and T are different!')
    _E = np.maximum(E, 0.0001)
    _T = np.maximum(T, 0.0001)
    _E = np.minimum(_E, 1.0e5)
    _T = np.minimum(_T, 1.0e5)
    LEP = np.log(_E)
    LTH2 = np.log(_T)
    if use_Bspline:
        # JSigmaV_EL_HH_P common block - will change later 
        EKnot_EL_H2_P=None
        TKnot_EL_H2_P=None
        order_EL_H2_P=None
        LogSigmaV_EL_H2_P_BSCoef=None
        if LogSigmaV_EL_H2_P_BSCoef is None:
            # haven't talked about how to restore yet
            pass 
        LEP = np.maximum(LEP, min(EKnot_EL_H2_P))
        LEP = np.minimum(LEP, max(EKnot_EL_H2_P))
        LTH2 = np.maximum(LTH2, min(TKnot_EL_H2_P))
        LTH2 = np.minimum(LTH2, max(TKnot_EL_H2_P))
        result = np.exp( ) # still haven't figured out what to do with BS2DR
    else:
        # SIGMAV_EL_H2_P_DATA common block - will change later 
        Ln_E_Particle=None
        Ln_T_Target=None
        SigmaV=None
        nEP=None
        nT=None
        if Ln_E_Particle is None:
            # haven't talked about how to restore yet
            pass
            nEP = np.size(Ln_E_Particle) - 1
            nT = np.size(Ln_T_Target) - 1
        LEP = np.maximum(LEP, min(Ln_E_Particle))
        LEP = np.minimum(LEP, max(Ln_E_Particle))
        LTH2 = np.maximum(LTH2, min(Ln_T_Target))
        LTH2 = np.minimum(LTH2, max(Ln_T_Target))
        iE = (LEP-Ln_E_Particle[0]) * nEP / (Ln_E_Particle[nEP] - Ln_E_Particle[0])
        iT = (LTH2-Ln_T_Target[0]) * nT / (Ln_T_Target[nT] - Ln_T_Target[0])
        result = np.interp(SigmaV, iE, iT)
    return result 
        