import numpy as np
from poly import poly
#   Returns momentum transfer cross section for elastic collisions of H2 onto H2
# for specififed energy of H2
# Data are taken from 
# 
#  Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
# Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 305.
# Input:
#   E       -   array or float, energy of H2 molecule (target H2 molecule is at rest)
#   VIS     -   if set, then return viscosity cross section instead of momentum 
#               transfer cross section 
# Gwendolyn Galleher 

def Sigma_EL_HH_HH( E, vis = 0):
    E = np.array([E])
    _E = E.astype(float)
    # ensures 3.03e0 < _E < 1.01e4
    _E = np.maximum(_E, 3.03e0)     
    _E = np.minimum(_E, 1.01e4)
    if vis: 
        # calculates viscosity cross section
        a = np.array([-3.430345e1, -2.960406e-1, -6.382532e-2, -7.557519e-3, 2.606259e-4])
        result = np.exp(poly(np.log(_E), a)) * 1e-4
    else: 
        # calculates momentum transfer cross section 
        a = np.array([-3.430345e1, -2.960406e-1, -6.382532e-2, -7.557519e-3, 2.606259e-4])
        result = np.exp(poly(np.log(_E), a)) * 1e-4 
    return result