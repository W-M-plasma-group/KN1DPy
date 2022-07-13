import numpy as np
from poly import poly

#   Returns momentum transfer cross section for elastic collisions of H onto H2
# for specified energy of H. 
# 
# Data taken 
# Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
# Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 305.
#
# Input: 
# E     -   array or float, energy of H atom (target H2 molecules is at rest)
#
# Gwendolyn Galleher

def Sigma_El_H_HH(E): 
    if np.size(E) == 1:
        _E = float(E) # this line was in the original code I don't 
                      # know if its neccessary but it means E must have only 1 element
    _E = np.maximum(_E, 3.03e0)     
    _E = np.minimum(_E, 1.01e4)
    a = np.array[-3.495671e1, -4.062257e-1, -3.820531e-2, -9.404486e-3, 3.963723e-4]
    result = np.exp(poly(np.log(_E), a)) * 1e-4
    if np.ndim(E) == 0:
        result = result[0]
    return result