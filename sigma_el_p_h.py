import numpy as np 
from poly import poly 
#   Returns momentum transfer cross section or elastic collisions of H+ onto H 
# for specified energy of H+. 
# Data are taken from 
#; Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
# Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 298.
#
# Inputs: 
#  E       -   array or float, energy of H2 molecule (target H2 molecule is at rest)
# Gwendolyn Galleher
def sigma_el_p_h(E):
    E = np.array([E])
    _E = E.astype(float)
    # ensures that 0.001e0 < E < 1.01e5
    _E = np.maximum(_E, 0.001e0)
    _E = np.minimum(_E, 1.01e5)
    
    result = _E ; result= np.zeros(result.shape)

    ilow = np.argwhere(_E < 10.0)
    if np.size(ilow) > 0:  
        a = np.array([-3.233966e1, -1.126918e-1, 5.287706e-3, -2.445017e-3, -1.044156e-3, 8.419691e-5, 3.824773e-5])
        for i in ilow: 
            result[tuple(i)] = np.exp(poly(np.log(_E[tuple(i)]), a)) * 1e-4

    ihigh = np.argwhere(_E > 10.0)
    if np.size(ihigh) > 0:
        a = np.array([-3.231141e1, -1.386002e-1])
        for j in ihigh:
            result[tuple(j)] = np.exp(poly(np.log(_E[tuple(j)]), a)) * 1e-4
    return result
