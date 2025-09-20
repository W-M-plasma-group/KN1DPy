import numpy as np

from .poly import poly

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

def sigma_el_h_hh(E): 
    E = np.array([E])
    _E = E.astype(float) 
    # ensures 0.03e0 < _E < 1.01e4
    _E = np.maximum(_E, 0.03e0)     
    _E = np.minimum(_E, 1.01e4)
    a = np.array([-3.495671e1, -4.062257e-1, -3.820531e-2, -9.404486e-3, 3.963723e-4]) # added missing parantheses - GG
    result = np.exp(poly(np.log(_E), a)) * 1e-4
    # deleted redundant if statement - GG
    return result