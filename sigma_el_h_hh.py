import numpy as np
from type_of.py import *
from poly.py import *

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
    ty = type_of(E) # I feel like there is a better way to do this is python I need to do some more thinking/research
    _E = float(E)
    _E = _E > 0.03e0 # I don't think these two lines are correct 
    _E = _E < 1.01e4
    a = np.array[-3.495671e1, -4.062257e-1, -3.820531e-2, -9.404486e-3, 3.963723e-4]
    result = np.exp(poly(np.log(_E), a)) * 1e-4
    if nDim == 0:
        result = result[0]
    return result