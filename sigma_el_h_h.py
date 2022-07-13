import numpy as np 
from poly import poly

#   Returns momentum transfer cross section for elastic collisions of H onto H
# for specified energy of H. Data are taken from 

#   Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
# Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 305.
#
# Input: 
#   E        -     array or float, energy of H atom (target H atom is at rest)
#   VIS      -     if set, then return viscosity cross section instead of momentum transfer cross section 

# Gwendolyn Galleher

def Sigma_EL_H_H(E, vis = 0):
    if np.size(E) == 1: 
        _E = float(E) # this line was in the original code I don't 
                      # know if its neccessary but it means E must have only 1 element
    _E = np.maximum(_E, 3.03e0)     
    _E = np.minimum(_E, 1.01e4)
    if vis: 
        a = np.array([ -3.344860e1, -4.238982e-1, -7.477873e-2, -7.915053e-3, -2.686129e-4])
        result = np.exp(poly(np.log(_E), a)) * 1e-4
    else: 
        a = np.array([ -3.330843e1, -5.738374e-1, -1.028610e-1, -3.920980e-3, 5.964135e-4])
        result = np.exp(poly(np.log(_E), a)) * 1e-4
    if np.ndim(E) == 0:      
        result = result[0]
    return result


