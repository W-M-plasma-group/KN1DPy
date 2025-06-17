import numpy as np
from poly import poly 
#   Returns momentum transfer cross section for elastic collisions of H+ onto H2
# for specified energy of H+.
# Data are taken from 
# Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
# Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
# Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
# Passovets, page 305.
# Inputs: 
#   E   -   array or float, energy of H2 molecule (target H2 molecule is at rest)
# Output:
#	Returns Sigma for 0.03 < E < 1e4. For E outside this range, 
#   the value of Sigma at the 0.03 or 1e4 eV boundary is returned.
#
#	units: m^-2
# Gwendolyn Galleher 
def sigma_el_p_hh(E):
    E = np.array([E])
    _E = E.astype(float)
    # ensures 0.03e0 < E < 1.01e4
    _E = np.maximum(_E, 0.03e0)
    _E = np.minimum(_E, 1.01e4)
    a = np.array([-3.355719e1, -5.696568e-1, -4.089556e-2, -1.143513e-2, 5.926596e-4])
    result = np.exp(poly(np.log(_E), a)) * 1e-4
    return result 