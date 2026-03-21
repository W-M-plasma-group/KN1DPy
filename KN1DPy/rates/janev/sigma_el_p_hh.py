import numpy as np

from ...utils import poly

def sigma_el_p_hh(E):
    '''
    Computes momentum transfer cross section for elastic collisions of H+ onto H2 
    for specified energy of H+. Data are taken from 

        Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
        Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
        Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
        Passovets, page 305.

    Parameters
    ----------
    E : ndarray or float
        energy of H+ ion (target H2 molecule is at rest)

    Returns
    -------
        ndarray
            Sigma for 0.03 < E < 1e4. For E outside this range, 
            the value of Sigma at the 0.03 or 1e4 eV boundary is returned. (m^-2)
    '''

    E = np.asarray(E, dtype=float)

    # Ensure 0.03e0 < E < 1.01e4
    E = np.clip(E, 0.03e0, 1.01e4)
    a = np.array([-3.355719e1, 
                  -5.696568e-1, 
                  -4.089556e-2, 
                  -1.143513e-2, 
                   5.926596e-4])
    
    result = np.exp(poly(np.log(E), a)) * 1e-4
    return result 