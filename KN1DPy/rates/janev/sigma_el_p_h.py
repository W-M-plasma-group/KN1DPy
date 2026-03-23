import numpy as np 

from ...utils import poly 

def sigma_el_p_h(E):
    '''
    Computes momentum transfer cross section for elastic collisions of H+ onto H 
    for specified energy of H+. Data are taken from 

        Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
        Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
        Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
        Passovets, page 298.

    Parameters
    ----------
    E : ndarray or float
        energy of H+ ion (target H atom is at rest

    Returns
    -------
        ndarray
            Sigma for 0.001 < E < 1e5. For E outside this range, 
            the value of Sigma at the 0.001 or 1e5 eV boundary is returned. (m^-2)
    '''

    E = np.asarray(E, dtype=float)
    
    # Ensure that 0.001e0 < E < 1.01e5
    E = np.clip(E, 0.001, 1.01e5)
    
    result = np.zeros_like(E)
    logE = np.log(E)

    low = E < 10.0
    high = ~low

    if np.any(low):
        a_low = np.array([
            -3.233966e1,
            -1.126918e-1,
             5.287706e-3,
            -2.445017e-3,
            -1.044156e-3,
             8.419691e-5,
             3.824773e-5
        ])
        result[low] = np.exp(poly(logE[low], a_low)) * 1e-4

    if np.any(high):
        a_high = np.array([
            -3.231141e1,
            -1.386002e-1
        ])
        result[high] = np.exp(poly(logE[high], a_high)) * 1e-4

    return result
