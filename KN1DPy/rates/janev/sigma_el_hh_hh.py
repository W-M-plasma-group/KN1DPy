import numpy as np

from ...utils import poly

def sigma_el_hh_hh(E, vis = False):
    '''
    Computes momentum transfer cross section for elastic collisions of H2 onto H2 
    for specified energy of H2. Data are taken from 

        Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
        Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
        Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
        Passovets, page 305.

    Parameters
    ----------
    E : ndarray or float
        energy of H2 molecule (target H2 molecule is at rest)
    vis : bool, defaul=False
        if true, then return viscosity cross section instead of momentum transfer cross section 
    
    Returns
    -------
        ndarray
            Sigma for 0.03 < E < 1e4. For E outside this range, 
            the value of Sigma at the 0.03 or 1e4 eV boundary is returned. (m^-2)
    '''
    E = np.asarray(E, dtype=float)

    # Ensure 0.03e0 < E < 1.01e4
    E = np.clip(E, 0.03e0, 1.01e4)

    if vis: #NOTE I am confused, why are these seperate statements, and why is this a warning?
        print('WARNING in SIGMA_EL_HH_HH => using momentum transfer as viscosity cross-section')
        # calculates viscosity cross section
        a = np.array([-3.430345e1, -2.960406e-1, -6.382532e-2, -7.557519e-3, 2.606259e-4])
    else: 
        # calculates momentum transfer cross section 
        a = np.array([-3.430345e1, -2.960406e-1, -6.382532e-2, -7.557519e-3, 2.606259e-4])
    
    result = np.exp(poly(np.log(E), a)) * 1e-4 
    return result