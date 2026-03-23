import numpy as np

from ...utils import poly

def sigmav_h1s_hn3_hh(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact dissociation
    of molecular hydrogen resulting in one H atom in the 1s state and one 
    H atom in the n=3 state. Coefficients are taken from 
    
        Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
        Springer-Verlag, 1987, p.259.

    Parameters
    ----------
    Te : ndarray or float
        electron temperature (eV)

    Returns
    -------
        ndarray
            Sigma V for 0.1 < Te < 2e4. (m^3/s)
    '''

    Te = np.asarray(Te)
    
    b = [-3.884976142596e+1,  
          1.520368281111e+1, 
         -6.078494762845e+0,
          1.535455119900e+0,
         -2.628667482712e-1,
          2.994456451213e-2, 
         -2.156175515382e-3,
          8.826547202670e-5, 
         -1.558890013181e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result