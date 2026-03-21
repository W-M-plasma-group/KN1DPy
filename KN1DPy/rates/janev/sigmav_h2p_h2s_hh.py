import numpy as np

from ...utils import poly

def sigmav_h2p_h2s_hh(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact dissociation
    of molecular hydrogen resulting in one H atom in the 2p state and one
    H atom in the 2s state. Coefficients are taken from 

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
    
    b = [-4.794288960529e+1, 
          2.629649351119e+1, 
         -1.151117702256e+1, 
          2.991954880790e+0, 
         -4.949305181578e-1, 
          5.236320848415e-2, 
         -3.433774290547e-3, 
          1.272097387363e-4, 
         -2.036079507592e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result