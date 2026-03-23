import numpy as np

from ...utils import poly

def sigmav_p_h1s_hp(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact dissociation 
    of molecular hydrogen ions resulting in  one proton and one H atom 
    in the 1s state. Coefficients are taken from 
    
        Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
        Springer-Verlag, 1987, p.260.

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
    
    b = [-1.781416067709e+1, 
          2.277799785711e+0, 
         -1.266868411626e+0, 
          4.296170447419e-1, 
         -9.609908013189e-2, 
          1.387958040699e-2, 
         -1.231349039470e-3, 
          6.042383126281e-5, 
         -1.247521040900e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result