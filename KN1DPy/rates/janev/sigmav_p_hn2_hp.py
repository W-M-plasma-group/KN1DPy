import numpy as np

from ...utils import poly

def sigmav_p_hn2_hp(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact dissociation 
    of molecular hydrogen ions resulting in one proton and one H atom 
    in the n=2 state. Coefficients are taken from 
    
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
    
    b = [-3.408905929046e+1, 
          1.573560727511e+1, 
         -6.992177456733e+0, 
          1.852216261706e+0, 
         -3.130312806531e-1, 
          3.383704123189e-2, 
         -2.265770525273e-3, 
          8.565603779673e-5, 
         -1.398131377085e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)
    
    return np.exp(poly(np.log(Te), b))*1e-6