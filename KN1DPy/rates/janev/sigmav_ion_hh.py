import numpy as np

from ...utils import poly

def sigmav_ion_hh(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact
    ionization of molecular hydrogen. Coefficients are taken from 
    
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

    b = [-3.568640293666e+1,
          1.733468989961e+1, 
         -7.767469363538e+0, 
          2.211579405415e+0, 
         -4.169840174384e-1, 
          5.088289820867e-2, 
         -3.832737518325e-3, 
          1.612863120371e-4, 
         -2.893391904431e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result
