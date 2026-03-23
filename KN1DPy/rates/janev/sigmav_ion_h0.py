import numpy as np

from ...utils import poly

def sigmav_ion_h0(Te):
    '''
    Returns maxwellian averaged <sigma V) for electron impact
    ionization of atomic hydrogen. Coefficients are taken from 
    
        Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
        Springer-Verlag, 1987, p.258.

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

    b = [-3.271396786375e+1, 
          1.353655609057e+1, 
         -5.739328757388e+0, 
          1.563154982022e+0, 
         -2.877056004391e-1,
          3.482559773737e-2, 
         -2.631976175590e-3, 
          1.119543953861e-4, 
         -2.039149852002e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result