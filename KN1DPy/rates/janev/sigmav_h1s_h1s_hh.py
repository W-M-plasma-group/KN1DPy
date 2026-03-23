import numpy as np

from ...utils import poly

def sigmav_h1s_h1s_hh(Te):
    '''
    Computes maxwellian averaged <sigma V) for electron impact
    dissociation of molecular hydrogen resulting in two H atoms in
    the 1s state. Coefficients are taken from 
    
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

    b = [-2.787217511174e+1, 
          1.052252660075e+1, 
         -4.973212347860e+0,
          1.451198183114e+0,
         -3.062790554644e-1, 
          4.433379509258e-2,
         -4.096344172875e-3, 
          2.159670289222e-4, 
         -4.928545325189e-6]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)

    result = np.exp(poly(np.log(Te), b))*1e-6
    return result 