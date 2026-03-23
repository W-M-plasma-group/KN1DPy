import numpy as np

from ...utils import poly

def sigmav_h1s_hn_hp(Te):
    '''
    Computes maxwellian averaged <sigma V) for electron impact dissociative 
    recombination of molecular hydrogen ions resulting in one H atom in the 
    1s state and one H atom in state n > or = 2. Coefficients are taken from 
    
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
    
    b = [-1.670435653561e+1, 
         -6.035644995682e-1, 
         -1.942745783445e-8, 
         -2.005952284492e-7, 
          2.962996104431e-8, 
          2.134293274971e-8, 
         -6.353973401838e-9, 
          6.152557460831e-10, 
         -2.025361858319e-11]
    
    # Ensure 0.1 < Te < 2.01e4
    Te = np.clip(Te, 0.1, 2.01e4)
    
    result = np.exp(poly(np.log(Te), b))*1e-6
    return result