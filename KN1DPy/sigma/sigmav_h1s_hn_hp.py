import numpy as np

from .poly import poly

def sigmav_h1s_hn_hp(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociative recombination of molecular hydrogen ions resulting in 
    # one H atom in the 1s state and one H atom in state n > or = 2. Coefficients 
    # are taken from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.260.

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-1.670435653561e+1, 
        -6.035644995682e-1, 
        -1.942745783445e-8, 
        -2.005952284492e-7, 
         2.962996104431e-8, 
         2.134293274971e-8, 
        -6.353973401838e-9, 
         6.152557460831e-10, 
        -2.025361858319e-11]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6