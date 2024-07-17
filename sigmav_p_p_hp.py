import numpy as np
from poly import poly

def sigmav_p_p_hp(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociation of molecular hydrogen ions resulting in 
    # two protons. Coefficients are taken from Janev, 
    # "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.260.

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-3.746192301092e+1, 
         1.559355031108e+1, 
        -6.693238367093e+0, 
         1.981700292134e+0, 
        -4.044820889297e-1, 
         5.352391623039e-2, 
        -4.317451841436e-3, 
         1.918499873454e-4, 
        -3.591779705419e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6