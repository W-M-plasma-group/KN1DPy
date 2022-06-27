import numpy as np
from poly import poly

def sigmav_h1s_hn3_hh(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # ionization and disociation of Molecular hydrogen resulting in one 
    # proton and one H atom in the 1s state. Coefficients are taken 
    # from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.260.

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-3.834597006782e+1, 
         1.426322356722e+1, 
        -5.826468569506e+0, 
         1.727940947913e+0, 
        -3.598120866343e-1, 
         4.822199350494e-2, 
        -3.909402993006e-3, 
         1.738776657690e-4, 
        -3.252844486351e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6