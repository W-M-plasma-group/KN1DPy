import numpy as np
from poly import poly

def sigmav_p_h1s_hp(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociation of molecular hydrogen ions resulting in 
    # one proton and one H atom in the 1s state. Coefficients are taken 
    # from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.260.

    # Also returns minimum, maximum, and average energy of the resultant proton and H(1s) atom.

    #the E0 variables are output keywords that contain energy values
    E0_ave = 4.3
    E0_max = 4.3
    E0_min = 4.3

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-1.781416067709e+1, 
         2.277799785711e+0, 
        -1.266868411626e+0, 
         4.296170447419e-1, 
        -9.609908013189e-2, 
         1.387958040699e-2, 
        -1.231349039470e-3, 
         6.042383126281e-5, 
        -1.247521040900e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6