import numpy as np
from poly import poly

def sigmav_p_hn2_hp(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociation of molecular hydrogen ions resulting in 
    # one proton and one H atom in the n=2 state. Coefficients are taken 
    # from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.260.

    # Also returns minimum, maximum, and average energy of the resultant proton and H(n=2) atom.

    #the E0 variables are output keywords that contain energy values
    E0_ave = 1.5
    E0_max = 1.5
    E0_min = 1.5

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-3.408905929046e+1, 
         1.573560727511e+1, 
        -6.992177456733e+0, 
         1.852216261706e+0, 
        -3.130312806531e-1, 
         3.383704123189e-2, 
        -2.265770525273e-3, 
         8.565603779673e-5, 
        -1.398131377085e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6