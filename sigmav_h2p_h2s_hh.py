import numpy as np
from poly import poly

def sigmav_h2p_h2s_hh(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociation of molecular hydrogen resulting in one H atom in
    # the 2p state and one H atom in the 2s state. Coefficients are taken 
    # from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.259.

    # Also returns minimum, maximum, and average energy of the resultant H(2p), H(2s) atoms.

    #the E0 variables are output keywords that contain energy values
    E0_ave = 4.85
    E0_max = 5.85
    E0_min = 2.85

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-4.794288960529e+1, 
         2.629649351119e+1, 
        -1.151117702256e+1, 
         2.991954880790e+0, 
        -4.949305181578e-1, 
         5.236320848415e-2, 
        -3.433774290547e-3, 
         1.272097387363e-4, 
        -2.036079507592e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6