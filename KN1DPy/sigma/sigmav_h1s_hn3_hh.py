import numpy as np

from .poly import poly

def sigmav_h1s_hn3_hh(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron impact
    # dissociation of molecular hydrogen resulting in one H atom in
    # the 1s state and one H atom in the n=3 state. Coefficients are taken 
    # from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.259.

    # Also returns minimum, maximum, and average energy of the resultant H(1s), H(n=3) atoms.

    #the E0 variables are output keywords that contain energy values
    E0_ave = 2.5
    E0_max = 1.25
    E0_min = 3.75

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    b = [-3.884976142596e+1,  1.520368281111e+1, -6.078494762845e+0,  1.535455119900e+0, -2.628667482712e-1,  2.994456451213e-2, -2.156175515382e-3,  8.826547202670e-5, -1.558890013181e-6]
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)
    return np.e**(poly(np.log(Te), b))*1e-6