import numpy as np

def sigmav_rec_h1s(Te):
    #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #Returns maxwellian averaged <sigma V) for electron-ion radiative
    # recombination to the atomic hydrogen in the 1s state.
    # Coefficients are taken from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    # Springer-Verlag, 1987, p.32.

    if type(Te) != np.ndarray:
        Te = np.array(Te) #if Te isn't an array it makes it one
    
    Te = np.maximum(Te, 0.1) #makes sure 0.1 < Te < 2.01e4
    Te = np.minimum(Te, 2.01e4)

    #data for nl = 1s
    n = 1
    Ry = 13.58
    Eion_n = Ry/n
    Anl = 3.92
    Xnl = 0.35

    Bn = Eion_n/Te
    return Anl*1e-14*np.sqrt(Eion_n/Ry)*(Bn**1.5/(Bn+Xnl))*1e-6