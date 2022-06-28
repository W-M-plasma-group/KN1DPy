import numpy as np
from poly import poly

def sigma_cx_hh(E):
    #E	- fltarr(*) or float, energy of molecule corresponding to the
    #relative velocity between molecule and molecular ion. (eV)

    # Returns charge exchange cross section for molecular hydrogen. Data are taken
    # the polynomial fit in

    # Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.253.
    
    if type(E) != np.ndarray:
        T = np.array(E) #if E isn't an array it makes it one

    E = np.maximum(E, 0.1) #makes sure 0.1 < E < 2.01e4
    E = np.minimum(E, 2.01e4)
    alpha = [-3.427958758517e+01, -7.121484125189e-02, 4.690466187943e-02,
          -8.033946660540e-03, -2.265090924593e-03,-2.102414848737e-04,
           1.948869487515e-04, -2.208124950005e-05, 7.262446915488e-07]

    return np.e**(poly(np.log(E), alpha))*1e-4