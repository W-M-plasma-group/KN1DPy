import numpy as np
from poly import poly

def sigmav_h1s_h2s_hh(Te):
  #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #   Returns maxwellian averaged <sigma V) for electron impact
    #   dissociation of molecular hydrogen resulting in one H atom in
    #   the 1s state and one H atom in the 2s state. Coefficients are taken 
    #   from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
    #   Springer-Verlag, 1987, p.259.

  if type(Te) != np.ndarray: #converts scalar or list input to np.array
    Te=np.array(Te)

  #   E0_ave, E0_max, and E0_min are output keywords

  E0_ave=3
  E0_max=.55
  E0_min=0

  b = [-3.454175591367e+1, 1.412655911280e+1, -6.004466156761e+0,
        1.589476697488e+0, -2.775796909649e-1, 3.152736888124e-2,
        -2.229578042005e-3, 8.890114963166e-5, -1.523912962346e-6]
  Te2=np.maximum(Te,0.1) #Sets values to minimum .1 and maximum 2.01e4
  Te2=np.minimum(Te2,2.01e4) 
  result = np.e**(poly(np.log(Te2),b))*1e-6  # - added the result variable to make it easier to read - GG
  return result 
