import numpy as np
from poly import poly

# Returns maxwellian averaged <sigma V) for electron impact
# ionization of molecular hydrogen. Coefficients are taken 
# from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
# Springer-Verlag, 1987, p.259.

def sigmav_ion_hh(Te):

#   Te - Scalar, list, or np.array. Gives electron temperature (eV)

  if type(Te) != np.ndarray: #converts scalar or list input to np.array
    Te=np.array(Te)
  b = [-3.568640293666e+1, 1.733468989961e+1, -7.767469363538e+0, 
         2.211579405415e+0, -4.169840174384e-1, 5.088289820867e-2, 
        -3.832737518325e-3, 1.612863120371e-4, -2.893391904431e-6]
  Te2=np.maximum(Te,0.1) #Sets values to minimum .1 and maximum 2.01e4
  Te2=np.minimum(Te2,2.01e4) 
  return np.e**(poly(np.log(Te2),b))*1e-6
