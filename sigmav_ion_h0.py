import numpy as np
from poly import poly

# Returns maxwellian averaged <sigma V) for electron impact
# ionization of atomic hydrogen. Coefficients are taken 
# from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
# Springer-Verlag, 1987, p.258.

def sigmav_ion_h0(Te):

#   Te - Scalar, list, or np.array. Gives electron temperature (eV)

  if type(Te) != np.ndarray: #converts scalar or list input to np.array
    Te=np.array(Te)
  b = [-3.271396786375e+1, 1.353655609057e+1, -5.739328757388e+0, 
         1.563154982022e+0, -2.877056004391e-1, 3.482559773737e-2, 
        -2.631976175590e-3, 1.119543953861e-4, -2.039149852002e-6]
  Te2=np.maximum(Te,0.1) #Sets values to minimum .1 and maximum 2.01e4
  Te2=np.minimum(Te2,2.01e4) 
  return np.e**(poly(np.log(Te2),b))*1e-6