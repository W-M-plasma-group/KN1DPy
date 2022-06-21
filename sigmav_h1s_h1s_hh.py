import numpy as np
from poly import poly

def sigmav_h1s_h1s_hh(Te):
  #   Te - Scalar, list, or np.array. Gives electron temperature (eV)

    #   Returns maxwellian averaged <sigma V) for electron impact
    #   dissociation of molecular hydrogen resulting in two H atoms in
    #   the 1s state. Coefficients are taken from Janev, 
    #   "Elementary Processes in Hydrogen-Helium Plasmas",
    #   Springer-Verlag, 1987, p.259.

  if type(Te) != np.ndarray: #converts scalar or list input to np.array
    Te=np.array(Te)

  #   E0_ave, E0_max, and E0_min are output keywords

  E0_ave=3
  E0_max=4.25
  E0_min=2

  b = [-2.787217511174e+1, 1.052252660075e+1, -4.973212347860e+0,
         1.451198183114e+0, -3.062790554644e-1, 4.433379509258e-2,
        -4.096344172875e-3, 2.159670289222e-4, -4.928545325189e-6]
  Te2=np.maximum(Te,0.1) #Sets values to minimum .1 and maximum 2.01e4
  Te2=np.minimum(Te2,2.01e4) 
  return np.e**(poly(np.log(Te2),b))*1e-6