import numpy as np

from .poly import poly

#Returns charge exchange cross-section for atomic hydrogen.

def sigma_cx_h0(E, freeman=0):

#	E - Scalar value, list, or np.array 
#       energy of proton corresponding to the relative velocity between proton and hydrogen atom. (eV)

#	Freeman - if set, then return CX based on Freeman and Jones' analytic fit in
#		Freeman, E.L., Jones, E.M., "Atomic Collision Processes in Plasma Physics
#       Experiments", UKAEA Report No. CLM-R137 (Culham Laboratory, Abington, England 1974)
#
#       Otherwise, return CX based on polynomial fit in
#		Janev, "Elementary Processes in Hydrogen-Helium Plasmas", 
#	    Springer-Verlag, 1987, p.250, other
                
  if type(E) != np.ndarray: #converts scalar or list input to np.array
    E=np.array(E)
  if np.any(freeman): # edit so that logic statement does not give error - GG
    E2=np.maximum(E,0.1) #Sets values to minimum .1 and maximum 1e5
    E2=np.minimum(E2,1e5) 
    return 1.0e-4 * 0.6937e-14*(1.0 - 0.155*np.log10(E2))**2/(1.0 + 0.1112e-14*E2**3.3)
  else:
    E2=np.maximum(E,0.1) #Sets values to minimum .1 and maximum 2.01e4
    E2=np.minimum(E2,2.01e4) 
    alpha=[-3.274123792568e+01, -8.916456579806e-02, -3.016990732025e-02,
              9.205482406462e-03,  2.400266568315e-03, -1.927122311323e-03,
              3.654750340106e-04, -2.788866460622e-05,  7.422296363524e-07]
    return np.e**(poly(np.log(E2),alpha))*1e-4
