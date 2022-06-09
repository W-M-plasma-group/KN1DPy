import numpy as np
from poly import poly

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
                
  is_scalar=False
  if type(E) not in [np.ndarray, list]: #converts scalar input to iterable list
    E=[E]
    is_scalar=True
  if freeman:
    E2=[e if e>.1 else .1 for e in E] #Sets values to minimum .1 and maximum 1e5
    E2=np.array([e if e<1e5 else 1e5 for e in E2]) #And converts to np.array
    result=1.0e-4 * 0.6937e-14*(1.0-0.155*np.log10(E2))**2/(1.0+0.1112e-14*E2**3.3)
  else:
    E2=[e if e>.1 else .1 for e in E] #Sets values to minimum .1 and maximum 2.01e4
    E2=np.array([e if e<2.01e4 else 2.01e4 for e in E2]) #And converts to np.array
    alpha=[-3.274123792568e+01, -8.916456579806e-02, -3.016990732025e-02,
              9.205482406462e-03,  2.400266568315e-03, -1.927122311323e-03,
              3.654750340106e-04, -2.788866460622e-05,  7.422296363524e-07]
    result=np.e**(poly(np.log(E2),alpha))*1e-4
  if(is_scalar): #if the input was a scalar, converts result back to a scalar as well
    result=result[0]
  return result
