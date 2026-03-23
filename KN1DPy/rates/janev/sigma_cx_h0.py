import numpy as np

from ...utils import poly

def sigma_cx_h0(E, freeman=False):
    '''
    Computes charge exchange cross section for atomic hydrogen. Data are taken either from the polynomial fit in
        Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.250.
    or from Freeman and Jone's analytic fit tabulated in 
        Freeman, E.L., Jones, E.M., "Atomic Collision Processes in Plasma Physics
        Experiments", UKAEA Report No. CLM-R137 (Culham Laboratory, Abington, England 1974)

    Parameters
    ----------
        E : ndarray or float
            energy of proton corresponding to the relative velocity between proton and hydrogen atom (eV)
        freeman: bool, default=false
            if true, then return CX based on Freeman and Jones' analytic fit in
            Freeman, E.L., Jones, E.M., "Atomic Collision Processes in Plasma Physics
            Experiments", UKAEA Report No. CLM-R137 (Culham Laboratory, Abington, England 1974).
            Otherwise, return CX based on polynomial fit in
		    Janev, "Elementary Processes in Hydrogen-Helium Plasmas", 
	        Springer-Verlag, 1987, p.250, other

    Returns
    -------
        ndarray
            sigma_CX for 0.1 < E < 2e4 (m^-2)
    '''
                
    E = np.asarray(E)
    if freeman:
        #Set values to minimum .1 and maximum 1e5
        E2 = np.clip(E, 0.1, 1e5)
        result = 1.0e-4 * 0.6937e-14*(1.0 - 0.155*np.log10(E2))**2/(1.0 + 0.1112e-14*E2**3.3)
    else:
        #Sets values to minimum .1 and maximum 2.01e4
        E2 = np.clip(E, 0.1, 2.01e4)
        alpha = [-3.274123792568e+01, -8.916456579806e-02, -3.016990732025e-02,
                  9.205482406462e-03,  2.400266568315e-03, -1.927122311323e-03,
                  3.654750340106e-04, -2.788866460622e-05,  7.422296363524e-07]
        result = np.exp(poly(np.log(E2), alpha))*1e-4

    return result
