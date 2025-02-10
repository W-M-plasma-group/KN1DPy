import numpy as np
import matplotlib.pyplot as plt

def sigma_cx_h0(E=None, coef='freeman'):
    """
    Compute the charge exchange cross-section for neutral hydrogen (H0).

    This function calculates the charge exchange cross-section for the reaction:
        H^+ + H(1s) → H(1s) + H^+
    as a function of the incident particle energy (E) using empirical coefficients.

    Parameters:
    -----------
    E : numpy.ndarray, optional
        Array of incident particle energies in eV. If None, a ValueError is raised.
    coef : str, optional
        Type of coefficient to use for the calculation. Options are:
        - 'freeman': Use coefficients from Freeman and Jones (1974).
        - 'janev': Use coefficients from Janev (1987).
        Default is 'freeman'.

    Returns:
    --------
    numpy.ndarray
        Array of charge exchange cross-sections in cm². The output has the same shape
        as the input array E.

    Raises:
    -------
    ValueError
        - If E is None.
        - If coef is not 'freeman' or 'janev'.

    Notes:
    ------
    - The input energy values (E) are clipped to ensure they fall within the valid range
      for the chosen coefficient set:
        - 'freeman': [0.1, 1.0e5] eV
        - 'janev': [0.1, 2.01e4] eV
    - The empirical coefficients are taken from:
        - Freeman, E. L., & Jones, E. M. (1974). "Atomic Collision Processes in Plasma Physics."
        - Janev, R. K., et al. (1987). "Elementary Processes in Hydrogen-Helium Plasmas:
          Cross Sections and Reaction Rate Coefficients." Springer-Verlag.
    """
    if E is None:
        # If E is not defined, raise an exception
        raise ValueError("The parameter E cannot be None.")
    
    # Convert E to a NumPy array with dtype=np.float64
    E = np.asarray(E, dtype=np.float64)  

    if coef == 'freeman':
        E = np.clip(E, np.float64(0.1), np.float64(1.0e5))
        ## Information: Janev 1987, p. 128
        return ((np.float64(1e-4) * np.float64(0.6937e-14)) * (np.float64(1.0) - np.float64(0.155) * np.log10(E))**2) / (np.float64(1.0) + np.float64(0.1112e-14) * (E**3.3))  
    elif coef == 'janev':
        E = np.clip(E, np.float64(0.1), np.float64(2.01e4))
        ## Information: Janev 1987, p. 250
        alpha = np.flip(np.array([-3.274123792568e+01, -8.916456579806e-02, -3.016990732025e-02,
                                   9.205482406462e-03,  2.400266568315e-03, -1.927122311323e-03,
                                   3.654750340106e-04, -2.788866460622e-05,  7.422296363524e-07], 
                                   dtype=np.float64))
        return np.exp(np.polyval(alpha, np.log(E))) * np.float64(1e4)
    else:
        raise ValueError("The coefficient must be 'freeman' or 'janev'.")

if __name__ == '__main__':
    E_values = np.logspace(-2, 6, num=100)
    # Compute the cross-section using Freeman's coefficients
    sigma_freeman = sigma_cx_h0(E_values, coef='freeman')
    # Compute the cross-section using Janev's coefficients
    sigma_janev   = sigma_cx_h0(E_values, coef='janev')
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8,6))
    plt.plot(E_values, sigma_freeman, label="Freeman")
    plt.plot(E_values, sigma_janev,   label="Janev")
    plt.title(r'Charge Exchange Cross Section for Atomic Hydrogen:$p + H(1s) \rightarrow H(1s) + p$')
    plt.axvline(0.1   ,color='k'    ,linestyle='--',alpha=0.4)
    plt.axvline(2.01e4,color='green',linestyle='--',alpha=0.7,label='Upper Limit Janev')
    plt.axvline(1.0e5 ,color='navy' ,linestyle='--',alpha=0.7,label='Upper Limit Freeman')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which='both')
    plt.xlabel("Energy (eV)")
    plt.ylabel(r"Cross-section ($m^{-2}$)")
    plt.legend(loc='best')
    plt.show()
