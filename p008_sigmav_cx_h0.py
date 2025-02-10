import numpy as np
import matplotlib.pyplot as plt

def sigma_v_cx_h0(E:np.ndarray,T:np.ndarray):
    """
    Calculate the charge exchange reaction rate coefficient for hydrogen.

    This function computes the reaction rate coefficient for the charge exchange process:
        H^+ + H(1s) → H(1s) + H^+
    as a function of the incident particle energy (E) and the plasma temperature (T).
    The calculation is based on empirical coefficients from Janev (1987), p. 272.

    Parameters:
    -----------
    E : numpy.ndarray
        Array of incident particle energies in eV. Must have the same shape as T.
    T : numpy.ndarray
        Array of plasma temperatures in eV. Must have the same shape as E.

    Returns:
    --------
    result : numpy.ndarray
        Array of reaction rate coefficients in m³/s. The output has the same shape
        as the input arrays E and T.

    Raises:
    -------
    ValueError
        If E and T do not have the same shape.

    Notes:
    ------
    - The input values of E and T are clipped to the range [0.1, 2.01e4] eV to ensure
      the validity of the empirical formula.
    - The empirical coefficients are taken from:
        Janev, R. K., et al. (1987). "Elementary Processes in Hydrogen-Helium Plasmas:
        Cross Sections and Reaction Rate Coefficients." Springer-Verlag.
    - The formula used is a polynomial expansion in log(E) and log(T), with coefficients
      stored in the matrix `alpha`.
    """

    # Ensure E and T have the same shape
    if E.shape != T.shape:
        raise ValueError("E & T must have the same length")

    # Clip E and T to the valid range [0.1, 2.01e4] eV
    E_clip = np.clip(E, 0.1, 2.01e4)
    T_clip = np.clip(T, 0.1, 2.01e4)

    # Empirical coefficients for the polynomial expansion
    alpha = np.zeros((9,9))
    alpha[:, 0:3] = np.array([
        [-1.829079581680E+01,	  1.640252721210E-01,	  3.364564509137E-02],
        [ 2.169137615703E-01,    -1.106722014459E-01,    -1.382158680424E-03],
        [ 4.307131243894E-02,	  8.948693624917E-03,    -1.209480567154E-02],
        [-5.754895093075E-04,     6.062141761233E-03,     1.075907881928E-03],
        [-1.552077120204E-03,    -1.210431587568E-03,     8.297212635856E-04],
        [-1.876800283030E-04,    -4.052878751584E-05,    -1.907025662962E-04],
        [ 1.125490270962E-04,     2.875900435985E-05,     1.338839628570E-05],
        [-1.238982763007E-05,    -2.616998139678E-06,    -1.171762874107E-07],
        [ 4.163596197181E-07,     7.558092849125E-08,    -1.328404104165E-08]
    ])
    alpha[:, 3:6] = np.array([
        [ 9.530225559189E-03,    -8.519413589968E-04,    -1.247583860943E-03],
        [ 7.348786286628E-03,    -6.343059502294E-04,    -1.919569450380E-04],
        [-3.675019470470E-04,     1.039643390686E-03,    -1.553840717902E-04],
        [-8.119301728339E-04,     8.911036876068E-06,     3.175388949811E-05],
        [ 1.361661816974E-04,    -1.008928628425E-04,     1.080693990468E-05],
        [ 1.141663041636E-05,     1.775681984457E-05,    -3.149286923815E-06],
        [-4.340802793033E-06,    -7.003521917385E-07,     2.318308730487E-07],
        [ 3.517971869029E-07,    -4.928692832866E-08,     1.756388998863E-10],
        [-9.170850253981E-09,     3.208853883734E-09,    -3.952740758950E-10]
    ])
    alpha[:, 6:9] = np.array([
        [ 3.014307545716E-04,    -2.499323170044E-05,     6.932627237765E-07],
        [ 4.075019351738E-05,    -2.850044983009E-06,     6.966822400446E-08],
        [ 2.670827249272E-06,     7.695300597935E-07,    -3.783302281524E-08],
        [-4.515123641755E-06,     2.187439283954E-07,    -2.911233951880E-09],
        [ 5.106059413591E-07,    -1.299275586093E-07,     5.117133050290E-09],
        [ 3.105491554749E-08,     2.274394089017E-08,    -1.130988250912E-09],
        [-6.030983538280E-09,    -1.755944926274E-09,     1.005189187279E-10],
        [-1.446756795654E-10,     7.143183138281E-11,    -3.989884105603E-12],
        [ 2.739558475782E-11,    -1.693040208927E-12,     6.388219930167E-14]
    ])
    
    # Compute log(E) and log(T)
    logE = np.log(E_clip)
    logT = np.log(T_clip)
    
    # Initialize the result array
    result = np.zeros_like(logE)

    # Compute the polynomial expansion
    for n in range(9):
        for m in range(9):
            ## Information: Janev 1987, p. 272
            result = result + alpha[n, m] * (logE**n) * (logT**m)
    
    # Convert the result to linear scale and adjust units
    result = np.exp(result) * 1E-6
    return result

if __name__ == '__main__':
    # Define E & T
    E = np.logspace(-2, 5, 100)
    T = np.logspace(-2, 5, 100)
    result = sigma_v_cx_h0(T,E)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(E, result,label=r'$\langle \sigma v \rangle$  vs E')
    plt.plot(T, result,label=r'$\langle \sigma v \rangle$  vs T')
    plt.yscale('log')
    plt.xscale('log')
    plt.axvline(0.1   ,color='k',linestyle='--',alpha=0.4)
    plt.axvline(2.01e4,color='k',linestyle='--',alpha=0.4)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel(r'E or T ($eV$)')
    plt.ylabel(r'$\langle \sigma v \rangle$ ($m^{3}/s$)')
    plt.legend()
    plt.title(r'$\langle \sigma v \rangle$  vs E, T')
    plt.show()