import numpy as np

def sigma_v_rec_h1s(Te):
    """
    Calculate the electron-ion radiative recombination rate coefficient
    for atomic hydrogen in the 1s state.

    This function computes the rate coefficient for the reaction:
        e + H^+ → H(n=1) + hν
    based on the empirical formula provided in Janev (1987), p. 32.

    Parameters:
    -----------
    Te : float or numpy.ndarray
        Electron temperature in eV. Can be a scalar or an array of temperatures.

    Returns:
    --------
    result : float or numpy.ndarray
        Recombination rate coefficient in cm³/s. The output type matches the input type:
        - If Te is a scalar, result is a scalar.
        - If Te is an array, result is an array of the same shape.

    Notes:
    ------
    - The temperature Te is clipped to the range [0.1, 2.01e4] eV to ensure
      the validity of the empirical formula.
    - The formula is based on the following parameters:
        - n = 1 (principal quantum number of the final state)
        - Ry = 13.58 eV (Rydberg energy for hydrogen)
        - A_{nl=1s} = 3.92 (empirical coefficient)
        - X_{nl=1s} = 0.35 (empirical coefficient)
    """
    # Clip Te to the valid range [0.1, 2.01e4] eV
    Te = np.clip(Te, 0.1, 2.01e4)

    # Constants for the recombination rate formula
    n       = 1  # Principal quantum number of the final state
    Ry      = 13.58  # Rydberg energy for hydrogen (eV)
    E_ion_n = Ry / n  # Ionization energy of the n-th level (eV)
    A_nl    = 3.92  # Empirical coefficient
    X_nl    = 0.35  # Empirical coefficient

    # Calculate the recombination rate coefficient based on the empirical formula provided in Janev (1987), p. 32.
    B_n = E_ion_n / Te
    result = (A_nl * (1e-14)) * np.sqrt(E_ion_n / Ry) * ((B_n**1.5) / (B_n + X_nl)) * 1e-6

    return result

if __name__ == '__main__':
    print('Te:scalar')
    Te = 1e4  # scalar
    result = sigma_v_rec_h1s(Te)
    print(f'Te:{Te},result:{result}')

    print('Te:array')
    Te = np.logspace(-2,6,1000)
    result = sigma_v_rec_h1s(Te)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    plt.plot(Te, result, label=r'$\text{e} + \text{H}^+ \rightarrow \text{H}(1s) + h\nu$')
    plt.yscale('log')
    plt.xscale('log')
    plt.axvline(0.1   ,color='k',linestyle='--',alpha=0.4)
    plt.axvline(2.01e4,color='k',linestyle='--',alpha=0.4)
    plt.grid(True, which='both', linestyle='-', linewidth=0.5)
    plt.xlabel('Te (eV)')
    plt.ylabel(r'$\langle \sigma v \rangle$ (m$^3$/s)')
    plt.legend()
    plt.title(r'$\text{e} + \text{H}^+ \rightarrow \text{H}(1s) + h\nu$')
    plt.show()