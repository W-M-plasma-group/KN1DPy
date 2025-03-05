import numpy as np

def sigma_v_p_h1s_hh(Te):
    """
    Returns Maxwellian-averaged <sigma V> for electron impact ionization and dissociation of molecular hydrogen,
    resulting in one proton and one H atom in the 1s state. Coefficients are taken from Janev,
    "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.260.

    Parameters:
        Te (float or array-like): Electron temperature in eV. Must be within the range [0.1, 2e4] eV.

    Returns:
        float or numpy.ndarray: Maxwellian-averaged <sigma V> in m^3/s.
                                For temperatures outside the valid range, the boundary values are returned.
    """
    # Ensure Te is a NumPy array and clip to valid temperature range
    Te = np.asarray(Te, dtype=np.float64)
    Te = np.clip(Te, 0.1, 2.01e4)

    # Coefficients for the polynomial expansion
    b = np.array([
        -3.834597006782e+1,         1.426322356722e+1,        -5.826468569506e+0,
         1.727940947913e+0,        -3.598120866343e-1,         4.822199350494e-2,
        -3.909402993006e-3,         1.738776657690e-4,        -3.252844486351e-6
    ])

    # Compute the Maxwellian-averaged <sigma V> using the polynomial expansion
    sigma_v = np.exp(np.polyval(b[::-1], np.log(Te))) * 1e-6

    return sigma_v

# Example usage
if __name__ == '__main__':
    # Test with a scalar value
    Te_scalar = 100  # Example electron temperature in eV
    print(f"<sigma V> for Te={Te_scalar} eV: {sigma_v_p_h1s_hh(Te_scalar)} m^3/s")

    # Test with an array of temperatures
    Te_array = np.logspace(-2, 5, 1000)  # Temperatures from 0.1 eV to 1e4 eV
    result = sigma_v_p_h1s_hh(Te_array)

    # Plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(Te_array, result, label=r'$e+H_2\rightarrow e+H(1s)+e$')
    plt.axvline(0.1, color='k', linestyle='--', alpha=0.4, label='Temperature Bounds')
    plt.axvline(2.01e4, color='k', linestyle='--', alpha=0.4)
    plt.title(r'Cross-Section for $H_2$ Ionization/Dissociation:$e+H_2\rightarrow e+H(1s)+e$', fontsize=14)
    plt.xlabel('Electron Temperature (eV)', fontsize=12)
    plt.ylabel(r'$<\sigma V>$ ($m^3/s$)', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()