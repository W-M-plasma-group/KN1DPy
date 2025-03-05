import numpy as np

def sigma_v_p_hn2_hp(Te):
    """
    Returns Maxwellian-averaged <sigma V> for electron impact dissociation of molecular hydrogen ions
    resulting in one proton and one H atom in the n=2 state. Coefficients are taken from Janev,
    "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.260.

    Also returns minimum, maximum, and average energy of the resultant proton and H(n=2) atom.

    Parameters:
        Te (float or array-like): Electron temperature in eV. Must be within the range [0.1, 2e4] eV.

    Returns:
        dict: A dictionary containing:
            - 'sigma_v': Maxwellian-averaged <sigma V> in m^3/s.
            - 'E0_ave': Average energy of P and H(n=2) atom (eV).
            - 'E0_min': Minimum energy of P and H(n=2) atom (eV).
            - 'E0_max': Maximum energy of P and H(n=2) atom (eV).
    """
    # Ensure Te is a NumPy array and clip to valid temperature range
    Te = np.asarray(Te, dtype=np.float64)
    Te = np.clip(Te, 0.1, 2.01e4)

    # Coefficients for the polynomial expansion
    b = np.array([
        -3.408905929046e+1,
         1.573560727511e+1,
        -6.992177456733e+0,
         1.852216261706e+0,
        -3.130312806531e-1,
         3.383704123189e-2,
        -2.265770525273e-3,
         8.565603779673e-5,
        -1.398131377085e-6
    ])

    # Compute the Maxwellian-averaged <sigma V> using the polynomial expansion
    sigma_v = np.exp(np.polyval(b[::-1], np.log(Te))) * 1e-6

    # Energy values for the resultant particles (constant as per the IDL code)
    E0_ave = 1.5  # Average energy (eV)
    E0_min = 1.5  # Minimum energy (eV)
    E0_max = 1.5  # Maximum energy (eV)

    # Return results as a dictionary
    return {
        'sigma_v': sigma_v,
        'E0_ave': E0_ave,
        'E0_min': E0_min,
        'E0_max': E0_max
    }

# Example usage
if __name__ == '__main__':
    # Test with a scalar value
    Te_scalar = 100  # Example electron temperature in eV
    result_scalar = sigma_v_p_hn2_hp(Te_scalar)
    print(f"Results for Te={Te_scalar} eV:")
    print(f"  <sigma V>: {result_scalar['sigma_v']} m^3/s")
    print(f"  E0_ave: {result_scalar['E0_ave']} eV")
    print(f"  E0_min: {result_scalar['E0_min']} eV")
    print(f"  E0_max: {result_scalar['E0_max']} eV")

    # Test with an array of temperatures
    Te_array = np.logspace(-2, 5, 1000)  # Temperatures from 0.1 eV to 1e4 eV
    result_array = sigma_v_p_hn2_hp(Te_array)

    # Plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(Te_array, result_array['sigma_v'], label=r'$e+H_2^{+}\rightarrow e+H^{+}+H^{*}(n=2)$')
    plt.axvline(0.1, color='k', linestyle='--', alpha=0.4, label='Temperature Bounds')
    plt.axvline(2.01e4, color='k', linestyle='--', alpha=0.4)
    plt.title(r'Cross-Section for $H_2^+$ Dissociation: $e+H_2^{+}\rightarrow e+H^{+}+H^{*}(n=2)$', fontsize=14)
    plt.xlabel('Electron Temperature (eV)', fontsize=12)
    plt.ylabel(r'$<\sigma V>$ ($m^3/s$)', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()