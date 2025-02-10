import numpy as np

def sigma_v_h1s_h1s_hh(Te):
    """
    Returns Maxwellian-averaged <sigma V> for electron impact dissociation of molecular hydrogen,
    resulting in two H atoms in the 1s state. Coefficients are taken from Janev,
    "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.259.

    Also returns minimum, maximum, and average energy of the resultant H(1s) atoms.

    Parameters:
        Te (float or array-like): Electron temperature in eV. Must be within the range [0.1, 2e4] eV.

    Returns:
        dict: A dictionary containing:
            - 'sigma_v': Maxwellian-averaged <sigma V> in m^3/s.
            - 'E0_ave': Average energy of H(1s) atoms (eV).
            - 'E0_min': Minimum energy of H(1s) atoms (eV).
            - 'E0_max': Maximum energy of H(1s) atoms (eV).
    """
    # Ensure Te is a NumPy array and clip to valid temperature range
    Te = np.asarray(Te, dtype=np.float64)
    Te = np.clip(Te, 0.1, 2.01e4)

    # Coefficients for the polynomial expansion (Data from ¿?)
    b = np.array([
        -2.787217511174e+1,         1.052252660075e+1,        -4.973212347860e+0,
         1.451198183114e+0,        -3.062790554644e-1,         4.433379509258e-2,
        -4.096344172875e-3,         2.159670289222e-4,        -4.928545325189e-6
    ])
    # Data from  p. 259, Reaction 2.2.5
    # b = np.array([
    #     -2.858072836568e+01,  1.038543976082e+01, -5.383825026583e+00,
    #      1.950636494405e+00, -5.393666392407e-01,  1.006916814453e-01,
    #     -1.160758573972e-02,  7.411623859122e-04, -2.00136918807e-05
    # ])

    # Compute the Maxwellian-averaged <sigma V> using the polynomial expansion
    sigma_v = np.exp(np.polyval(b[::-1], np.log(Te))) * 1e-6

    # Energy values for the resultant particles (constant as per the IDL code)
    E0_ave = 3.0   # Average energy (eV)
    E0_min = 2.0   # Minimum energy (eV)
    E0_max = 4.25  # Maximum energy (eV)

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
    result_scalar = sigma_v_h1s_h1s_hh(Te_scalar)
    print(f"Results for Te={Te_scalar} eV:")
    print(f"  <sigma V>: {result_scalar['sigma_v']} m^3/s")
    print(f"  E0_ave: {result_scalar['E0_ave']} eV")
    print(f"  E0_min: {result_scalar['E0_min']} eV")
    print(f"  E0_max: {result_scalar['E0_max']} eV")

    # Test with an array of temperatures
    Te_array = np.logspace(-2, 5, 1000)  
    result_array = sigma_v_h1s_h1s_hh(Te_array)

    # Plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(Te_array, result_array['sigma_v'], label=r'$<\sigma V>$ for $H_2$ dissociation')
    plt.axvline(0.1, color='k', linestyle='--', alpha=0.4, label='Temperature Bounds')
    plt.axvline(2.01e4, color='k', linestyle='--', alpha=0.4)
    plt.title('Maxwellian-Averaged Cross-Section for $H_2$ Dissociation', fontsize=14)
    plt.xlabel('Electron Temperature (eV)', fontsize=12)
    plt.ylabel(r'$<\sigma V>$ ($m^3/s$)', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    # plt.ylim([1e-16,1e-6])
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show() 