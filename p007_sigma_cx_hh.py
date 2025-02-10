import numpy as np

def sigma_cx_hh(E):
    """
    Returns charge exchange cross section for molecular hydrogen. Data are taken from the polynomial fit in:
    Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.253.

    Parameters:
        E (float or array-like): Energy of the molecule corresponding to the relative velocity between the
                                 molecule and the molecular ion (eV). Must be within the range [0.1, 2e4] eV.

    Returns:
        float or numpy.ndarray: Charge exchange cross section in m^-2.
                                For energies outside the valid range, the boundary values are returned.
    """
    # Ensure E is a NumPy array and clip to valid energy range
    E = np.asarray(E, dtype=np.float64)
    E = np.clip(E, 0.1, 2.01e4)

    # Coefficients for the polynomial expansion
    alpha = np.array([
        -3.427958758517e+01,        -7.121484125189e-02,         4.690466187943e-02,
        -8.033946660540e-03,        -2.265090924593e-03,        -2.102414848737e-04,
         1.948869487515e-04,        -2.208124950005e-05,         7.262446915488e-07
    ])

    # Compute the charge exchange cross section using the polynomial expansion
    sigma_cx = np.exp(np.polyval(alpha[::-1], np.log(E))) * 1e-4

    return sigma_cx

# Example usage
if __name__ == '__main__':
    # Test with a scalar value
    E_scalar = 100  # Example energy in eV
    print(f"Charge exchange cross section for E={E_scalar} eV: {sigma_cx_hh(E_scalar)} m^-2")

    # Test with an array of energies
    E_array = np.logspace(-2, 5, 1000)  # Energies from 0.1 eV to 1e4 eV
    result = sigma_cx_hh(E_array)

    # Plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(E_array, result, label=r'$H_2^{+}+H_2\rightarrow H_2+H_2^{+}$')
    plt.axvline(0.1, color='k', linestyle='--', alpha=0.4, label='Energy Bounds')
    plt.axvline(2.01e4, color='k', linestyle='--', alpha=0.4)
    plt.title(r'Charge Exchange Cross Section for Molecular Hydrogen:$H_2^{+}+H_2\rightarrow H_2+H_2^{+}$', fontsize=14)
    plt.xlabel('Energy (eV)', fontsize=12)
    plt.ylabel(r'$\sigma_{CX}$ ($m^{-2}$)', fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()