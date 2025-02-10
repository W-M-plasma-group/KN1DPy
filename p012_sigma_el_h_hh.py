import numpy as np

def sigma_el_h_hh(E: np.ndarray) -> np.ndarray:
    """
    Compute the momentum transfer cross-section for elastic collisions of H with H2.

    This function calculates the cross-section for elastic scattering between a hydrogen atom (H)
    and a hydrogen molecule (H2) as a function of the incident hydrogen atom's energy (E).

    Parameters:
    -----------
    E : float or numpy.ndarray
        Energy of the incident hydrogen atom in eV. If a scalar is provided, it will be
        converted to a numpy array.

    Returns:
    --------
    numpy.ndarray or float
        The momentum transfer cross-section in m². If the input `E` is a scalar, the output
        will also be a scalar. Otherwise, it will be an array of the same shape as `E`.

    Notes:
    ------
    - The energy values (E) are clipped to the range [0.03, 1.01e4] eV to ensure the validity
      of the empirical formula.
    - The empirical coefficients are taken from:
        Janev, R. K., et al. (1995). "Atomic and Molecular Processes in Fusion Edge Plasmas",
        Chapter 11, page 305.
    - The formula used is a polynomial expansion in log(E), with coefficients stored in the list `a`.
    """
    # Convert E to a numpy array and clip to the valid range
    E = np.asarray(E, dtype=np.float64)
    E = np.clip(E, 0.03, 1.01e4)

    # Coefficients for the polynomial expansion
    a = np.flip(np.array([
        -3.495671e1, -4.062257e-1, -3.820531e-2, -9.404486e-3, 3.963723e-4
    ]))

    # Compute the cross-section using the polynomial expansion
    result = np.exp(np.polyval(a, np.log(E))) * 1e-4

    return result

if __name__ == '__main__':
    E_scalar = 100.0
    sigma_scalar = sigma_el_h_hh(E_scalar)
    print(f"Sección transversal de transferencia de momento (E = {E_scalar} eV): {sigma_scalar} m²")

    E = np.logspace(-2, 5, 100)  # Desde 0.1 eV hasta 10,000 eV
    sigma = sigma_el_h_hh(E)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    # Plot
    plt.plot(E, sigma, label=r'$H + H_2 \rightarrow H + H_2$ (Momentum Transfer)')
    plt.axvline(0.03  ,color='k',linestyle='--',alpha=0.4)
    plt.axvline(1.01e4,color='k',linestyle='--',alpha=0.4)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Energía del átomo de hidrógeno (eV)')
    plt.ylabel(r'$\sigma$ (m$^2$)')
    plt.title('Sección transversal de transferencia de momento para $H + H_2$')
    plt.legend()
    plt.show()