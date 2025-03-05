import numpy as np

def sigma_el_p_h(E: np.ndarray) -> np.ndarray:
    """
    Compute the momentum transfer cross-section for elastic collisions of H+ with H.

    This function calculates the cross-section for elastic scattering between a proton (H+)
    and a hydrogen atom (H) as a function of the incident proton's energy (E).

    Parameters:
    -----------
    E : float or numpy.ndarray
        Energy of the incident proton in eV. If a scalar is provided, it will be
        converted to a numpy array.

    Returns:
    --------
    numpy.ndarray or float
        The momentum transfer cross-section in m². If the input `E` is a scalar, the output
        will also be a scalar. Otherwise, it will be an array of the same shape as `E`.

    Notes:
    ------
    - The energy values (E) are clipped to the range [0.001, 1.01e5] eV to ensure the validity
      of the empirical formula.
    - The empirical coefficients are taken from:
        Janev, R. K., et al. (1995). "Atomic and Molecular Processes in Fusion Edge Plasmas",
        Chapter 11, page 298.
    - The formula uses two sets of coefficients depending on the energy range:
        - For E ≤ 10.0 eV: A 7-term polynomial expansion in log(E).
        - For E > 10.0 eV: A 2-term polynomial expansion in log(E).
    """
    # Convert E to a numpy array and clip to the valid range
    E = np.asarray(E, dtype=np.float64)
    E = np.clip(E, 0.001, 1.01e5)

    # Initialize the result array with zeros
    result = np.zeros_like(E)

    # Coefficients for the low-energy range (E ≤ 10.0 eV)
    a_low = np.array([
        -3.233966e1, -1.126918e-1, 5.287706e-3, -2.445017e-3, -1.044156e-3, 8.419691e-5, 3.824773e-5
        ])

    # Coefficients for the high-energy range (E > 10.0 eV)
    a_high = np.array([
        -3.231141e1, -1.386002e-1
        ])

    # Compute the cross-section for the low-energy range
    mask_low = E <= 10.0
    if np.any(mask_low):
        result[mask_low] = np.exp(np.polyval(a_low[::-1], np.log(E[mask_low]))) * 1e-4

    # Compute the cross-section for the high-energy range
    mask_high = E > 10.0
    if np.any(mask_high):
        result[mask_high] = np.exp(np.polyval(a_high[::-1], np.log(E[mask_high]))) * 1e-4

    return result

if __name__ == '__main__':
    # E as a scalar
    E_scalar_low = 5.0  # Energía en el rango bajo (E ≤ 10.0 eV)
    E_scalar_high = 100.0  # Energía en el rango alto (E > 10.0 eV)

    # Calcular la sección transversal
    sigma_low = sigma_el_p_h(E_scalar_low)
    sigma_high = sigma_el_p_h(E_scalar_high)

    print(f"Sección transversal (E = {E_scalar_low} eV): {sigma_low} m²")
    print(f"Sección transversal (E = {E_scalar_high} eV): {sigma_high} m²")
    
    # E as an array
    E = np.logspace(-3.5, 5.5, 1000)  # Desde 0.001 eV hasta 100,000 eV
    sigma = sigma_el_p_h(E)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt  
    # Plot  
    plt.plot(E, sigma, label=r'$H^+ + H \rightarrow H^+ + H$ (Momentum Transfer)')
    plt.axvline(0.001 , color='k',linestyle='--',alpha=0.4)
    plt.axvline(1.01e5, color='k',linestyle='--',alpha=0.4)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which='both', linestyle='-', linewidth=0.5)
    plt.xlabel(r'$E$ (eV)')
    plt.ylabel(r'$\langle\sigma\rangle$ (m$^2$)')
    plt.title('Sección transversal para colisiones elásticas de $H^+ + H$')
    plt.legend()
    plt.show()