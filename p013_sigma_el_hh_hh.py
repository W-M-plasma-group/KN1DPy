import numpy as np

def sigma_el_hh_hh(E, vis: bool = False) -> np.ndarray:
    """
    Compute the momentum transfer or viscosity cross-section for elastic collisions of H2 with H2.

    This function calculates the cross-section for elastic scattering between two hydrogen molecules (H2)
    as a function of the incident hydrogen molecule's energy (E). The cross-section can be either the
    momentum transfer cross-section or the viscosity cross-section, depending on the `vis` flag.

    Parameters:
    -----------
    E : float or numpy.ndarray
        Energy of the incident hydrogen molecule in eV. If a scalar is provided, it will be
        converted to a numpy array.
    vis : bool, optional
        If True, the viscosity cross-section is returned. If False (default), the momentum
        transfer cross-section is returned.

    Returns:
    --------
    numpy.ndarray or float
        The cross-section in m². If the input `E` is a scalar, the output will also be a scalar.
        Otherwise, it will be an array of the same shape as `E`.

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

    # Coefficients for the polynomial expansion (in reverse order for np.polyval)
    a = np.array([
            -3.430345e1, -2.960406e-1, -6.382532e-2, -7.557519e-3, 2.606259e-4
        ])

    # Compute the cross-section using the polynomial expansion
    result = np.exp(np.polyval(a[::-1], np.log(E))) * 1e-4

    # If vis is True, print a warning (as in the IDL code)
    if vis:
        print("WARNING in sigma_el_hh_hh => using momentum transfer as viscosity cross-section")

    return result

if __name__ == '__main__':
    # E as a scalar
    E_scalar = 100.0
    # Calcular la sección transversal de transferencia de momento
    sigma_momentum_scalar = sigma_el_hh_hh(E_scalar)
    print(f"Sección transversal de transferencia de momento (E = {E_scalar} eV): {sigma_momentum_scalar} m²")

    # Calcular la sección transversal de viscosidad
    sigma_viscosity_scalar = sigma_el_hh_hh(E_scalar, vis=True)
    print(f"Sección transversal de viscosidad (E = {E_scalar} eV): {sigma_viscosity_scalar} m²")

    # E as an array
    E = np.logspace(-2, 5, 1000)  
    # Calcular la sección transversal de transferencia de momento
    sigma_momentum = sigma_el_hh_hh(E)
    # Calcular la sección transversal de viscosidad
    sigma_viscosity = sigma_el_hh_hh(E, vis=True)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    # Plot
    plt.plot(E, sigma_momentum, label=r'$H_2 + H_2 \rightarrow H_2 + H_2$ (Momentum Transfer)')
    plt.plot(E, sigma_viscosity, label=r'$H_2 + H_2 \rightarrow H_2 + H_2$ (Viscosity)', linestyle='--')
    plt.axvline(0.03  ,color='k',linestyle='--',alpha=0.4)
    plt.axvline(1.01e4,color='k',linestyle='--',alpha=0.4)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel('Energía de la molécula de hidrógeno (eV)')
    plt.ylabel(r'$\langle\sigma\rangle$ (m$^2$)')
    plt.title('Sección transversal para colisiones elásticas de $H_2 + H_2$')
    plt.legend()
    plt.show()


