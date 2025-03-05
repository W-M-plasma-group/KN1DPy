import numpy as np

def sigma_el_h_h(E, vis=False):
    """
    Compute the elastic scattering cross-section for hydrogen-hydrogen collisions.

    This function calculates either the momentum transfer cross-section or the viscosity
    cross-section for elastic scattering between hydrogen atoms, depending on the `vis` flag.

    Parameters:
    -----------
    E : float or numpy.ndarray
        Energy of the incident hydrogen atom in eV. If a scalar is provided, it will be
        converted to a numpy array.
    vis : bool, optional
        If True, the viscosity cross-section is returned. If False (default), the momentum
        transfer cross-section is returned.

    Returns:
    --------
    numpy.ndarray or float
        The elastic scattering cross-section in m². If the input `E` is a scalar, the output
        will also be a scalar. Otherwise, it will be an array of the same shape as `E`.

    Notes:
    ------
    - The energy values (E) are clipped to the range [0.03, 1.01e4] eV to ensure the validity
      of the empirical formula.
    - The empirical coefficients are taken from:
        - Momentum transfer cross-section: Unknown source (similar to Janev or Freeman).
        - Viscosity cross-section: Unknown source (similar to Janev or Freeman).
    """
    # Convert E to a numpy array and clip to the valid range
    E = np.asarray(E, dtype=np.float64)
    E = np.clip(E, 0.03, 1.01e4)
    print(E)
    

    if vis: ### vis = True
        # Coefficients for the viscosity cross-section
        a = np.array([-3.344860e1, -4.238982e-1, -7.477873e-2, -7.915053e-3, -2.686129e-4])
    else:   ### vis = False
        # Coefficients for the momentum transfer cross-section
        a = np.array([-3.330843e1, -5.738374e-1, -1.028610e-1, -3.920980e-3, 5.964135e-4])

    # Compute the cross-section using the polynomial expansion
    result = np.exp(np.polyval(a[::-1], np.log(E))) * 1e-4
    return result


if __name__ == '__main__':
    # Ejemplo con E como escalar
    E_scalar = 100.0  # Energía en eV
    # Calcular la sección transversal de transferencia de momento
    sigma_momentum_scalar = sigma_el_h_h(E_scalar)
    print(f"Sección transversal de transferencia de momento (E = {E_scalar} eV): {sigma_momentum_scalar} m²")
    # Calcular la sección transversal de viscosidad
    sigma_viscosity_scalar = sigma_el_h_h(E_scalar, vis=True)
    print(f"Sección transversal de viscosidad (E = {E_scalar} eV): {sigma_viscosity_scalar} m²")

    # Ejemplo con E como array
    E_array = np.logspace(-2,5,100)  # Energías en eV
    # Calcular la sección transversal de transferencia de momento
    sigma_momentum_array = sigma_el_h_h(E_array)
    # Calcular la sección transversal de viscosidad
    sigma_viscosity_array = sigma_el_h_h(E_array, vis=True)
    #------------------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    plt.plot(E_array,sigma_momentum_array,label='Momentum')
    plt.plot(E_array,sigma_viscosity_array,label='Viscosity')
    plt.title(r'Elastic collision: $H-H$')
    plt.axvline(0.03  ,color='k',linestyle='--',alpha=0.4)
    plt.axvline(1.01e4,color='k',linestyle='--',alpha=0.4)
    plt.xlabel('Energy (eV)')
    plt.ylabel(r'$\langle\sigma\rangle$ ($m^2$)')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True,'both')
    plt.legend()
    plt.show()