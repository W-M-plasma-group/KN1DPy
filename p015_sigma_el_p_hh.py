import numpy as np

def sigma_el_p_hh(E):
    """
    Returns momentum transfer cross section for elastic collisions of H⁺ onto H₂
    for specified energy of H⁺. Data are taken from:

    Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 -
    Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and
    Helium Ions, Neutrals, and Isotopes by D.R. Schultz, S. Yu. Ovchinnikov, and S.V.
    Passovets, page 305.

    Parameters:
    -----------
    E : float or np.ndarray
        Energy of H⁺ ion (target H₂ molecule is at rest), in eV.

    Returns:
    --------
    result : float or np.ndarray
        Momentum transfer cross section (Sigma) in m² for 0.03 < E < 1e4.
        For E outside this range, the value of Sigma at the 0.03 or 1e4 eV boundary is returned.
    """
    # Convert input to a numpy array and ensure it is of type float64
    E = np.asarray(E, dtype=np.float64)
    
    # Clip energy values to the range [0.03, 1.01e4] eV
    E = np.clip(E, 0.03, 1.01e4)
    
    # Coefficients for the polynomial fit (reversed order for np.polyval)
    a = np.array([-3.355719e1, -5.696568e-1, -4.089556e-2, -1.143513e-2, 5.926596e-4])
    
    # Calculate the cross section using the polynomial fit
    result = np.exp(np.polyval(a[::-1], np.log(E))) * 1e-4
    
    return result

if __name__ == '__main__':
    # E as a scalar
    E = 50
    print(sigma_el_p_hh(E))

    # E as an array
    E = np.logspace(-2, 5, 1000)
    result = sigma_el_p_hh(E)

    import matplotlib.pyplot as plt
    plt.plot(E, result, label=r'$H^{+} - H_2$ Elastic Cross Section')
    plt.title(r'Elastic Collision: $H^{+} - H_2$')
    plt.axvline(0.03, color='k', linestyle='--', alpha=0.4)#, label='Lower Bound (0.03 eV)')
    plt.axvline(1.01e4, color='k', linestyle='--', alpha=0.4)#, label='Upper Bound (1.01e4 eV)')
    plt.xlabel('Energy (eV)')
    plt.ylabel(r'Elastic Scattering Cross-Section (m$^2$)')
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which='both')
    plt.legend()
    plt.show()