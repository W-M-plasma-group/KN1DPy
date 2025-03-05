import numpy as np
from scipy.interpolate import interp1d

from p000_variables import VARIABLES
from p002_Create_vr_vx_Mesh import create_vr_vx_mesh
from p006_sigmav_ion_hh import sigma_v_ion_hh
from p008_sigmav_cx_hh import sigma_v_cx_hh
from p022_sigma_v_h1s_h1s_hh import sigma_v_h1s_h1s_hh
from p023_sigma_v_h1s_h2s_hh import sigma_v_h1s_h2s_hh


def create_kinetic_H2_mesh(
    mu: float,
    x: np.ndarray,
    Ti: np.ndarray,
    Te: np.ndarray,
    n: np.ndarray,
    PipeDia: np.ndarray,
    nv: int = 20,
    fctr: float = 1.0,
    JH: int = 0,
    Use_collrad_ionization: int = 1,
    E0: np.ndarray = np.array([0.0])
) -> tuple:
    """
    Creates a refined mesh for hydrogen molecule (H2) kinetic calculations.

    Inputs:
        mu (float): Reduced mass of the system.
        x (np.ndarray): Original spatial coordinates.
        Ti (np.ndarray): Ionic temperature at the positions given by `x`.
        Te (np.ndarray): Electronic temperature at the positions given by `x`.
        n (np.ndarray): Density at the positions given by `x`.
        PipeDia (np.ndarray): Tube diameter at the positions given by `x`.
        nv (int, optional): Number of desired elements in the radial velocity mesh. Default is 20.
        fctr (float, optional): Scaling factor to adjust the maximum grid spacing. Default is 1.0.
        JH (int, optional): Flag to choose ionization coefficients based on Johnson-Hinnov. Default is 0.
        Use_collrad_ionization (int, optional): Flag to use COLLRAD ionization data. Default is 1.
        E0 (np.ndarray, optional): Energy array for velocity mesh initialization. Default is [0.0].

    Outputs:
        Tuple: A tuple with the following values:
            - xH2 (np.ndarray): New refined spatial mesh.
            - TiH2 (np.ndarray): Interpolated ionic temperature on the new mesh.
            - TeH2 (np.ndarray): Interpolated electronic temperature on the new mesh.
            - neH2 (np.ndarray): Interpolated density on the new mesh.
            - PipeDiaH2 (np.ndarray): Interpolated tube diameter on the new mesh.
            - vx (np.ndarray): Velocities in the x-direction.
            - vr (np.ndarray): Radial velocities.
            - Tnorm (float): Normalized temperature.
            - ixE0 (Optional[int]): Index corresponding to E0 in the velocity mesh vx.
            - irE0 (Optional[int]): Index corresponding to E0 in the radial velocity mesh vr.
    """

    # Physical constants
    mH = VARIABLES.mH       # Proton mass [kg]
    q = VARIABLES.q         # Elementary charge [C]
    k_b = VARIABLES.k_b     # Boltzmann constant [J/K]

    Twall = 293.0 * (k_b / q)  # Room temperature [eV]
    v0_bar = np.sqrt(8.0 * Twall * q / (np.pi * 2 * mu * mH))  # Average thermal speed of a diatomic molecule [m/s]

    nx = len(x)  # Number of points in the original mesh

    # Estimate interaction rate with side walls
    gamma_wall = np.zeros(nx, dtype=np.float64)

    # Total reaction rate for destruction of molecules
    RR = (
        (n * sigma_v_ion_hh(Te)) +
        (n * sigma_v_h1s_h1s_hh(Te)['sigma_v']) +
        (n * sigma_v_h1s_h2s_hh(Te)['sigma_v'])
    )

    # Calculate Y from RR and v0_bar
    y = np.zeros(nx, dtype=np.float64)
    for k in range(1, nx):
        y[k] = y[k - 1] - 0.5 * (x[k] - x[k - 1]) * (RR[k] + RR[k - 1]) / v0_bar

    # Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
    f_interp = interp1d(y, x, kind='linear', fill_value="extrapolate")
    xmaxH2_interpolated = f_interp(-10.0)
    xmaxH2 = min(xmaxH2_interpolated, max(x))
    xminH2 = x[0]

    # Interpolate Ti, Te, n, and PipeDia onto a fine mesh between xminH2 and xmaxH2
    xfine = np.linspace(xminH2, xmaxH2, 1001)
    Tifine = interp1d(x, Ti, kind='linear', fill_value="extrapolate")(xfine)
    Tefine = interp1d(x, Te, kind='linear', fill_value="extrapolate")(xfine)
    nfine = interp1d(x, n, kind='linear', fill_value="extrapolate")(xfine)
    PipeDiafine = interp1d(x, PipeDia, kind='linear', fill_value="extrapolate")(xfine)

    # Setup a vx, vr mesh based on raw data to get typical vx, vr values
    vx, vr, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nv, Tifine, E0=E0)
    vth = np.sqrt((2 * q * Tnorm) / (mu * mH))  # Thermal velocity

    # Estimate interaction rate with side walls
    nxfine = len(xfine)
    gamma_wall = np.zeros_like(xfine)

    max_vr = np.max(vr)
    for k in range(0, nxfine):
        if PipeDiafine[k] > 0:
            gamma_wall[k] = (2 * max_vr * vth) / PipeDiafine[k]

    # Total reaction rate, including charge exchange, elastic scattering, and interaction with side walls
    RR = (
        (nfine * sigma_v_ion_hh(Tefine)) +
        (nfine * sigma_v_h1s_h1s_hh(Tefine)['sigma_v']) +
        (nfine * sigma_v_h1s_h2s_hh(Tefine)['sigma_v']) +
        (nfine * sigma_v_cx_hh(Tifine, Tifine) * 0.1) +
        gamma_wall
    )

    # Compute local maximum grid spacing from dx_max = 2 * min(vr) / RR
    big_dx = 0.02 * fctr
    dx_max = np.minimum(fctr * 0.8 * (2 * vth * np.min(vr) / RR), big_dx)

    # Construct xH2 axis
    xpt = xmaxH2
    xH2 = [xpt]
    f_dx_max = interp1d(xfine, dx_max, kind='linear', fill_value="extrapolate")
    # dxh_max = 0.0005   # This line was added to got more points on this new mesh
    while xpt > xminH2:
        dxpt1 = f_dx_max(xpt)
        dxpt2 = dxpt1
        xpt_test = xpt - dxpt1

        if xpt_test > xminH2:
            dxpt2 = f_dx_max(xpt_test)
        dxpt = min(dxpt1, dxpt2)#, dxh_max)
        xpt = xpt - dxpt
        xH2.append(xpt)

    xH2 = xH2[::-1]
    xH2 = [xminH2] + xH2[:-1]

    # Interpolate Ti, Te, n, and PipeDia onto the new xH2 mesh
    TiH2 = interp1d(xfine, Tifine, kind='linear', fill_value="extrapolate")(xH2)
    TeH2 = interp1d(xfine, Tefine, kind='linear', fill_value="extrapolate")(xH2)
    neH2 = interp1d(xfine, nfine, kind='linear', fill_value="extrapolate")(xH2)
    PipeDiaH2 = interp1d(xfine, PipeDiafine, kind='linear', fill_value="extrapolate")(xH2)

    # Recreate the vx, vr mesh based on the interpolated data
    vx, vr, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nv, TiH2, E0=E0)

    print(len(xH2), len(x))
    print(len(TiH2), len(Ti))
    return xH2, TiH2, TeH2, neH2, PipeDiaH2, vx, vr, Tnorm, ixE0, irE0


if __name__ == '__main__':
    import time
    from scipy.io import readsav

    start_time = time.time()
    data_file = readsav('1090904024_950to1050.sav')

    # Extract variables from the data file
    mu = data_file['mu']
    x = data_file['x']
    Ti = data_file['t_i']
    Te = data_file['t_e']
    n = data_file['n_e']
    PipeDia = data_file['d_pipe']
    E0 = [0.003,0.01,0.03,0.1,0.3,1.0,3.0]
    # Call the function
    xH2, TiH2, TeH2, neH2, PipeDiaH2, vx, vr, Tnorm, ixE0, irE0 = create_kinetic_H2_mesh(
        mu=mu,
        x=x,
        Ti=Ti,
        Te=Te,
        n=n,
        PipeDia=PipeDia,
        E0=E0
    )

    import matplotlib.pyplot as plt

    # Create scatter plot
    scatter1 = plt.scatter(x, x, s=10, label='x')
    scatter2 = plt.scatter(xH2, xH2, s=10, label='xH')

    # Get automatically assigned colors by scatter
    color_x = scatter1.get_facecolor()[0]  # Color of the first scatter
    color_xH = scatter2.get_facecolor()[0]  # Color of the second scatter

    # Add grids at the positions of x and xH with the same color as the points
    plt.vlines(x, ymin=min(x), ymax=max(x), colors=color_x, linestyles='dotted', alpha=0.5)
    plt.hlines(x, xmin=min(x), xmax=max(x), colors=color_x, linestyles='dotted', alpha=0.5)

    plt.vlines(xH2, ymin=min(xH2), ymax=max(xH2), colors=color_xH, linestyles='dotted', alpha=0.5)
    plt.hlines(xH2, xmin=min(xH2), xmax=max(xH2), colors=color_xH, linestyles='dotted', alpha=0.5)

    # Show legend and plot
    plt.legend()

    print("Tnorm:", Tnorm)
    end_time = time.time()
    print(end_time - start_time)
    plt.show()

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))  # Create a figure with a 2x2 grid

    # Plot TiH2 vs xH2 and compare with original Ti vs x
    axs[0, 0].plot(xH2, TiH2, label="TiH", color="b", alpha=0.5)
    axs[0, 0].plot(x, Ti, label="TiH", color="k", linestyle='--', linewidth=0.5, alpha=0.9)
    axs[0, 0].set_title("TiH")
    axs[0, 0].set_xlabel("xH")
    axs[0, 0].set_ylabel("Values")
    axs[0, 0].grid(True)

    # Plot TeH2 vs xH2 and compare with original Te vs x
    axs[0, 1].plot(xH2, TeH2, label="TeH", color="r", alpha=0.5)
    axs[0, 1].plot(x, Te, label="TeH", color="k", linestyle='--', linewidth=0.5, alpha=0.9)
    axs[0, 1].set_title("TeH")
    axs[0, 1].set_xlabel("xH")
    axs[0, 1].set_ylabel("Values")
    axs[0, 1].grid(True)

    # Plot neH2 vs xH2 and compare with original n vs x
    axs[1, 0].plot(xH2, neH2, label="neH", color="g", alpha=0.5)
    axs[1, 0].plot(x, n, label="neH", color="k", linestyle='--', linewidth=0.5, alpha=0.9)
    axs[1, 0].set_title("neH")
    axs[1, 0].set_xlabel("xH")
    axs[1, 0].set_ylabel("Values")
    axs[1, 0].grid(True)

    # Plot PipeDiaH2 vs xH2 and compare with original PipeDia vs x
    axs[1, 1].plot(xH2, PipeDiaH2, label="PipeDiaH", color="m", alpha=0.5)
    axs[1, 1].plot(x, PipeDia, label="PipeDiaH", color="k", linestyle='--', linewidth=0.5, alpha=0.9)
    axs[1, 1].set_title("PipeDiaH")
    axs[1, 1].set_xlabel("xH")
    axs[1, 1].set_ylabel("Values")
    axs[1, 1].grid(True)

    plt.tight_layout()  # Automatically adjust spaces between subplots
    plt.show()