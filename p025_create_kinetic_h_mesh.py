import numpy as np
from scipy.interpolate import interp1d

from p000_variables             import VARIABLES
from p002_Create_vr_vx_Mesh     import create_vr_vx_mesh
from p004_collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from p005_jhs_coef              import JHS_Coef
from p006_sigmav_ion_h0         import sigma_v_ion_h0
from p008_sigmav_cx_h0          import sigma_v_cx_h0


def create_kinetic_H_mesh(    
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
    Creates a refined mesh for hydrogen kinetic calculations.

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

    Outputs:
        Tuple: A tuple with the following values:
            - xH (np.ndarray): New refined spatial mesh.
            - TiH (np.ndarray): Interpolated ionic temperature on the new mesh.
            - TeH (np.ndarray): Interpolated electronic temperature on the new mesh.
            - neH (np.ndarray): Interpolated density on the new mesh.
            - PipeDiaH (np.ndarray): Interpolated tube diameter on the new mesh.
            - vx (np.ndarray): Velocities in the x-direction.
            - vr (np.ndarray): Radial velocities.
            - Tnorm (float): Normalized temperature.
            - ixE0 (Optional[int]): Index corresponding to E0 in the velocity mesh vx.
            - irE0 (Optional[int]): Index corresponding to E0 in the radial velocity mesh vr.
    """
    # Initialize optional indices
    ixE0, irE0 = None, None

    # Physical constants
    mH = VARIABLES.mH       # Proton mass [kg]
    q  = VARIABLES.q        # Elementary charge [C]
        
    # Number of points in the original mesh			
    nx = len(x)

    # These definitions have NOT been used; in fact, they were redefined later
    # # # Estimate interaction rate with side walls
    # # gamma_wall = np.zeros(nx, dtype=np.float64)

    # # # Estimate total reaction rate for destruction of hydrogen atoms and for interaction with side walls
    RR = n*sigma_v_ion_h0(Te)

    # Set v0 to thermal speed of 10 eV neutral
    v0 = np.sqrt((2*10*q)/(mu*mH))
    
    # Calculate Y from RR and v0
    y = np.zeros(nx, dtype=np.float64)
    for k in range(1,nx):
        y[k] = y[k-1] - 0.5*(x[k]-x[k-1])*(RR[k]+RR[k-1])/v0
    
    # Find x location where y = -5, i.e., where nH should be down by exp(-5)
    f_interp = interp1d(y, x, kind='linear', fill_value="extrapolate")
    xmaxH_interpolated = f_interp(-5.0)
    xmaxH = min(xmaxH_interpolated, max(x))
    xminH = x[0]


    # Interpolate Ti and Te onto a fine mesh between xminH and xmaxH
    xfine = np.linspace(xminH, xmaxH, 1001)

    Tifine      = interp1d(x, Ti,       kind='linear', fill_value="extrapolate")(xfine)
    Tefine      = interp1d(x, Te,       kind='linear', fill_value="extrapolate")(xfine)
    nfine       = interp1d(x, n,        kind='linear', fill_value="extrapolate")(xfine)
    PipeDiafine = interp1d(x, PipeDia,  kind='linear', fill_value="extrapolate")(xfine)

    # Setup a vx,vr mesh based on raw data to get typical vx, vr values
    vx, vr, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nv, Tifine,E0=E0)
    vth     = np.sqrt((2*q*Tnorm)/(mu*mH))
    min_Vr  = vth*min(vr)
    minE0   = 0.5*mH*min_Vr*min_Vr/q

    # Estimate interaction rate with side walls
    nxfine = len(xfine)
    gamma_wall = np.zeros_like(xfine)

    max_vr  = np.max(vr)
    for k in range(0,nxfine):
        if PipeDiafine[k] > 0:
            gamma_wall[k] = (2*max_vr*vth)/PipeDiafine[k]
    # # gamma_wall = np.where(PipeDiafine > 0, (2 * max_vr * vth) / PipeDiafine, 0)  
    # # Using np.where instead of the loop doesn't significantly reduce computation time,  
    # # so I chose to keep the loop.
    # print(gamma_wall)
    
    # Estimate total reaction rate, including charge exchange and elastic scattering, and interaction with side walls
    if Use_collrad_ionization:
        ioniz_rate = collrad_sigmav_ion_h0(nfine,Tefine)[0]
    else:
        if JH:
            ioniz_rate = JHS_Coef(nfine,Tefine)
        else:
            ioniz_rate = sigma_v_ion_h0(Tefine)

    minE0_array = np.full_like(xfine, minE0)
    RR = (  (nfine * ioniz_rate) +
            (nfine * sigma_v_cx_h0(Tifine, minE0_array)) +
            (gamma_wall)
        )  

    # Compute local maximum grid spacing from dx_max = 2 min(vr) / RR
    big_dx = 0.02 * fctr
    dx_max = np.minimum(fctr * 0.8 * (2 * vth * np.min(vr) / RR), big_dx)

    # Construct xH axis 
    xpt = xmaxH
    xH = [xpt]
    dxh_max=0.0005
    f_dx_max = interp1d(xfine, dx_max, kind='linear', fill_value="extrapolate")

    while xpt > xminH:
        dxpt1 = f_dx_max(xpt)
        dxpt2 = dxpt1
        xpt_test = xpt - dxpt1
        
        if xpt_test > xminH:
            dxpt2 = f_dx_max(xpt_test)
        dxpt = min(dxpt1, dxpt2, dxh_max)
        xpt  = xpt - dxpt
        xH.append(xpt)        

    xH = xH[::-1]
    xH = [xminH] + xH[:-1]

    TiH      = interp1d(xfine, Tifine,      kind='linear', fill_value="extrapolate")(xH)
    TeH      = interp1d(xfine, Tefine,      kind='linear', fill_value="extrapolate")(xH)
    neH      = interp1d(xfine, nfine,       kind='linear', fill_value="extrapolate")(xH)
    PipeDiaH = interp1d(xfine, PipeDiafine, kind='linear', fill_value="extrapolate")(xH)


    vx, vr, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nv, TiH,E0=E0)
    print(len(xH),len(x))
    print(len(TiH),len(Ti))
    return xH, TiH, TeH, neH, PipeDiaH, vx, vr, Tnorm, ixE0, irE0

if __name__=='__main__':
    import time
    from scipy.io import readsav

    start_time = time.time()
    data_file = readsav('1090904024_950to1050.sav')

    # nv = 20
    mu = data_file['mu']
    x = data_file['x']
    Ti = data_file['t_i']
    Te = data_file['t_e']
    n = data_file['n_e']
    PipeDia = data_file['d_pipe']

    # Function call
    xH, TiH, TeH, neH, PipeDiaH, vx, vr, Tnorm, ixE0, irE0 = create_kinetic_H_mesh(
        mu=mu,
        x=x,
        Ti=Ti,
        Te=Te,
        n=n,
        PipeDia=PipeDia
    )
    import matplotlib.pyplot as plt

    # Create scatter plot
    scatter1 = plt.scatter(x, x, s=10, label='x')
    scatter2 = plt.scatter(xH, xH, s=10, label='xH')

    # Get automatically assigned colors by scatter
    color_x = scatter1.get_facecolor()[0]  # Color of the first scatter
    color_xH = scatter2.get_facecolor()[0]  # Color of the second scatter

    # Add grids at the positions of x and xH with the same color as the points
    plt.vlines(x, ymin=min(x), ymax=max(x), colors=color_x, linestyles='dotted', alpha=0.5)
    plt.hlines(x, xmin=min(x), xmax=max(x), colors=color_x, linestyles='dotted', alpha=0.5)

    plt.vlines(xH, ymin=min(xH), ymax=max(xH), colors=color_xH, linestyles='dotted', alpha=0.5)
    plt.hlines(xH, xmin=min(xH), xmax=max(xH), colors=color_xH, linestyles='dotted', alpha=0.5)

    # Show legend and plot
    plt.legend()

    print("Tnorm:", Tnorm)
    end_time = time.time()
    print(end_time-start_time)
    plt.show()

    fig, axs = plt.subplots(2, 2, figsize=(10,8))  # Create a figure with a 2x2 grid

    axs[0, 0].plot(xH, TiH, label="TiH", color="b",alpha=0.5)
    axs[0, 0].plot(x, Ti, label="TiH", color="k",linestyle='--',linewidth=0.5,alpha=0.9)
    axs[0, 0].set_title("TiH")
    axs[0, 0].set_xlabel("xH")
    axs[0, 0].set_ylabel("Values")
    axs[0, 0].grid(True)

    axs[0, 1].plot(xH, TeH, label="TeH", color="r",alpha=0.5)
    axs[0, 1].plot(x, Te, label="TeH", color="k",linestyle='--',linewidth=0.5,alpha=0.9)
    axs[0, 1].set_title("TeH")
    axs[0, 1].set_xlabel("xH")
    axs[0, 1].set_ylabel("Values")
    axs[0, 1].grid(True)

    axs[1, 0].plot(xH, neH, label="neH", color="g",alpha=0.5)
    axs[1, 0].plot(x, n, label="neH", color="k",linestyle='--',linewidth=0.5,alpha=0.9)
    axs[1, 0].set_title("neH")
    axs[1, 0].set_xlabel("xH")
    axs[1, 0].set_ylabel("Values")
    axs[1, 0].grid(True)

    axs[1, 1].plot(xH, PipeDiaH, label="PipeDiaH", color="m",alpha=0.5)
    axs[1, 1].plot(x, PipeDia, label="PipeDiaH", color="k",linestyle='--',linewidth=0.5,alpha=0.9)
    axs[1, 1].set_title("PipeDiaH")
    axs[1, 1].set_xlabel("xH")
    axs[1, 1].set_ylabel("Values")
    axs[1, 1].grid(True)

    plt.tight_layout()  # Automatically adjust spaces between subplots
    plt.show()


    # Print results
    # print("xH:", xH)