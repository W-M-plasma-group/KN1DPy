from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

from .utils import get_config, interp_1d, reverse
from .rates.janev.sigmav_ion_h0 import sigmav_ion_h0
from .rates.janev.sigmav_cx_h0 import sigmav_cx_h0
from .rates.janev.sigmav_ion_hh import sigmav_ion_hh
from .rates.janev.sigmav_h1s_h1s_hh import sigmav_h1s_h1s_hh
from .rates.janev.sigmav_h1s_h2s_hh import sigmav_h1s_h2s_hh
from .rates.janev.sigmav_cx_hh import sigmav_cx_hh
from .rates.collrad.collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from .rates.johnson_hinnov.johnson_hinnov import Johnson_Hinnov

from .common import constants as CONST


class KineticMesh:
    '''
    Mesh data for kinetic_neutrals procedure

    Attributes
    ----------
        mesh_type : str
            Type of mesh data
            - 'h' for atomic
            - 'h2' for molecular
        x : ndarray
            spatial coordinates
        vx : ndarray
            axial velocity coordinates
        vr : ndarray
            radial velocity coordinates
        Ti : ndarray
            ion temperature profile interpolated over x
        Te : ndarray
            electron temperature profile interpolated over x
        ne : ndarray
            density profile interpolated over x  
        PipeDia : ndarray
            effective pipe diameter interpolated over x
        Tnorm : float
            Average ion temperature
    '''

    def __init__( #NOTE Simplify this later, consider using class inheritance
            self,
            mesh_type   : str, #'h' for kinetic_h_mesh, 'h2' for kinetic_h2_mesh
            mu          : int,
            x           : NDArray,
            Ti          : NDArray,
            Te          : NDArray,
            n           : NDArray,
            PipeDia     : NDArray,
            jh          : Johnson_Hinnov = None,
            E0          : NDArray = np.array([0.0]),
            fctr        : float   = None,
            config_path : str     = './config.json'):

        print("generating kinetic_" + mesh_type + "_mesh")

        #Get mesh size from config file
        cfg = get_config(config_path)
        nv = cfg["kinetic_" + mesh_type]["mesh_size"]
        if fctr is None:
            fctr = cfg["kinetic_" + mesh_type]["grid_fctr"]

        # estimate Interaction rate with side walls
        #NOTE Commented gamma_wall calculations here, revisit later

        if mesh_type == 'h':
            # Estimate total reaction rate for destruction of hydrogen atoms and for interation with side walls
            react_rate = n*sigmav_ion_h0(Te) 
            # Set v0 to thermal speed to 10 eV neutral 
            v0 = np.sqrt(2*10*CONST.Q / (mu*CONST.H_MASS))

        elif mesh_type == 'h2':
            #Estimate total reaction rate for destruction of molecules and for interation with side walls
            react_rate = n*sigmav_ion_hh(Te) + n*sigmav_h1s_h1s_hh(Te) + n*sigmav_h1s_h2s_hh(Te)
            #directed random velocity of diatomic molecule
            v0 = np.sqrt(8.0*CONST.TWALL*CONST.Q / (np.pi*2*mu*CONST.H_MASS))

        else:
            raise Exception("ERROR: Mesh type invalid:", mesh_type)


        # Determine x range for atoms by finding distance into plasma where density persists.
        y = np.zeros(np.size(x), float)
        for k in range(1, np.size(x)): 
            y[k] = y[k-1] - ((x[k] - x[k-1])*0.5*(react_rate[k] + react_rate[k-1]))/v0
        if mesh_type == 'h':
            # Find x location where Y = -5, i.e. where nH should be down by exp(-5)
            expdown = max(-5, np.min(y))
            xmax = np.minimum(interp_1d(y, x, expdown, fill_value="extrapolate"), max(x))
        elif mesh_type == 'h2':
            #Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
            xmax = np.minimum(interp_1d(y, x, -10.0), max(x))
        xmin = x[0]


        # Interpolate Ti and Te onto a fine mesh between xmin and xmax 
        xfine = xmin + (xmax - xmin)*np.arange(1001)/1000

        Tifine = interp_1d(x, Ti, xfine, fill_value="extrapolate")
        Tefine = interp_1d(x, Te, xfine)
        nfine = interp_1d(x, n, xfine)
        PipeDiafine = interp_1d(x, PipeDia, xfine)


        # Set up a vx, vr mesh based on raw data to get typical vx, vr values 
        vx, vr, Tnorm = self.create_vr_vx_mesh(nv, Tifine, E0=E0)

        vth = np.sqrt( (2*CONST.Q*Tnorm) / (mu*CONST.H_MASS))
        # Estimate interaction rate with side walls
        nxfine = np.size(xfine)
        gamma_wall = np.zeros(nxfine, float)
        for k in range(nxfine):
            if PipeDiafine[k] > 0:
                gamma_wall[k] = 2 * max(vr) * vth / PipeDiafine[k]
        
        # Estimate total reaction rate, including charge exchange and elastic scattering, and interaction with side walls 
        if mesh_type == 'h':
            minVr = vth*min(vr)
            minE0 = 0.5*CONST.H_MASS*(minVr**2) / CONST.Q

            # Hardwired to JH for comparison test with IDL
            if jh is None:
                jh = Johnson_Hinnov()
            ioniz_rate = jh.jhs_coef(nfine, Tefine, no_null=True)
            react_rate = nfine*(ioniz_rate + sigmav_cx_h0(Tifine, np.full(xfine.shape, minE0))) + gamma_wall

        elif mesh_type == 'h2':
            react_rate = nfine*(sigmav_ion_hh(Tefine) + sigmav_h1s_h1s_hh(Tefine) + sigmav_h1s_h2s_hh(Tefine) + 0.1*sigmav_cx_hh(Tifine,Tifine)) + gamma_wall
        
        # Compute local maximum grid spacing dx_max = 2
        dx_max = np.minimum(fctr*0.8*(2*vth*min(vr)/react_rate), 0.02*fctr)

        # Construct xH Axis 
        xpt = xmax
        xH = np.array([xpt])
        
        while xpt > xmin:
            xH = np.concatenate([np.array([xpt]), xH])
            dxpt1 = interp_1d(xfine, dx_max, xpt, fill_value="extrapolate")
            dxpt2 = np.copy(dxpt1)
            xpt_test = xpt - dxpt1
            if xpt_test > xmin:
                dxpt2 = interp_1d(xfine, dx_max, xpt_test, fill_value="extrapolate")

            dxpt = min([dxpt1, dxpt2])

            xpt -= dxpt 
        xH = np.concatenate([np.array([xmin]), xH[0:np.size(xH) - 1]])


        TiH = np.interp(xH, xfine, Tifine)
        TeH = interp_1d(xfine, Tefine, xH)
        neH = interp_1d(xfine, nfine, xH)
        PipeDiaH = interp_1d(xfine, PipeDiafine, xH)

        vx, vr, Tnorm = self.create_vr_vx_mesh(nv, TiH, E0=E0)


        self.mesh_type : str = mesh_type
        self.x : NDArray = xH
        self.Ti : NDArray = TiH
        self.Te : NDArray = TeH
        self.ne : NDArray = neH
        self.PipeDia : NDArray = PipeDiaH
        self.vx : NDArray = vx
        self.vr : NDArray = vr
        self.Tnorm : float = Tnorm


    def create_vr_vx_mesh(self, nv: int, Ti: NDArray, E0: NDArray = np.array([0.0]), Tmax: float = 0.0) -> tuple[NDArray, NDArray, float] :
        # Gwendolyn Galleher 
        '''
        Sets up optimum Vr and Vx velocity space mesh for Kinetic_Neutrals procedure 

        Parameters
        ----------
            nv : int
                number of elements desired in vr mesh
            Ti : ndarray
                ion temperature profile
            E0 : ndarray
                energy where a velocity is desired (optional)
            Tmax : float
                maximum temperature, ignore Ti above this value
                
        Returns
        -------
            vr: ndarray
                radial velocities
            vx: ndarray
                axial velocities
            Tnorm
                average of Ti
        '''

        Ti = np.array(Ti) 
        Ti = np.concatenate([Ti, E0[E0>0]])
        if Tmax > 0:
            ii = np.where(Ti < Tmax)
            Ti = Ti[ii]
        
        maxTi = Ti.max()
        minTi = Ti.min()
        Tnorm = np.nanmean(Ti)
        vmax = 3.5
        if (maxTi-minTi) <= (0.1*maxTi):
            v = (np.arange(nv+1)*vmax) / nv
        else:
            g = 2*nv*np.sqrt(minTi/maxTi) / (1 - np.sqrt(minTi/maxTi))
            b = vmax / (nv*(nv + g))
            v = (g*b)*np.arange(nv+1) + b*(np.arange(nv+1)**2)

        # Option: add velocity bins corresponding to E0     
        v0 = 0
        for k in range(np.size(E0)):
            if E0[k] > 0.0:
                v0 = np.sqrt(E0[k]/Tnorm)
                ii = np.argwhere(v > v0).T[0]
                if np.size(ii) > 0:
                    v = np.concatenate([v[0:ii[0]], [v0], v[ii[0]:]])
                else: 
                    v = np.concatenate([v, v0])
            
        vr = v[1:]
        vx = np.concatenate([-reverse(vr), vr]) 

        return vx,vr,Tnorm


    #Setup string conversion for printing
    def __str__(self):
        string = "Kinetic Mesh:\n"
        string += "    x: " + str(self.x) + "\n"
        string += "    Ti: " + str(self.Ti) + "\n"
        string += "    Te: " + str(self.Te) + "\n"
        string += "    ne: " + str(self.ne) + "\n"
        string += "    PipeDia: " + str(self.PipeDia) + "\n"
        string += "    vx: " + str(self.vx) + "\n"
        string += "    vr: " + str(self.vr) + "\n"
        string += "    Tnorm: " + str(self.Tnorm) + "\n"
        return string