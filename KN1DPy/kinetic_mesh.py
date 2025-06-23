import numpy as np 
from numpy.typing import NDArray
from scipy import interpolate

from .sigma.sigmav_ion_h0 import sigmav_ion_h0
from .sigma.sigma_cx_h0 import sigma_cx_h0
from .sigma.sigmav_ion_hh import sigmav_ion_hh
from .sigma.sigmav_h1s_h1s_hh import sigmav_h1s_h1s_hh
from .sigma.sigmav_h1s_h2s_hh import sigmav_h1s_h2s_hh
from .sigma.sigmav_cx_hh import sigmav_cx_hh
from .collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from .jh_related.jhs_coef import jhs_coef
from .create_vr_vx_mesh import create_vr_vx_mesh

from .common import constants as CONST
from .common.JH_Coef import JH_Coef

class kinetic_mesh:

    def __init__(
            self, 
            mesh_type  : str, #'h' for kinetic_h_mesh, 'h2' for kinetic_h2_mesh
            nv         : int, 
            mu         : int, 
            x          : NDArray,
            Ti         : NDArray,
            Te         : NDArray, 
            n          : NDArray, 
            PipeDia    : NDArray,
            jh_coeffs  : JH_Coef = None,
            E0 = 0, ixE0 = 0, irE0 = 0, fctr = 1.0
        ):
        
        self.mesh_type = mesh_type

        print("generating kinetic_" + mesh_type + "_mesh")

        # add these lines to improve velocity space resolution (use with care):
        #nv = 20
        # E0[.1]

        nx = np.size(x)

        # estimate Interaction rate with side walls
        gamma_wall = np.zeros(nx) # fixed typo - GG
        # for k in range( 0, nx):
        #   if PipeDia[k] > 0:
        #        gamma_wall[k] = 2*sqrt(2*Ti(k)*CONST.Q/(2*CONST.H_MASS))/PipeDia(k)

        if mesh_type == 'h':
            # Estimate total reaction rate for destriction of hydrogen atoms and for interation with side walls
            # RR = n*sigmav_ion_H0(Te)+gamma_wall
            RR = n * sigmav_ion_h0(Te) 

            # Set v0 to thermal speed to 10 eV neutral 
            v0 = np.sqrt( 2 * 10 * CONST.Q / (mu * CONST.H_MASS))

        elif mesh_type == 'h2':
            #Estimate total reaction rate for destruction of molecules and for interation with side walls
            RR=n*sigmav_ion_hh(Te)+n*sigmav_h1s_h1s_hh(Te)+n*sigmav_h1s_h2s_hh(Te)

            #directed random velocity of diatomic molecule
            v0 = np.sqrt(8.0*CONST.TWALL*CONST.Q/(np.pi*2*mu*CONST.H_MASS))	

        else:
            raise Exception("ERROR: Mesh type invalid:", mesh_type)


        # Determine x range for atoms by finding distance into plasma where density persists.
        #  dGamma/dx=-nH*RR = v0 dnH/dx = -nH*RR - these two lines are commented in the original code I dont understand them
        #  d ln(nH)/dx = -RR/v0 = dY/dx

        Y = np.zeros(nx) # fixed typo
        for k in range(1, nx): 
            Y[k] = Y[k-1] - (x[k] - x[k-1] ) * 0.5 * (RR[k] + RR[k-1]) / v0
        
        if mesh_type == 'h':
            # Find x location where Y = -5, i.e. where nH should be down by exp(-5)
            interpfunc = interpolate.interp1d(Y, x, kind = 'linear', bounds_error=False, fill_value="extrapolate") # previous version caused error on next line
            xmaxH = np.minimum(interpfunc(-5), max(x))
        elif mesh_type == 'h2':
            #Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
            interpfunc = interpolate.interp1d(Y, x) # fixed error with interpolation - GG
            xmaxH=np.minimum(interpfunc(-10.0), max(x))

        xminH = x[0]

        # Interpolate Ti and Te onto a fine mesh between xminH and xmaxH 
        xfine = xminH + (xmaxH - xminH) * np.arange( 1001 )/1000

        interpfunc = interpolate.interp1d(x, Ti, kind = 'linear')
        Tifine = interpfunc(xfine)

        interpfunc = interpolate.interp1d(x, Te, kind = 'linear')
        Tefine = interpfunc(xfine)

        interpfunc = interpolate.interp1d(x, n, kind = 'linear')
        nfine = interpfunc(xfine)

        interpfunc = interpolate.interp1d(x, PipeDia)
        PipeDiafine = interpfunc(xfine)

        # Set up a vx, vr mesh based on raw data to get typical vx, vr values 
        vx, vr, Tnorm, ixE0, ixE0 = create_vr_vx_mesh(nv, Tifine) # fixed error from not assigning all outputs - GG
        
        vth = np.sqrt( (2 * CONST.Q * Tnorm) / (mu * CONST.H_MASS))

        # Estimate interaction rate with side walls
        nxfine = np.size(xfine)
        gamma_wall = np.zeros(nxfine)
        #print(PipeDiafine)
        for k in range(nxfine): # fixed typo on next two lines
            if PipeDiafine[k] > 0:
                gamma_wall[k] = 2 * max(vr) * vth / PipeDiafine[k]
        
        # Estimate total reaction rate, including charge exchange and elastic scattering, and interaction with side walls 
        
        if mesh_type == 'h':
            minVr = vth * min(vr)
            minE0 = 0.5 * CONST.H_MASS * minVr * minVr / CONST.Q
            if CONST.USE_COLLRAD_IONIZATION:
                ioniz_rate = collrad_sigmav_ion_h0(nfine, Tefine)
            elif CONST.USE_JH:
                #Checks that the JH_Coef class has been passed into the function before calling jhs_coef
                if (jh_coeffs == None):
                    raise Exception("kinetic_h_mesh generated using JH, but no JH coefficients given")
                ioniz_rate = jhs_coef(nfine, Tefine, jh_coeffs, no_null = True) # deleted unecessary variable - GG
            else:
                ioniz_rate = sigmav_ion_h0(Tefine)
            RR = nfine * ioniz_rate + nfine * sigma_cx_h0(Tifine, np.array([minE0] * nxfine)) + gamma_wall # replaced size(nxfine) with nxfine
            
            # Compute local maximum grid spacing dx_max = 2 
            big_dx = 0.02 * fctr 
            dx_max = np.maximum(fctr * 0.8 * ( 2 * vth * min(vr) / RR), big_dx)

        elif mesh_type == 'h2':
            RR=nfine*sigmav_ion_hh(Tefine)+nfine*sigmav_h1s_h1s_hh(Tefine)+nfine*sigmav_h1s_h2s_hh(Tefine)+0.1*nfine*sigmav_cx_hh(Tifine,Tifine) + gamma_wall
             #Compute local maximum grid spacing from dx_max = 2 min(vr) / RR
            big_dx=0.02*fctr
            dx_max=np.minimum(fctr*0.8*(2*vth*min(vr)/RR), big_dx) # fixed typo - GG
            #TODO Check this calculation for dx_max here and above, they are the same in the idl code
            #TODO Note the usage of minimum vs maximum

        # NOTE See note in above statement, may be able to be combined here
        # # Compute local maximum grid spacing dx_max = 2 
        # big_dx = 0.02 * fctr 
        # dx_max = np.maximum(fctr * 0.8 * ( 2 * vth * min(vr) / RR), big_dx)

        # Construct xH Axis 
        xpt = xmaxH
        xH = np.array([xpt])
        while xpt > xminH:
            xH = np.concatenate([np.array([xpt]), xH]) # put xpt in array to fix concatenation error
            interpfunc = interpolate.interp1d(xfine, dx_max)
            dxpt1 = interpfunc(xpt)
            dxpt2 = dxpt1
            xpt_test = xpt - dxpt1
            if xpt_test > xminH:
                interpfunc = interpolate.interp1d(xfine, dx_max)
                dxpt2 = interpfunc(xpt_test)
            # dxpt = min([dxpt1, dxpt2]) * 0.5 ; FS: reduce spacing 
            # FS: force a preset maximum grid spacing 

            if mesh_type == 'h':
                dxh_max = 0.0004 # JWH: 0.0015 should be sufficient for D3D because scale lengths are 2.5x larger
                # lowered dxh_max from 5e-4 to 4e-4; original was giving mesh size errors in kinetic_h - nh
                dxpt = min([dxpt1, dxpt2, dxh_max])
            elif mesh_type == 'h2':
                dxpt=min([dxpt1,dxpt2])

            xpt = xpt - dxpt 
            
        # NOTE Old version, changed -2 to fit idl code xH = np.concatenate([np.array([xminH]), xH[0:np.size(xH) - 1]]) # put xminH in array to fix concatenation error
        xH = np.concatenate([np.array([xminH]), xH[0:np.size(xH) - 2]])
        # if xH[1] - xH[0] > 0.5 * big_dx:
        #   xH = np.concatenate(xH[0], xH[2:])

        interpfunc = interpolate.interp1d(xfine, Tifine)
        TiH = interpfunc(xH)
        interpfunc = interpolate.interp1d(xfine, Tefine)
        TeH = interpfunc(xH)
        interpfunc = interpolate.interp1d(xfine, nfine)
        neH = interpfunc(xH)
        interpfunc = interpolate.interp1d(xfine, PipeDiafine)
        PipeDiaH = interpfunc(xH)
        vx, vr, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nv, TiH) # fixed error from not assigning all outputs - GG

        self.x : NDArray = xH
        self.Ti : NDArray = TiH
        self.Te : NDArray = TeH
        self.ne : NDArray = neH
        self.PipeDia : NDArray = PipeDiaH
        self.vx : NDArray = vx
        self.vr : NDArray = vr
        self.Tnorm : float = Tnorm

        
def create_kinetic_h_mesh(
        nv         : int, 
        mu         : int, 
        x          : NDArray,
        Ti         : NDArray,
        Te         : NDArray, 
        n          : NDArray, 
        PipeDia    : NDArray,
        jh_coeffs  : JH_Coef = None,
        E0 = 0, ixE0 = 0, irE0 = 0, fctr = 1.0
    ) -> kinetic_mesh:
    
    mesh = kinetic_mesh('h', nv, mu, x, Ti, Te, n, PipeDia, jh_coeffs = jh_coeffs,
                        E0 = E0, ixE0 = ixE0, irE0 = irE0, fctr = fctr)
    return mesh

def create_kinetic_h2_mesh(
        nv         : int, 
        mu         : int, 
        x          : NDArray,
        Ti         : NDArray,
        Te         : NDArray, 
        n          : NDArray, 
        PipeDia    : NDArray,
        E0 = 0, ixE0 = 0, irE0 = 0, fctr = 1.0
    ) -> kinetic_mesh:
    
    mesh = kinetic_mesh('h2', nv, mu, x, Ti, Te, n, PipeDia, 
                        E0 = E0, ixE0 = ixE0, irE0 = irE0, fctr = fctr)
    return mesh