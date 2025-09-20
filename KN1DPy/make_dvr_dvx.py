# Make_dVr_dVx
#   Constructs velocity space differentials for distrobution functions 
# used by Kinetic_Neutrals, Kinetic_H2, Kinetic_H, and other related
# procedures 
#
# Gwendolyn Galleher 

import numpy as np
from numpy.typing import NDArray
import copy

''' 
Authors: Julio Balbin, Carlo Becerra
Date: August 17th, 2024
'''

def make_dvr_dvx(vr : NDArray, vx : NDArray):

    """
    Constructs velocity space differentials for distribution functions.

    Parameters:
    -----------
    vr : np.ndarray
        Array of radial velocities.
    vx : np.ndarray
        Array of axial velocities.

    Returns:
    --------
    Vr2pidVr : np.ndarray
        Differential volume element for radial velocities.
    VrVr4pidVr : np.ndarray
        Differential volume element for radial velocities (higher order).
    dVx : np.ndarray
        Differential for axial velocities.
    vrL, vrR : np.ndarray
        Left and right boundaries for radial velocities.
    vxL, vxR : np.ndarray
        Left and right boundaries for axial velocities.
    vol : np.ndarray
        Volume elements in velocity space.
    vth_Deltavx, vx_Deltavx, vr_Deltavr : np.ndarray
        Auxiliary quantities for kinetic equations.
    vr2vx2 : np.ndarray
        Squared magnitude of the velocity.
    jpa, jpb, jna, jnb : int
        Indices for positive and negative axial velocities.
    """

    # nvr & nvx are taking the shape of vr & vx respectively
    nvr = vr.size
    nvx = vx.size

    # Calculations for r-dimension
    _vr    = np.append(vr, 2 * vr[-1] - vr[-2])
    vr_mid = np.concatenate(([0.0], 0.5 * (_vr + np.roll(_vr, -1))))

    vrR = np.roll(vr_mid, -1)[0:nvr]
    vrL = copy.copy(vr_mid)[0:nvr]

    Vr2pidVr    =         np.pi * (vrR**2 - vrL**2)
    VrVr4pidVr  = (4/3) * np.pi * (vrR**3 - vrL**3)

    # Calculations for x-dimension
    _vx = np.concatenate(([2 * vx[0] - vx[1]], vx, [2 * vx[-1] - vx[-2]]))

    vxR = 0.5 * (np.roll(_vx, -1) + _vx)[1:nvx+1]
    vxL = 0.5 * (np.roll(_vx,  1) + _vx)[1:nvx+1]

    dVx = vxR - vxL

    # Calc. volume
    # --------------------------------------------------------------------
    vol = np.zeros((nvr, nvx), dtype=np.float64)
    # # for i in tqdm(range(nvr), desc=f"Make_dVr_dVx: calc. vol"): 
    # for i in range(nvr):     
    #     vol[i, :] = Vr2pidVr[i] * dVx
    vol = Vr2pidVr[:, np.newaxis] * dVx
    # --------------------------------------------------------------------
    Deltavx = vxR - vxL
    Deltavr = vrR - vrL
    # --------------------------------------------------------------------
    vth_Deltavx = np.zeros((nvr+2,nvx+2))
    vx_Deltavx  = np.zeros((nvr+2,nvx+2))
    vr_Deltavr  = np.zeros((nvr+2,nvx+2))
    # for j in tqdm(range(1,nvr+1),desc=f"Make_dVr_dVx: calc. vth & vx"):
    # for j in range(1,nvr+1):
        # vth_Deltavx[j,1:nvx+1] = 1.0/Deltavx     # vth_Deltavx(i,1:nvx)=1.0/Deltavx
        # vx_Deltavx[ j,1:nvx+1] =  vx/Deltavx     #  vx_Deltavx(i,1:nvx)= vx/Deltavx
    # for k in tqdm(range(1,nvx+1),desc=f"Make_dVr_dVx: calc. vr"):
    # for k in range(1,nvx+1):
    #     vr_Deltavr[ 1:nvr+1,k] =  vr/Deltavr     #  vr_Deltavr(1:nvr,j)=vr/Deltavr
    vth_Deltavx[1:nvr+1, 1:nvx+1] = 1.0 / Deltavx
    vx_Deltavx[ 1:nvr+1, 1:nvx+1] = vx  / Deltavx 
    vr_Deltavr[ 1:nvr+1, 1:nvx+1] = vr[:, np.newaxis] / Deltavr[:, np.newaxis]
    # --------------------------------------------------------------------

    # Compute v^2
    vr2vx2=np.zeros((nvr,nvx), dtype=np.float64)
    # for l in range(0,nvr):
    #     vr2vx2[l,:] = vr[l]**2 + vx**2
    vr2vx2 = vr[:, np.newaxis]**2 + vx**2 

    # vx's positive index
    jp = np.where(vx>0)[0]   # This saves the positives index of vx
    jpa = int(jp[0])                    # This saves the first index of jp
    jpb = int(jp[len(jp)-1])            # This saves the last index of jp

    # vx's negative index   
    jn = np.where(vx<0)[0]   # This saves the negatives index of vx
    jna = int(jn[0])                    # This saves the first index of jp
    jnb = int(jn[len(jn)-1])            # This saves the last index of jp

    # # In case there are only positive or negative values, we can use this option
    jpa = int(jp[ 0]) if len(jp) > 0 else None
    jpb = int(jp[-1]) if len(jp) > 0 else None
    jna = int(jn[ 0]) if len(jn) > 0 else None
    jnb = int(jn[-1]) if len(jn) > 0 else None

    return Vr2pidVr,VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,\
           vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,\
           jpa,jpb,jna,jnb