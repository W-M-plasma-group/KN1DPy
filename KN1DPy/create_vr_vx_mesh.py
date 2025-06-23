#from cmath import sqrt     replaced with np.sqrt - nh
import numpy as np
from numpy.typing import NDArray

from .reverse import * # fixed import

# sets up optimum Vr and Vx velocity space mesh for Kinetic_Neutrals procedure 
# Input: 
#   nv - Integer, number of elements desired in vr mesh
#   Ti - arrray, Ti profile
#   E0 - array, energy where a velocity is desired ( optional )
#   Tmax - float, ignore Ti above this value
#
# Gwendolyn Galleher 

def create_vr_vx_mesh(nv : int, Ti : NDArray, E0 : NDArray = np.array([0.0]), Tmax : float = 0.0
                      ) -> tuple[NDArray, NDArray, float, NDArray, NDArray] : # removed unused input parameters
    _Ti = np.array(Ti) 
    _Ti = np.concatenate([_Ti, E0[E0>0]]) 
    if Tmax > 0:
        _Ti = _Ti[_Ti<Tmax] # simplified previous lines by replacing loops - nh
    maxTi = np.max(_Ti)
    minTi = min(_Ti)
    Tnorm = np.mean(_Ti)
    vmax = 3.5
    if maxTi - minTi <= 0.1 * maxTi: # changed < to <= like from IDL code - nh
        v = np.arange(nv+1)*vmax/nv # added *vmax/nv from IDL code - nh // deleted float type because it caused bug - GG
    else:
        G = 2 * nv * np.sqrt(minTi / maxTi) / (1- np.sqrt(minTi / maxTi))
        b = vmax / (nv * (nv + G))
        a = G * b 
        v = a * np.arange(nv +1) + b * (np.arange(nv+1) ** 2)

        # Option: add velocity bins corresponding to E0 
        
    v0=0 # this line was missing - nh
    for k in range(np.size(E0)): # fixed range - nh
        if E0[k] > 0.0:
            v0 = np.sqrt(E0[k]/Tnorm)
            ii = np.argwhere(v > v0).T[0] # argwhere has a weird output, but this should put it in a usable format - nh
            if np.size(ii) > 0: # removed count variable
                v = np.concatenate([v[0 : ii[0]], [v0], v[ii[0] : ]])
            else: 
                v = np.concatenate([v, v0]) # changed lines to np.concatenate - nh
        
    vr = v[1 : ] # fixed indexing - nh
    vx = np.concatenate([-reverse(vr), vr]) 

    ixE0 = np.argwhere(abs(vx) == v0).T[0] # modified argwhere output - nh
    if np.size(ixE0) == 1: # removed count variable; same changes made to following lines - nh
        ixE0 = ixE0[0]
    irE0 = np.argwhere(vr == v0).T[0]
    if np.size(irE0) == 1:
        irE0 = irE0[0]
    return vx,vr,Tnorm,ixE0,irE0 # changed to return outputs as a variables - GG
