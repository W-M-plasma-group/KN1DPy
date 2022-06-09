from cmath import sqrt
from reverse.py import *
import numpy as np

# sets up optimum Vr and Vx velocity space mesh for Kinetic_Neutrals procedure 
# Input: 
#   nv - Integer, number of elements desired in vr mesh
#   Ti - arrray, Ti profile
#   E0 - array, energy where a velocity is desired ( optional )
#   Tmax - float, ignore Ti above this value
#
# Gwendolyn Galleher 

def create_VrVxMesh(nv, Ti, vx, Tnorm, ixE0, irE0, E0 = np.array([0.0]), Tmax = 0.0):
    _Ti = Ti 
    for k in range(0, np.size(E0) - 1):
        if E0[k] > 0:
            _Ti = np.array([_Ti, E0[k]])
    if Tmax > 0:
        ii = np.argwhere(_Ti < Tmax)
        count = np.size(ii)
        if count > 0:
            _Ti = _Ti[ii]
    maxTi = np.max(_Ti)
    minTi = min(_Ti)
    Tnorm = np.mean(_Ti)
    vmax = 3.5
    if maxTi - minTi < 0.1 * maxTi:
        v = np.arange(nv+1, float)
    else:
        G = 2 * nv * sqrt(minTi / maxTi) / (1- sqrt(minTi / maxTi))
        b = vmax / (nv * (nv + G))
        a = G * b 
        v = a * np.arange(nv +1) + b * (np.arange(nv+1) ** 2)

        # Option: add velocity bins corresponding to E0 

        for k in range(np.size(E0) - 1):
            if E0[k] > 0.0:
                v0 = sqrt(E0[k]/Tnorm)
                ii = np.argwhere(v > v0)
                count = np.size(ii)
                if count > 0:
                    v = np.array([v[0 : ii[0]], v0, v[ii[0] : v[np.size(v) - 1]]])
                else: 
                    v = [v, v0]
        
        vr = v[1 : np.size(v) - 1]
        vx = [-reverse(vr), vr]

        ixE0 = np.argwhere(abs(vx) == v0)
        count = np.size(ixE0)
        if count == 1:
            ixE0 = ixE0[0]
        irE0 = np.argwhere(vr == v0)
        count = np.size(irE0)
        if count == 1:
            irE0 = irE0[0]
        return 