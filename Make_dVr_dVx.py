# Make_dVr_dVx
#   Constructs velocity space differentials for distrobution functions 
# used by Kinetic_Neutrals, Kinetic_H2, Kinetic_H, and other related
# procedures 
#
# Gwendolyn Galleher 

import numpy as np

def Make_dVr_dVx(Vr, Vx): # For this to work inputs must be np arrays so I might need to ammend this later to make sure all inputs are arrays

    Vr=np.array(Vr)
    Vx=np.array(Vx) # sets inputs to np.array if they aren't already; resolves earlier comment - nh
    
    # Determine velocity space differentials 
    nVr = np.size(Vr)
    nVx = np.size(Vx)

    # for Vr first 
    _Vr = np.concatenate([Vr, [2 * Vr[nVr-1] - Vr[nVr-2]]])  # this is an array and Vr(nVr-1) is calling the last cell of Vr
    Vr_mid = np.concatenate([[0.0], 0.5 * (_Vr + np.roll(_Vr, -1))]) # changed to np.concatenate - nh
    
    VrR = np.roll(Vr_mid, -1)
    VrL = Vr_mid

    Vr2pidVr = np.pi * ((VrR ** 2) - (VrL ** 2))
    Vr2pidVr = Vr2pidVr[0 : nVr]  # makes it the same length as Vr 

    VrVr4pidVr = (4/3) * np.pi * ((VrR ** 3) - (VrL ** 3))
    VrVr4pidVr = VrVr4pidVr[0 : nVr]
    VrR = VrR[0 : nVr]
    VrL = VrL[0 : nVr] 

    # now for Vx
    _Vx = np.concatenate([[2 * Vx[0] - Vx[1]], Vx, [2 * Vx[nVx - 1] - Vx[nVx - 2]]]) # changed to np.concatenate - nh
    VxR = 0.5 * (np.roll(_Vx, -1) + _Vx)
    VxL = 0.5 * (np.roll(_Vx, 1) + _Vx)
    dVx = VxR[1: nVx+1] - VxL[1:nVx+1]
    VxR = VxR[1: nVx+1]
    VxL = VxL[1 : nVx+1]

    # compute volume elements 
    vol = np.zeros((nVx, nVr), float)
    for i in range(0, nVr):
        vol[:,i] = Vr2pidVr[i] * dVx # fixed minor indexing bugs - nh

    #compute DeltaVx, DeltaVr
    DeltaVx = VxR - VxL
    DeltaVr = VrR - VrL

    # compute vth_Deltavx, vx_Deltavx, vr_Deltavr, padded with zeros
    Vth_DeltaVx = np.zeros((nVx + 2, nVr + 2), float)
    Vx_DeltaVx = np.zeros((nVx + 2, nVr + 2), float)
    Vr_DeltaVr = np.zeros((nVx + 2, nVr + 2), float) # replaced np.array with np.zeros - nh
    for i in range(1, nVr+1):
        Vth_DeltaVx[1 : nVx+1,i] = 1.0/DeltaVx
        Vx_DeltaVx[1 : nVx+1,i] = Vx/DeltaVx
    for j in range(1, nVx+1):
        Vr_DeltaVr[j,1 : nVr+1] = Vr/DeltaVr 
    
    #compute v^2
    Vr2Vx2 = np.zeros((nVx, nVr), float)
    for i in range(0, nVr):
        Vr2Vx2[:,i] = (Vr[i] ** 2) + (Vx ** 2)

    # Determine indice range of positive and negative Vx 
    jpa=jpb=jna=jnb=-1
    jp = np.argwhere(Vx > 0)
    if jp.size>0:
        jpa = jp[0][0]; jpb = jp[np.size(jp) - 1][0]
    jn = np.argwhere(Vx < 0)
    if jn.size>0:
        jna = jn[0][0]; jnb = jn[np.size(jn) - 1][0] # modified section to return -1 if jp or jn is empty (previously this raised an error) - nh
    
    # changed return line to provide an output as a list - nh
    return [Vr2pidVr,VrVr4pidVr,dVx,VrL,VrR,VxL,VxR,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2,jpa,jpb,jna,jnb]
