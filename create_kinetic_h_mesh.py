import numpy as np 
from scipy import interpolate
from sigmav_ion_h0 import sigmav_ion_h0
from create_vrvxmesh import create_VrVxMesh
from collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from sigma_cx_h0 import sigma_cx_h0
from JHS_coef import JHS_coef

def create_kinetic_h_mesh(nv, mu, x, Ti, Te, n, PipeDia, E0 = 0, ixE0 = 0 ,irE0 = 0,fctr = 1, g=None): # added common block input
    # add these lines to improve velocity space resolution (use with care):
    nv = 20
    # E0[.1]
    JH = 0 
    Use_Collrad_Ionization = 1
    mH = 1.6726231e-27
    q=1.602177e-19
    nx = np.size(x)

    # estimate Interaction rate with side walls
    gamma_wall = np.zeros(np.shape(nv))
    # for k in range( 0, nx):
    #   if PipeDia[k] > 0:
    #        gamma_wall[k] = 2*sqrt(2*Ti(k)*q/(2*mH))/PipeDia(k)

    # Estimate total reaction rate for destriction of hydrogen atoms and for interation with side walls
    # RR = n*sigmav_ion_H0(Te)+gamma_wall
    RR = n * sigmav_ion_h0(Te) 

    # Set v0 to thermal speed to 10 eV neutral 
    v0 = np.sqrt( 2 * 10 * q / (mu * mH))

    # Determine x range for atoms by finding distance into plasma where density persists.
    #  dGamma/dx=-nH*RR = v0 dnH/dx = -nH*RR - these two lines are commented in the original code I dont understand them
    #  d ln(nH)/dx = -RR/v0 = dY/dx

    Y = np.zeros(nx) # fixed typo
    for k in range(1, nx): 
        Y[k] = Y[k-1] - (x[k] - x[k-1] ) * 0.5 * (RR[k] + RR[k-1]) / v0
    
    # Find x location where Y = -5, i.e. where nH should be down by exp(-5)
    interpfunc = interpolate.interp1d(Y, x, kind = 'linear', bounds_error=False, fill_value="extrapolate") # previous version caused error on next line
    xmaxH = np.minimum(interpfunc(-5), max(x))
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
    vx, vr, Tnorm, ixE0, ixE0 = create_VrVxMesh(nv, Tifine)[:3] # fixed error from not assigning all outputs - GG
    vth = np.sqrt( 2 * q * Tnorm / (mu * mH))
    minVr = vth * min(vr)
    minE0 = 0.5 * mH * minVr * minVr / q

    # Estimate interaction rate with side walls
    nxfine = np.size(xfine)
    gamma_wall = np.zeros(nxfine)
    #print(PipeDiafine)
    for k in range(0, nxfine): # fixed typo on next two lines
        if PipeDiafine[k] > 0:
            gamma_wall[k] = 2 * max(vr) * vth / PipeDiafine[k]
    
    # Estimate total reaction rate, including charge exchange and elastic scattering, and interaction with side walls 
    if Use_Collrad_Ionization:
        ioniz_rate = collrad_sigmav_ion_h0(nfine, Tefine)
    else:
        if JH:
            ioniz_rate = JHS_coef(nfine, Tefine, no_null = True) # deleted unecessary variable - GG
        else:
            ioniz_rate = sigmav_ion_h0(Tefine)
    RR = nfine * ioniz_rate + nfine * sigma_cx_h0(Tifine, np.array([minE0] * nxfine)) + gamma_wall # replaced size(nxfine) with nxfine

    # Compute local maximum grid spacing dx_max = 2 
    big_dx = 0.02 * fctr 
    dx_max = np.maximum(fctr * 0.8 * ( 2 * vth * min(vr) / RR), big_dx)

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
        dxh_max = 0.0005 # JWH: 0.0015 should be sufficient for D3D because scale lengths are 2.5x larger
        dxpt = min([dxpt1, dxpt2, dxh_max])
        xpt = xpt - dxpt 
    xH = np.concatenate([np.array([xminH]), xH[0:np.size(xH) - 2]]) # put xminH in array to fix concatenation error
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
    vx, vr, Tnorm, ixE0, ixE0 = create_VrVxMesh(nv, TiH)[:3] # fixed error from not assigning all outputs - GG
    return xH, TiH, TeH, neH, PipeDiaH, vx, vr, Tnorm
