import numpy as np
from sigmav_ion_hh import sigmav_ion_hh
from sigmav_h1s_h1s_hh import sigmav_h1s_h1s_hh
from sigmav_h1s_h2s_hh import sigmav_h1s_h2s_hh
from sigmav_cx_hh import sigmav_cx_hh
from create_vrvxmesh import create_vrvxmesh
from scipy import interpolate

def create_kinetic_h2_mesh(nv, mu, x, Ti, Te, n, PipeDia, E0 = 0, ixE0 = 0, irE0 = 0, fctr = 1.0): # - removed output variables from input - GG
    mH = 1.6726231e-27		
    q = 1.602177e-19				
    k_boltz = 1.380658e-23				#Boltzmann's constant, J K^-1
    Twall = 293.0*k_boltz/q			    #room temperature (eV)
    v0_bar = np.sqrt(8.0*Twall*q/(np.pi*2*mu*mH))	#directed random velocity of diatomic molecule
    nx=np.size(x)

    gamma_wall = [0] * nx

    #Estimate total reaction rate for destruction of molecules and for interation with side walls
    RR=n*sigmav_ion_hh(Te)+n*sigmav_h1s_h1s_hh(Te)+n*sigmav_h1s_h2s_hh(Te)

    Y = [0] * nx
    for k in range(1, nx-1):
        Y[k]=Y[k-1]-(x[k]-x[k-1])*0.5*(RR[k]+RR[k-1])/v0_bar

    #Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
    interpfunc = interpolate.interp1d(Y, x) # fixed error with interpolation - GG
    xmaxH2=np.minimum(interpfunc(-10.0), max(x))

    xminH2=x[0]

    #Interpolate Ti and Te onto a fine mesh between xminH2 and xmaxH2
    xfine=xminH2+(xmaxH2-xminH2)*np.arange(1001)/1000 # fixed error with findgen - GG

    interpfunc = interpolate.interp1d(x, Ti, kind = 'linear') # fixed errors with interpolation
    Tifine=interpfunc(xfine)
    
    interpfunc = interpolate.interp1d(x, Te, kind = 'linear')
    Tefine=interpfunc(xfine)
    
    interpfunc = interpolate.interp1d(x, n, kind = 'linear')
    nfine=interpfunc(xfine)
    
    interpfunc = interpolate.interp1d(x, PipeDia)
    PipeDiafine=interpfunc(xfine)

    #Setup a vx,vr mesh based on raw data to get typical vx, vr values
    #probably need to do stuff about the namespace because of the differences between IDL and python
    vrvxmesh = create_vrvxmesh(nv, Tifine) # pulled necessary variables from the return of create_vrvxmesh - GG
    vx = vrvxmesh[0]
    vr = vrvxmesh[1]
    Tnorm = vrvxmesh[2]

    vth = np.sqrt(2*q*Tnorm/(mu*mH))

    #Estimate interaction rate with side walls
    nxfine=np.size(xfine)
    gamma_wall = [0] * nxfine
    for k in range(nxfine-1):
        if PipeDiafine(k) > 0:
            gamma_wall[k]=2*max(vr)*vth/PipeDiafine[k]

    #Estimate total reaction rate, including charge exchange, elastic scattering, and interaction with side walls
    RR=nfine*sigmav_ion_hh(Tefine)+nfine*sigmav_h1s_h1s_hh(Tefine)+nfine*sigmav_h1s_h2s_hh(Tefine)+0.1*nfine*sigmav_cx_hh(Tifine,Tifine) + gamma_wall

    #Compute local maximum grid spacing from dx_max = 2 min(vr) / RR
    big_dx=0.02*np.fctr
    dx_max=np.minimum(np.fctr*0.8*(2*vth*min(vr)/RR), big_dx)

    #Construct xH2 axis
    xpt=xmaxH2
    xH2=np.array([xpt])
    while xpt > xminH2:
        xH2=np.concatenate([xpt,xH2])
        interpfunc = interpolate.interp1d(xfine, dx_max)
        dxpt1=interpfunc(xpt)
        dxpt2=dxpt1
        xpt_test=xpt-dxpt1
        if xpt_test > xminH2:
            interpfunc = interpolate.interp1d(xfine, dx_max)
            dxpt2=interpfunc(xpt_test)
        dxpt=min([dxpt1,dxpt2])
        xpt=xpt-dxpt
    xH2=np.concatenate(np.array([xminH2]), xH2[0:np.size(xH2) - 2]) 

    interpfunc = interpolate.interp1d(xfine, Tifine)
    print(xH2)
    TiH2=interpfunc(xH2)

    interpfunc = interpolate.interp1d(xfine, Tefine)
    TeH2=interpfunc(xH2)

    interpfunc = interpolate.interp1d(xfine, nfine)
    neH2=interpfunc(xH2)

    interpfunc = interpolate.interp1d(xfine, PipeDiafine)
    PipeDiaH2=interpfunc(xH2)

    vrvxmesh = create_VrVxMesh(nv,TiH2)
    vx = vrvxmesh[0]
    vr = vrvxmesh[1]
    Tnorm = vrvxmesh[2]

    return xH2, TiH2, TeH2, neH2, PipeDiaH2, vx, vr, Tnorm # returned necessary varibales in a list - GG