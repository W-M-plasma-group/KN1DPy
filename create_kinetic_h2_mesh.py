import numpy as np
from sigmav_ion_hh import sigmav_ion_hh
from sigmav_h1s_h1s_hh import sigmav_h1s_h1s_hh
from sigmav_h1s_h2s_hh import sigmav_h1s_h2s_hh
from sigmav_cx_hh import sigmav_cx_hh
from create_vrvxmesh import create_vrvxmesh
from scipy import interpolate

def create_kinetic_h2_mesh(nv, mu, x, Ti, Te, n, PipeDia, E0=0,ixE0=0,irE0=0,fctr=0): # - removed output variables from input - GG
    mH = 1.6726231e-27		
    q = 1.602177e-19				
    k_boltz = 1.380658e-23				#Boltzmann's constant, J K^-1
    Twall = 293.0*k_boltz/q			#room temperature (eV)
    v0_bar = np.sqrt(8.0*Twall*q/(np.pi*2*mu*mH))	#directed random velocity of diatomic molecule
    nx=len(x)

    gamma_wall = [0] * nx

    #Estimate total reaction rate for destruction of molecules and for interation with side walls
    RR=n*sigmav_ion_hh(Te)+n*sigmav_h1s_h1s_hh(Te)+n*sigmav_h1s_h2s_hh(Te)

    Y = [0] * nx
    for k in range(1, nx-1):
        Y[k]=Y[k-1]-(x[k]-x[k-1])*0.5*(RR[k]+RR[k-1])/v0_bar

    #Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
    xmaxH2=np.minimum(interpolate.interp1d(x,Y,-10.0), max(x))
    xminH2=x[0]

    #Interpolate Ti and Te onto a fine mesh between xminH2 and xmaxH2
    findgen = [None] * 1000
    for k in range(1001):
        findgen[k] = float(k)
    xfine=xminH2+(xmaxH2-xminH2)*findgen/1000
    Tifine=interpolate.interp1d(Ti,x,xfine)
    Tefine=interpolate.interp1d(Te,x,xfine)
    nfine=interpolate.interp1d(n,x,xfine)
    PipeDiafine=interpolate.interp1d(PipeDia,x,xfine)

    #Setup a vx,vr mesh based on raw data to get typical vx, vr values
    #probably need to do stuff about the namespace because of the differences between IDL and python
    vrvxmesh = create_vrvxmesh(nv, Tifine) # pulled necessary variables from the return of create_vrvxmesh - GG
    vx = vrvxmesh[0]
    vr = vrvxmesh[1]
    Tnorm = vrvxmesh[2]
    vth = np.sqrt(2*q*Tnorm/(mu*mH))

    #Estimate interaction rate with side walls
    nxfine=len(xfine)
    gamma_wall = [0] * nxfine
    for k in range(nxfine-1):
        if PipeDiaFine(k) > 0:
            gamma_wall[k]=2*max(vr)*vth/PipeDiaFine[k]

    #Estimate total reaction rate, including charge exchange, elastic scattering, and interaction with side walls
    RR=nfine*sigmav_ion_hh(Tefine)+nfine*sigmav_h1s_h1s_hh(Tefine)+nfine*sigmav_h1s_h2s_hh(Tefine)+0.1*nfine*sigmav_cx_hh(Tifine,Tifine) + gamma_wall

    #Compute local maximum grid spacing from dx_max = 2 min(vr) / RR
    big_dx=0.02*np.fctr
    dx_max=np.minimum(np.fctr*0.8*(2*vth*min(vr)/RR), big_dx)

    #Construct xH2 axis
    xpt=xmaxH2
    xH2=[xpt]
    while xpt > xminH2:
        xH2=[xpt,xH2]
        dxpt1=interpolate.interp1d(dx_max,xfine,xpt)
        dxpt2=dxpt1
        xpt_test=xpt-dxpt1
        if xpt_test > xminH2:
            dxpt2=interpolate.interp1d(dx_max,xfine,xpt_test)
        dxpt=min([dxpt1,dxpt2])
        xpt=xpt-dxpt
    xH2=[xminH2,xH2[0:len(xH2)-2]]

    TiH2=interpolate.interp1d(Tifine,xfine,xH2)
    TeH2=interpolate.interp1d(Tefine,xfine,xH2)
    neH2=interpolate.interp1d(nfine,xfine,xH2)
    PipeDiaH2=interpolate.interp1d(PipeDiafine,xfine,xH2)
    vrvxmesh = create_vrvxmesh(nv,TiH2) # pulled necessary variables from the return of create_vrvxmesh - GG
    vx = vrvxmesh[0]
    vr = vrvxmesh[1]
    Tnorm = vrvxmesh[2]
    return [xH2, TiH2, TeH2, neH2, PipeDiaH2, vx, vr, Tnorm] # returned necessary varibales in a list - GG