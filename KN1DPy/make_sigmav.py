import numpy as np

from .create_vr_vx_mesh import *
from .make_dvr_dvx import *

def make_sigmav(E_particle, mu_particle, T_target, mu_target, sigma_function):

    #   Input:
    #       E_Particle	    - array, eV - I changed the input type to use python vocab - GG
    #       mu_particle     - floatn
    #       T_Target	    - float, eV
    #       mu_target	    - float
    #       sigma_function  - function to return Sigma (m^2) for inputted E_Particle (eV); note this is no longer a string

  if type(E_particle)!=np.ndarray:
    E_particle=np.array(E_particle)
  
  trange=np.concatenate([E_particle,[T_target]])
  trange=np.sort(trange)
  nvb=100

  vxb, vrb, Tnorm, ixE0, irE0 = create_vr_vx_mesh(nvb,trange)
  Vr2pidVrb, VrVr4pidVrb, dVxb = make_dvr_dvx(vrb,vxb)[:3] 
  # fixed lines to properly use create_VrVxMesh and Make_dVr_dVx to set variables - nh

  mH=1.6726231e-27
  q=1.602177e-19
  vth=np.sqrt(2*q*Tnorm/(mu_target*mH))

    #   Set normalized particle velocities

  vxa=np.array([np.sqrt(2*q*i/(mu_particle*mH))/vth for i in E_particle]) 
    # original version was causing errors; this fixed it during testing, but it may need to be reviewed later - nh
  nvxa=vxa.size
  nvrb=vrb.size
  nvxb=vxb.size

    #   Set Normalized Target Distribution Function

  fi_hat=np.zeros((nvxb,nvrb))
  for i in range(nvrb):
    fi_hat[:,i]=np.e**(-(vrb[i]**2+vxb**2)*Tnorm/T_target)
  fi_hat=fi_hat/(Vr2pidVrb*np.matmul(dVxb,fi_hat)).sum() # fixed indexing errors and order of matmul arguments - nh

    #   Compute relative velocity at each mesh point

  vrel=np.zeros((nvxa,nvxb,nvrb))
  for k in range(nvxa):
    for i in range(nvrb):
      vrel[k,:,i]=np.sqrt(vrb[i]**2+(vxb-vxa[k])**2)

    #   Get sigma for inputted E_Particle
    
  sig=sigma_function(0.5*vth*vth*vrel*vrel*mu_particle*mH/q) # changed to no longer use exec()

    #   Compute Sigmav by integrating sigma x Vrel x Fi_hat over velocity space

  sigv=np.zeros(nvxa)  
  for k in range(nvxa):
    sigv[k]=vth*(Vr2pidVrb*np.matmul(dVxb,sig[k,:,:]*vrel[k,:,:]*fi_hat)).sum() # fixed indexing errors and order of matmul arguments - nh

  return sigv
