import numpy as np
from create_vrvxmesh import *
from Make_dVr_dVx import *

def make_sigmav(E_particle, mu_particle, T_target, mu_target, sigma_function):

    #   Input:
    #       E_Particle	    - fltarr(*), eV
    #       mu_particle     - float
    #       T_Target	    - float, eV
    #       mu_target	    - float
    #       sigma_function  - string, name of function to return Sigma (m^2)
    #                           for inputted E_Particle (eV)

  if type(E_particle)!=np.ndarray:
    E_particle=np.array(E_particle)
  
  trange=np.concatenate([E_particle,[T_target]])
  trange=np.sort(trange)
  nvb=100

  vxb, vrb, Tnorm, ixE0, irE0 = create_VrVxMesh(nvb,trange)
  Vr2pidVrb, VrVr4pidVrb, dVxb = Make_dVr_dVx(vrb,vxb)[:3] 
  # fixed lines to properly use create_VrVxMesh and Make_dVr_dVx to set variables - nh

  mH=1.6726231e-27
  q=1.602177e-19
  vth=np.sqrt(2*q*Tnorm/(mu_target*mH))

    #   Set normalized particle velocities

  vxa=np.sqrt(2*q*E_particle/(mu_particle*mH))/vth
  nvxa=vxa.size
  nvrb=vrb.size
  nvxb=vxb.size

    #   Set Normalized Target Distribution Function

  fi_hat=np.zeros((nvrb,nvxb))
  for i in range(nvrb):
    fi_hat[i,:]=np.e**(-(vrb[i]**2+vxb**2)*Tnorm/T_target)
  fi_hat=fi_hat/(Vr2pidVrb*np.matmul(fi_hat,dVxb)).sum()

    #   Compute relative velocity at each mesh point

  vrel=np.zeros((nvrb,nvxb,nvxa))
  for k in range(nvxa):
    for i in range(nvrb):
      vrel[i,:,k]=np.sqrt(vrb[i]**2+(vxb-vxa[k])**2)

    #   Get sigma for inputted E_Particle
    
  exec('from '+sigma_function.lower()+' import *; global sig; sig='+sigma_function+'(0.5*vth*vth*vrel*vrel*mu_particle*mH/q)')

    #   Compute Sigmav by integrating sigma x Vrel x Fi_hat over velocity space

  sigv=np.zeros(nvxa)  
  for k in range(nvxa):
    sigv[k]=vth*(Vr2pidVrb*np.matmul(sig[:,:,k]*vrel[:,:,k]*fi_hat,dVxb)).sum()

  return sigv
