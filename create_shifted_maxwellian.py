import numpy as np
from sval import sval
from Make_dVr_dVx import Make_dVr_dVx
from create_shifted_maxwellian_include import create_shifted_maxwellian_include
from global_vars import mH, q

def create_shifted_maxwellian(vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm,numero=0): # fixed function name - nh
  variables = {
    'vr': vr,
    'vx': vx, 
    'Tmaxwell': Tmaxwell, 
    'vx_shift': vx_shift, 
    'mu': mu, 
    'mol': mol, 
    'Tnorm': Tnorm
    }
  # Abrir el archivo en modo de escritura
  with open(f'variables_created_shifted_maxwellian_{numero}.txt', 'w') as file:
      for name, value in variables.items():
          if hasattr(value, '__len__'):
              file.write(f'{name}: {value}, tipo: {type(value)}, len({name}): {value.shape}\n')
          else:
              file.write(f'{name}: {value}, tipo: {type(value)}\n')
  nx=vx_shift.size
  nvr,nvx=vr.size,vx.size
  maxwell=np.zeros((nvr,nvx,nx)).T
  vr2vx2_ran2=np.zeros((nvr,nvx)).T

  vth=np.sqrt(2*q*Tnorm/(mu*mH))
  vth2=vth**2
  vth3=vth**3

  shifted_maxwellian_debug=0

  Vr2pidVr,VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb = Make_dVr_dVx(vr,vx)

  return create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_maxwellian_debug,mu,mol,nx,nvx,nvr,vth,vth2,
                                           maxwell,vr2vx2_ran2,Vr2pidVr,dVx,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb)