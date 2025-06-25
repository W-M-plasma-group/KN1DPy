import numpy as np

from .make_dvr_dvx import make_dvr_dvx
from .create_shifted_maxwellian_include import create_shifted_maxwellian_include
from .common import constants as CONST

def create_shifted_maxwellian(vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm): # fixed function name - nh
  nx=vx_shift.size
  nvr,nvx=vr.size,vx.size
  maxwell=np.zeros((nvr,nvx,nx)).T
  vr2vx2_ran2=np.zeros((nvr,nvx)).T

  vth=np.sqrt(2*CONST.Q*Tnorm/(mu*CONST.H_MASS))
  vth2=vth**2
  vth3=vth**3

  shifted_maxwellian_debug=0

  Vr2pidVr,VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb = make_dvr_dvx(vr,vx)

  return create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_maxwellian_debug,mu,mol,nx,nvx,nvr,vth,vth2,
                                           maxwell,vr2vx2_ran2,Vr2pidVr,dVx,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb)