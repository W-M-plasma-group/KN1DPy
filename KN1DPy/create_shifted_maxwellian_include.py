import numpy as np

from .sval import sval
from .common import constants as CONST
from .make_dvr_dvx import make_dvr_dvx
import copy
from tqdm import tqdm

#   This INCLUDE file is used by Kinetic_H2 and Kinetic_H
#   The code is also written within the create_shifted_maxwellian function

def create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                      nx,nvx,nvr,vth,vth2,maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb): # added vr,vx,Tnorm,mu to arguments; should not have been left out originally

    #   Input:
	#       Vx_shift  - dblarr(nx), (m s^-1)
    #       Tmaxwell  - dblarr(nx), (eV)
    #       Shifted_Maxwellian_Debug - if set, then print debugging information
    #       mol       - 1=atom, 2=diatomic molecule
 
    #   Output:
	#       Maxwell   - dblarr(nvr,nvx,nx) a shifted Maxwellian distribution function
	#	        having numerically evaluated vx moment close to Vx_shift and
	#	        temperature close to Tmaxwell

    #   Notes on Algorithm:

    #   One might think that Maxwell could be simply computed by a direct evaluation of the EXP function:

    #       for i=0,nvr-1 do begin
    #           arg=-(vr(i)^2+(vx-Vx_shift/vth)^2) * mol*Tnorm/Tmaxwell
    #           Maxwell(i,*,k)=exp(arg > (-80))
    #       endfor

    #   But owing to the discrete velocity space bins, this method does not necessarily lead to a digital representation 
    #   of a shifted Maxwellian (Maxwell) that when integrated numerically has the desired vx moment of Vx_shift
    #   and temperature, Tmaxwell.

    #   In order to insure that Maxwell has the desired vx and T moments when evaluated numerically, a compensation
    #   scheme is employed - similar to that used in Interp_fVrVxX
  

  Vr2pidVr, VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb = make_dvr_dvx(vr,vx)

  AN=np.zeros((nvr,nvx,2), float)
  BN=np.zeros((nvr,nvx,2), float) # fixed array creation - nh
  sgn=[1,-1]
  maxwell = np.zeros((nvr,nvx,nx),dtype=np.float32)
  # print("nx", nx)
  # print("Tmaxwell", Tmaxwell)
  # input()
  for k in tqdm(range(nx),desc=f'k'):
    if Tmaxwell[k] > 0:
      arg = -((vr[:, np.newaxis]**2 + (vx - (vx_shift[k] / vth))**2) * mol * Tnorm / Tmaxwell[k])
      # print("arg", arg)
      # input()
      arg = np.where(np.logical_and((-80 < arg), (arg < 0.0)), arg, -80)
      maxwell[:, :, k] = np.exp(arg)
      # input()

      # print("Vr2pidVr", Vr2pidVr)
      # print("dVx", dVx)
      # print("maxwell", maxwell.T)
      # print("maxwell#dvx", np.matmul(maxwell[:,:,k],dVx))
      # input()
      variable = np.matmul(maxwell[:, :, k], dVx)
      maxwell[:, :, k] = maxwell[:, :, k] / np.nansum(Vr2pidVr * variable)
      
      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out1=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out1/vth)**2
        T_out1=(mol*mu*CONST.H_MASS)*vth2*np.sum(Vr2pidVr*(np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:])))/(3*CONST.Q) # fixed matmul argument order - nh
        vth_local=0.1*np.sqrt(2*Tmaxwell[k]*CONST.Q/(mol*mu*CONST.H_MASS))
        Terror=abs(Tmaxwell[k]-T_out1)/Tmaxwell[k]
        Verror=abs(vx_out1-vx_shift[k])/vth_local

      # Compute desired moments

      WxD = copy.copy(vx_shift[k])
      ED = (WxD**2)+3*CONST.Q*Tmaxwell[k]/(mol*mu*CONST.H_MASS)

      # Compute present moments of Maxwell, WxMax, and EMax 

      max_2d_xd = np.matmul(maxwell[:, :, k],(vx * dVx))            
      WxMax       = vth  * (np.nansum(Vr2pidVr * (max_2d_xd) ))
      EMax        = vth2 * (np.nansum(Vr2pidVr*(np.matmul((vr2vx2_2D*maxwell[:,:,k]),dVx))))
      # print("Include Test", (vx*dVx))
      # print("Include Test", maxwell[k,:,:])
      # print("Vx", vx)
      # print("Include Test", np.matmul((vx*dVx),maxwell[k,:,:]))
      # print(vth)
      # print("Include Test", WxMax)
      # print("Include Test", EMax)
      # input()

      # Compute Nij from Maxwell, padded with zeros
      Nij = np.zeros((nvr+2,nvx+2), dtype=np.float64)  # The order of the error is maintained with float64, but not for float 32
      Nij[1:nvr+1, 1:nvx+1]   = maxwell[:,:,k]*Vol
      auxNij                  = Nij*Vx_DVx
      Nijp1_vx_Dvx            = np.roll(auxNij,  shift=-1, axis=1)
      Nij_vx_Dvx              = Nij*Vx_DVx
      Nijm1_vx_Dvx            = np.roll(auxNij,  shift= 1, axis=1)
      auxNij2                 = Nij*Vr_DVr
      Nip1j_vr_Dvr            = np.roll(auxNij2, shift=-1, axis=0)
      Nij_vr_Dvr              = Nij*Vr_DVr
      Nim1j_vr_Dvr            = np.roll(auxNij2, shift= 1, axis=0)
      # Compute Ap, Am, Bp, and Bm (0=p 1=m)
      aux_AN = np.zeros((nvx+2,nvr+2),dtype=np.float64)
      aux_AN = Nij*Vth_DVx

      _AN         = np.zeros((nvx+2,nvr+2), dtype=np.float64)
      _AN         = np.roll(aux_AN, shift=1, axis=1) - aux_AN
      AN[:,:,0]   = copy.copy(_AN[1:nvr+1,1:nvx+1])
      
      _AN         = np.zeros((nvx+2,nvr+2), dtype=np.float64)
      _AN         = -np.roll(aux_AN, shift=-1, axis=1) + aux_AN
      AN[:,:,1]   = copy.copy(_AN[1:nvr+1,1:nvx+1])

      BN[:,jpa+1:jpb+1,0] =  Nijm1_vx_Dvx[1:nvr+1,jpa+2:jpb+2] - Nij_vx_Dvx[1:nvr+1,jpa+2:jpb+2]
      BN[:,jpa,0]         = -Nij_vx_Dvx[1:nvr+1,jpa+1]
      BN[:,jnb,0]         =  Nij_vx_Dvx[1:nvr+1,jnb+1]
      BN[:,jna:jnb,0]     = -Nijp1_vx_Dvx[1:nvr+1,jna+1:jnb+1] + Nij_vx_Dvx[1:nvr+1,jna+1:jnb+1]
      BN[:,:,0]           =  BN[:,:,0] + Nim1j_vr_Dvr[1:nvr+1,1:nvx+1] - Nij_vr_Dvr[1:nvr+1,1:nvx+1]

      BN[:,jpa+1:jpb+1,1] = -Nijp1_vx_Dvx[1:nvr+1,jpa+2:jpb+2] + Nij_vx_Dvx[1:nvr+1,jpa+2:jpb+2]
      BN[:,jpa,1]         = -Nijp1_vx_Dvx[1:nvr+1,jpa+1]
      BN[:,jnb,1]         =  Nijm1_vx_Dvx[1:nvr+1,jnb+1]
      BN[:,jna:jnb,1]     =  Nijm1_vx_Dvx[1:nvr+1,jna+1:jnb+1] - Nij_vx_Dvx[1:nvr+1,jna+1:jnb+1]
      BN[1:nvr,:,1]       =  BN[1:nvr,:,1] - Nip1j_vr_Dvr[2:nvr+1,1:nvx+1] + Nij_vr_Dvr[2:nvr+1,1:nvx+1]
      BN[0,:,1]           =  BN[0,:,1] - Nip1j_vr_Dvr[1,1:nvx+1]

      # print("AN", AN.T)
      # print("BN", BN.T)
      # input()

      # Remove padded zeros in Nij

      Nij = Nij[1:nvr+1,1:nvx+1]

      # Cycle through 4 possibilies of sign(a_Max),sign(b_Max)

      TB1=np.zeros(2, float)
      TB2=np.zeros(2, float)
      ia=0
      while ia<2:

        # Compute TA1, TA2

        aux_TA1 = np.sum(np.matmul(AN[:, :, ia], vx))
        TA1 = vth*aux_TA1
        TA2 = vth2*np.sum(vr2vx2_2D*AN[:,:,ia])
        # print("maxwell", maxwell)
        # print("TA1", TA1)
        # print("TA2", TA2)
        ib=0
        while ib<2:

          # Compute TB1, TB2

          if TB1[ib]==0:
            aux_TB1 = np.sum(np.matmul(BN[:, :, ib], vx))
            TB1[ib] = vth*aux_TB1
          if TB2[ib]==0:
            TB2[ib] = vth2*np.sum(vr2vx2_2D*BN[:,:,ib])
          denom=TA2*TB1[ib]-TA1*TB2[ib]
          #print("TB1", TB1) NOTE This value is off

          b_max=0
          a_max=0
          if denom!=0 and TA1!=0:
            b_max=(TA2*(WxD-WxMax)-TA1*(ED-EMax))/denom
            a_max=(WxD-WxMax-TB1[ib]*b_max)/TA1
            # print("a_max", a_max)
            # print("WxD", WxD)
            # print("WxMax", WxMax)
            # print("TB1[ib]", TB1[ib])
            # print("b_max", b_max)
            # print("TA1", TA1)
            # input()
            # NOTE Some of these values are still off, but maxwell seems to be working for now
          if a_max*sgn[ia]>0 and b_max*sgn[ib]>0:
            maxwell[:,:,k] = (Nij + AN[:,:,ia]*a_max + BN[:,:,ib]*b_max)/Vol
            #print("maxwell", maxwell[k,:,:])
            ia=2
            ib=2
          ib+=1
        ia+=1
      aux_maxwell = np.sum(Vr2pidVr * (np.matmul(maxwell[:, :, k], dVx)))
      maxwell[:,:,k] = maxwell[:,:,k]/aux_maxwell
      
      # print("Maxwell iter", maxwell.T)

      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out2=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out2/vth)**2
        T_out2=(mol*mu*CONST.H_MASS)*vth2*np.sum(Vr2pidVr*np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:]))/(3*CONST.Q) # fixed matmul argument order - nh
        Terror2=abs(Tmaxwell[k]-T_out2)/Tmaxwell[k]
        Verror2=abs(vx_shift[k]-vx_out2)/vth_local
        print('CREATE_SHIFTED_MAXWELLIAN=> Terror:'+sval(Terror)+'->'+sval(Terror2)+'  Verror:'+sval(Verror)+'->'+sval(Verror2))

  # print("maxwell final", maxwell.T)
  # input()
  return maxwell.T
