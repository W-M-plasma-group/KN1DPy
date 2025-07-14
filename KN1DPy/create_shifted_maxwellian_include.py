import numpy as np

from .sval import sval
from .common import constants as CONST
import copy

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
  

  AN=np.zeros((2,nvx,nvr), float)
  BN=np.zeros((2,nvx,nvr), float) # fixed array creation - nh
  sgn=[1,-1]
  print(nx)
  for k in range(nx):
    if Tmaxwell[k]>0:
      for i in range(nvr):
        arg = -(vr[i]**2 + (vx - vx_shift[k]/vth)**2)*mol*Tnorm/Tmaxwell[k]
        arg = np.minimum(arg,0)
        maxwell[k,:,i] = np.e**(np.maximum(arg,-80)) # replaced < and > with np.minimum and np.maximum - nh

      maxwell[k,:,:] /= np.sum(Vr2pidVr*(np.matmul(dVx,maxwell[k,:,:]))) # fixed matmul argument order - nh

      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out1=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out1/vth)**2
        T_out1=(mol*mu*CONST.H_MASS)*vth2*np.sum(Vr2pidVr*(np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:])))/(3*CONST.Q) # fixed matmul argument order - nh
        vth_local=0.1*np.sqrt(2*Tmaxwell[k]*CONST.Q/(mol*mu*CONST.H_MASS))
        Terror=abs(Tmaxwell[k]-T_out1)/Tmaxwell[k]
        Verror=abs(vx_out1-vx_shift[k])/vth_local

      # Compute desired moments

      WxD=vx_shift[k]
      ED=WxD**2+3*CONST.Q*Tmaxwell[k]/(mol*mu*CONST.H_MASS)

      # Compute present moments of Maxwell, WxMax, and EMax 

      WxMax=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
      EMax=vth2*np.sum(Vr2pidVr*np.matmul(dVx,(vr2vx2_2D*maxwell[k,:,:]))) # fixed matmul argument order - nh
      # print("Include Test", (vx*dVx))
      # print("Include Test", maxwell[k,:,:])
      # print("Vx", vx)
      # print("Include Test", np.matmul((vx*dVx),maxwell[k,:,:]))
      # print(vth)
      # print("Include Test", WxMax)
      # print("Include Test", EMax)

      # Compute Nij from Maxwell, padded with zeros

      nij = np.zeros((nvr+2,nvx+2), float).T
      nij[1:nvx+1,1:nvr+1] = maxwell[k,:,:]*Vol

      #NOTE I Think a lot of these can be removed and/or combined
      #There are a lot of padded zeros, surely they are unnecessary?
      Nijp1_vx_Dvx = np.roll(nij*Vx_DVx, -1, axis=0)
      Nij_vx_Dvx   = nij*Vx_DVx
      Nijm1_vx_Dvx = np.roll(nij*Vx_DVx, 1, axis=0)
      Nip1j_vr_Dvr = np.roll(nij*Vr_DVr, -1, axis=1)
      Nij_vr_Dvr   = nij*Vr_DVr
      Nim1j_vr_Dvr = np.roll(nij*Vr_DVr, 1, axis=1)

      # Compute Ap, Am, Bp, and Bm (0=p 1=m)

      _AN = np.roll(nij*Vth_DVx,1,0)-nij*Vth_DVx
      AN[0,:,:] = _AN[1:nvx+1,1:nvr+1]
      _AN = -np.roll(nij*Vth_DVx,-1,0)+nij*Vth_DVx
      AN[1,:,:] = _AN[1:nvx+1,1:nvr+1] #NOTE Check these values later, once values hit ~-e40, round to 0, otherwise seems about right

      #NOTE This works, but check the positionings of the various jpa, jpb, nvr, etc for array indexing
      BN[0,jpa+1:jpb+1,:] = Nijm1_vx_Dvx[jpa+2:jpb+2,1:nvr+1] - Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[0,jpa,:] = -Nij_vx_Dvx[jpa+1,1:nvr+1]
      BN[0,jnb,:] = Nij_vx_Dvx[jnb+1,1:nvr+1]
      BN[0,jna:jnb,:] = -Nijp1_vx_Dvx[jna+1:jnb+1,1:nvr+1] + Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[0,:,:] = BN[0,:,:] + Nim1j_vr_Dvr[1:nvx+1,1:nvr+1] - Nij_vr_Dvr[1:nvx+1,1:nvr+1]

      BN[1,jpa+1:jpb+1,:] = -Nijp1_vx_Dvx[jpa+2:jpb+2,1:nvr+1] + Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[1,jpa,:] = -Nijp1_vx_Dvx[jpa+1,1:nvr+1]
      BN[1,jnb,:] = Nijm1_vx_Dvx[jnb+1,1:nvr+1]
      BN[1,jna:jnb,:] = Nijm1_vx_Dvx[jna+1:jnb+1,1:nvr+1] - Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[1,:,1:nvr] = BN[1,:,1:nvr] - Nip1j_vr_Dvr[1:nvx+1,2:nvr+1] + Nij_vr_Dvr[1:nvx+1,2:nvr+1]
      BN[1,:,0] = BN[1,:,0] - Nip1j_vr_Dvr[1:nvx+1,1]

      # Remove padded zeros in Nij

      nij=nij[1:nvx+1,1:nvr+1]

      # Cycle through 4 possibilies of sign(a_Max),sign(b_Max)

      TB1=np.zeros(2, float)
      TB2=np.zeros(2, float)
      ia=0
      while ia<2:

        # Compute TA1, TA2

        TA1=vth*np.sum(np.matmul(vx,AN[ia,:,:])) # fixed matmul argument order - nh
        TA2=vth2*np.sum(vr2vx2_2D*AN[ia,:,:])
        # print("maxwell", maxwell)
        # print("TA1", TA1)
        # print("TA2", TA2)
        ib=0
        while ib<2:

          # Compute TB1, TB2

          if TB1[ib]==0:
            TB1[ib]=vth*np.sum(np.matmul(vx,BN[ib,:,:])) # fixed matmul argument order - nh
          if TB2[ib]==0:
            TB2[ib]=vth2*np.sum(vr2vx2_2D*BN[ib,:,:])
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
            maxwell[k,:,:]=(nij+AN[ia,:,:]*a_max+BN[ib,:,:]*b_max)/Vol
            #print("maxwell", maxwell[k,:,:])
            ia=ib=2
          ib+=1
        ia+=1
      maxwell[k,:,:]=maxwell[k,:,:]/np.sum(Vr2pidVr*np.matmul(dVx,maxwell[k,:,:])) # fixed matmul argument order - nh

      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out2=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out2/vth)**2
        T_out2=(mol*mu*CONST.H_MASS)*vth2*np.sum(Vr2pidVr*np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:]))/(3*CONST.Q) # fixed matmul argument order - nh
        Terror2=abs(Tmaxwell[k]-T_out2)/Tmaxwell[k]
        Verror2=abs(vx_shift[k]-vx_out2)/vth_local
        print('CREATE_SHIFTED_MAXWELLIAN=> Terror:'+sval(Terror)+'->'+sval(Terror2)+'  Verror:'+sval(Verror)+'->'+sval(Verror2))

  # print("maxwell", maxwell)
  # input()
  return maxwell
