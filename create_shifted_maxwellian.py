import numpy as np
from sval import sval
from Make_dVr_dVx import Make_dVr_dVx

def create_shifted_maxwellain(vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm):
  nx=vx_shift.size
  nvr,nvx=vr.size,vx.size
  maxwell=np.zeros((nvr,nvx,nx)).T
  vr2vx2_ran2=np.zeros((nvr,nvx)).T

  mH=1.6726231e-27
  q=1.602177e-19

  vth=np.sqrt(2*q*Tnorm/(mu*mH))
  vth2=vth**2
  vth3=vth**3

  shifted_maxwellian_debug=0

  Vr2pidVr,VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb = Make_dVr_dVx(vr,vx)

  # create_shifted_maxwellian.include

  AN=BN=np.zeros((nvx,nvr,2)).T
  sgn=[-1,1]
  for k in range(nx):
    if Tmaxwell[k]>0:
      for i in range(nvr):
        arg=-(vr[i]**2+(vx-vx_shift[k]/vth)**2) * mol*Tnorm/Tmaxwell[k]
        arg=arg[arg<0]
        maxwell[k,:,i]=np.e**(arg[arg>-80])

      maxwell[k,:,:]/=np.sum(Vr2pidVr*(np.matmul(maxwell[k,:,:],dVx)))

      if shifted_maxwellian_debug:
        vx_out1=vth*np.sum(Vr2pidVr*np.matmul(maxwell[k,:,:],(vx*dVx)))
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out1/vth)**2
        T_out1=(mol*mu*mH)*vth2*np.sum(Vr2pidVr*(np.matmul(vr2vx2_ran2*maxwell[k,:,:],dVx)))/(3*q)
        vth_local=0.1*np.sqrt(2*Tmaxwell[k]*q/(mol*mu*mH))
        Terror=abs(Tmaxwell[k]-T_out1)/Tmaxwell[k]
        Verror=abs(vx_out1-vx_shift[k])/vth_local

      # Compute desired moments

      WxD=vx_shift[k]
      ED=WxD**2+3*q*Tmaxwell[k]/(mol*mu*mH)

      # Compute present moments of Maxwell, WxMax, and EMax 

      WxMax=vth*np.sum(Vr2pidVr*np.matmul(maxwell[k,:,:],(vx*dVx)))
      EMax=vth2*np.sum(Vr2pidVr*np.matmul((vr2vx2_2D*maxwell[k,:,:]),dVx))

      # Compute Nij from Maxwell, padded with zeros

      nij=np.zeros((nvr+2,nvx+2)).T
      nij[1:nvx+1,1:nvr+1]=maxwell[k,:,:]*Vol

      Nijp1_vx_Dvx=np.roll(nij*Vx_DVx,-1,0)
      Nij_vx_Dvx  =nij*Vx_DVx
      Nijm1_vx_Dvx=np.roll(nij*Vx_DVx,1,0)
      Nip1j_vr_Dvr=np.roll(nij*Vr_DVr,-1,1)
      Nij_vr_Dvr  =nij*Vr_DVr
      Nim1j_vr_Dvr=np.roll(nij*Vr_DVr,1,1)

      # Compute Ap, Am, Bp, and Bm (0=p 1=m)

      _AN=np.roll(nij*Vth_DVx,1,0)-nij*Vth_DVx
      AN[0,:,:]=_AN[1:nvx+1,1:nvr+1]
      _AN=-np.roll(nij*Vth_DVx,-1,0)+nij*Vth_DVx
      AN[0,:,:]=_AN[1:nvx+1,1:nvr+1]

      BN[0,jpa+1:jpb+1,:]=Nijm1_vx_Dvx[jpa+2:jpb+2,1:nvr+1]-Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[0,jpa,:]=-Nij_vx_Dvx[jpa+2,1:nvr+1]
      BN[0,jnb,:]=Nij_vx_Dvx[jnb+2,1:nvr+1]
      BN[0,jna:jnb,:]=-Nijp1_vx_Dvx[jna+1:jnb+1,1:nvr+1]+Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[0,:,:]=BN[0,:,:]+ Nim1j_vr_Dvr[1:nvx+1,1:nvr+1]-Nij_vr_Dvr[1:nvx+1,1:nvr+1]

      BN[1,jpa+1:jpb+1,:]=-Nijp1_vx_Dvx[jpa+2:jpb+2,1:nvr+1]+Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[1,jpa,:]=-Nijp1_vx_Dvx[jpa+1,1:nvr+1]
      BN[1,jnb,:]=Nijm1_vx_Dvx[jnb+1,1:nvr+1]
      BN[1,jna:jnb,:]=Nijm1_vx_Dvx[jna+1:jnb+1,1:nvr+1]-Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[1,:,1:nvr]=BN[1,:,1:nvr] - Nip1j_vr_Dvr[1:nvx+1,2:nvr+1]+Nij_vr_Dvr[1:nvx+1,2:nvr+1]
      BN[1,:,0]=BN[1,:,0] - Nip1j_vr_Dvr[1:nvx+1,1]

      # Remove padded zeros in Nij

      nij=nij[1:nvx+1,1:nvr+1]

      # Cycle through 4 possibilies of sign(a_Max),sign(b_Max)

      TB1=TB2=np.zeros(2)
      ia=0
      while ia<2:

        # Compute TA1, TA2

        TA1=vth*np.sum(np.matmul(AN[ia,:,:],vx))
        TA2=vth2*np.sum(vr2vx2_2D*AN[ia,:,:])
        ib=0
        while ib<2:

          # Compute TB1, TB2

          if TB1[ib]==0:
            TB1[ib]=vth*np.sum(np.matmul(BN[ib,:,:],vx))
          if TB2[ib]==0:
            TB2[ib]=vth2*np.sum(vr2vx2_2D*BN[ib,:,:])
          denom=TA2*TB1[ib]-TA1*TB2[ib]
          b_max=a_max=0
          if denom!=0 and TA1!=0:
            b_Max=(TA2*(WxD-WxMax)-TA1*(ED-EMax))/denom
            a_Max=(WxD-WxMax-TB1[ib]*b_Max)/TA1
          if a_max*sgn[ia]>0 and b_max*sgn[ib]>0:
            maxwell[k,:,:]=(nij+AN[ia,:,:]*a_Max+BN[ib,:,:]*b_Max)/Vol
            ia=ib=2
          ib+=1
        ia+=1
      maxwell[k,:,:]=maxwell[k,:,:]/np.sum(Vr2pidVr*np.matmul(maxwell[k,:,:],dVx))

      if shifted_maxwellian_debug:
        vx_out2=vth*np.sum(Vr2pidVr*np.matmul(maxwell[k,:,:],(vx*dVx)))
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out2/vth)**2
        T_out2=(mol*mu*mH)*vth2*np.sum(Vr2pidVr*np.matmul(vr2vx2_ran2*maxwell[k,:,:],dVx))/(3*q)
        Terror2=abs(Tmaxwell[k]-T_out2)/Tmaxwell[k]
        Verror2=abs(vx_shift[k]-vx_out2)/vth_local
        print('CREATE_SHIFTED_MAXWELLIAN=> Terror:'+sval(Terror)+'->'+sval(Terror2)+'  Verror:'+sval(Verror)+'->'+sval(Verror2))

  return maxwell