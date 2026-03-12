;
; CREATE_SHIFTED_MAXWELLIAN.PRO
;
pro create_shifted_Maxwellian,vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm,Maxwell
   nx=n_elements(vx_shift)
   nvr=n_elements(vr)
   nvx=n_elements(vx)
   Maxwell=dblarr(nvr,nvx,nx)
   vr2vx2_ran2=dblarr(nvr,nvx)
   mH=1.6726231D-27
   q=1.602177D-19				
   Vth=sqrt(2*q*Tnorm/(mu*mH))
   vth2=vth*vth
   vth3=vth2*vth
   shifted_maxwellian_debug=0

   Make_dVr_dVx,vr,vx,Vr2pidVr,VrVr4pidVr,dVx,vrL=vrL,vrR=vrR,vxL=vxL,vxR=vxR,$
                Vol=Vol,Vth_DeltaVx=Vth_DVx,Vx_DeltaVx=Vx_DVx,Vr_DeltaVr=Vr_DVr,Vr2Vx2=Vr2Vx2_2D,$
                jpa=jpa,jpb=jpb,jna=jna,jnb=jnb

@create_shifted_maxwellian.include
   return
   end
