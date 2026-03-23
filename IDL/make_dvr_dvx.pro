;
; Make_dVr_dVx.pro
;
;   Constructs velocity space differentials for distribution functions
; used by Kinetic_Neutrals.pro, Kinetic_H2.pro, Kinetic_H2.pro, and other 
; related procedures.
;
;   03/25/2004 Bug in computation of Vr2pidVr and VrVr4pidVr found by Jerry Hughes and corrected
;          
;
pro Make_dVr_dVx,vr,vx,Vr2pidVr,VrVr4pidVr,dVx,vrL=vrL,vrR=vrR,vxL=vxL,vxR=vxR,$
                 Vol=Vol,Vth_DeltaVx=Vth_DeltaVx,Vx_DeltaVx=Vx_DeltaVx,Vr_DeltaVr=Vr_DeltaVr,Vr2Vx2=Vr2Vx2,$
                 jpa=jpa,jpb=jpb,jna=jna,jnb=jnb
;
; Determine velocity space differentials
;
   nvr=n_elements(vr)
   nvx=n_elements(vx)
   _vr=[vr,2*vr(nvr-1)-vr(nvr-2)]
   vr_mid=[0.0,0.5*(_vr+shift(_vr,-1))]
   vrR=shift(vr_mid,-1)
   vrL=vr_mid
;   Vr2pidVr=2*!pi*(vrR^2-vrL^2) - original line 3/25/2004 B. LaBombard
   Vr2pidVr=!pi*(vrR^2-vrL^2)
   Vr2pidVr=Vr2pidVr(0:nvr-1)
;   VrVr4pidVr=4*!pi*(vrR^3-vrL^3) - original line 3/25/2004 B. LaBombard
   VrVr4pidVr=(4/3.)*!pi*(vrR^3-vrL^3)
   VrVr4pidVr=VrVr4pidVr(0:nvr-1)
   vrR=vrR(0:nvr-1)
   vrL=vrL(0:nvr-1)

   _vx=[2*vx(0)-vx(1),vx,2*vx(nvx-1)-vx(nvx-2)]
   vxR=0.5*(shift(_vx,-1)+_vx)
   vxL=0.5*(shift(_vx,1)+_vx)
   dVx=vxR(1:nvx)-vxL(1:nvx)
   vxR=vxR(1:nvx)
   vxL=vxL(1:nvx)
;
; Compute volume elements
;
   vol=dblarr(nvr,nvx)
   for i=0,nvr-1 do vol(i,*)=Vr2pidVr(i)*dVx
;
; Compute Deltavx, Deltavr
;
       Deltavx=vxR-vxL
       Deltavr=vrR-vrL
;
; Compute vth_Deltavx, vx_Deltavx, vr_Deltavr, padded with zeros
;
      vth_Deltavx=fltarr(nvr+2,nvx+2)
      vx_Deltavx=fltarr(nvr+2,nvx+2)
      vr_Deltavr=fltarr(nvr+2,nvx+2)
      for i=1,nvr do begin
         vth_Deltavx(i,1:nvx)=1.0/Deltavx
         vx_Deltavx(i,1:nvx)=vx/Deltavx
      endfor
      for j=1,nvx do begin
         vr_Deltavr(1:nvr,j)=vr/Deltavr
      endfor
;
; Compute v^2
;
      vr2vx2=dblarr(nvr,nvx)
      for i=0,nvr-1 do vr2vx2(i,*)=vr(i)^2+vx^2
;
; Determine indice range of positive and negative vx
;
      jp=where(vx gt 0)
      jpa=jp(0) & jpb=jp(n_elements(jp)-1)
      jn=where(vx lt 0)
      jna=jn(0) & jnb=jn(n_elements(jn)-1)

   return
end
