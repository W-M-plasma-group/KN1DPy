;
; Create_Kinetic_H2_Mesh.pro
;
pro Create_Kinetic_H2_Mesh,nv,mu,x,Ti,Te,n,PipeDia,xH2,TiH2,TeH2,neH2,PipeDiaH2,vx,vr,Tnorm,E0=E0,ixE0=ixE0,irE0=irE0,fctr=fctr

   mH=1.6726231e-27		
   q=1.602177e-19				
   k_boltz=1.380658e-23				&;Bolzmann's constant, J K^-1
   Twall=293.0*k_boltz/q			&;room temperature (eV)
   v0_bar=sqrt(8.0*Twall*q/(!pi*2*mu*mH))	&;directed random velocity of diatomic molecule
   key_default,fctr,1.0
   nx=n_elements(x)
;
;  Estimate interaction rate with side walls
;
   gamma_wall=dblarr(nx)
;   for k=0,nx-1 do begin
;      if PipeDia(k) gt 0 then gamma_wall(k)=2*sqrt(2*Ti(k)*q/(2*mH))/PipeDia(k)
;   endfor
;
;  Estimate total reaction rate for destruction of molecules and for interation with side walls
;
;   RR=n*sigmav_ion_HH(Te)+n*sigmav_H1s_H1s_HH(Te)+n*sigmav_H1s_H2s_HH(Te)+gamma_wall

   RR=n*sigmav_ion_HH(Te)+n*sigmav_H1s_H1s_HH(Te)+n*sigmav_H1s_H2s_HH(Te)
;
; Determine x range for molecules by finding distance into plasma where density persists.
;
;  dGamma/dx=-nH2*RR = v0_bar dnH2/dx = -nH2*RR
;  d ln(nH2)/dx = -RR/v0_bar = dY/dx
;
;  Compute Y from RR and v0_bar
;
   Y=dblarr(nx)
   for k=1,nx-1 do Y(k)=Y(k-1)-(x(k)-x(k-1))*0.5*(RR(k)+RR(k-1))/v0_bar
;
; Find x location where Y = -10, i.e., where nH2 should be down by exp(-10)
;
   xmaxH2=interpol(x,Y,-10.0) < max(x)
   xminH2=x(0)
;
; Interpolate Ti and Te onto a fine mesh between xminH2 and xmaxH2
;
   xfine=xminH2+(xmaxH2-xminH2)*findgen(1001)/1000
   Tifine=interpol(Ti,x,xfine)
   Tefine=interpol(Te,x,xfine)
   nfine=interpol(n,x,xfine)
   PipeDiafine=interpol(PipeDia,x,xfine)
;
; Setup a vx,vr mesh based on raw data to get typical vx, vr values
;
   Create_VrVxMesh,nv,Tifine,vx,vr,Tnorm,E0=E0,ixE0=ixE0,irE0=irE0
   Vth=sqrt(2*q*Tnorm/(mu*mH))

;
;  Estimate interaction rate with side walls
;
   nxfine=n_elements(xfine)
   gamma_wall=dblarr(nxfine)
   for k=0,nxfine-1 do begin
      if PipeDiaFine(k) gt 0 then gamma_wall(k)=2*max(vr)*vth/PipeDiaFine(k)
   endfor
;
;
; Estimate total reaction rate, including charge exchange, elastic scattering, and interaction with side walls
;
   RR=nfine*sigmav_ion_HH(Tefine)+nfine*sigmav_H1s_H1s_HH(Tefine)+nfine*sigmav_H1s_H2s_HH(Tefine)+$
      0.1*nfine*sigmav_cx_HH(Tifine,Tifine) + gamma_wall

;
; Compute local maximum grid spacing from dx_max = 2 min(vr) / RR
;
   big_dx=0.02*fctr
   dx_max=fctr*0.8*(2*vth*min(vr)/RR) < big_dx
;
; Construct xH2 axis
;
   xpt=xmaxH2
   xH2=[xpt]
   while xpt gt xminH2 do begin
      xH2=[xpt,xH2]
      dxpt1=interpol(dx_max,xfine,xpt)
      dxpt2=dxpt1
      xpt_test=xpt-dxpt1
      if xpt_test gt xminH2 then dxpt2=interpol(dx_max,xfine,xpt_test)
      dxpt=min([dxpt1,dxpt2])
      xpt=xpt-dxpt
   endwhile
   xH2=[xminH2,xH2(0:n_elements(xH2)-2)]
;   if xH2(1)-xH2(0) > 0.5*big_dx then xH2=[xH2(0),xH2(2:*)]

   TiH2=interpol(Tifine,xfine,xH2)
   TeH2=interpol(Tefine,xfine,xH2)
   neH2=interpol(nfine,xfine,xH2)
   PipeDiaH2=interpol(PipeDiafine,xfine,xH2)
   Create_VrVxMesh,nv,TiH2,vx,vr,Tnorm,E0=E0,ixE0=ixE0,irE0=irE0
   return
   end
