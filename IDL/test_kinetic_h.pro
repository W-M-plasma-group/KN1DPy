;
; Test_kinetic_H.pro
;
   tek
   key_default,name,'test'
   file='Kinetic_H_'+name+'.dat'
   key_default,plot,1
   key_default,debug,0
   key_default,debrief,1
   key_default,pause,1
   key_default,truncate,1.e-4
   key_default,test,1
   key_default,simple_cx,1 ;1
   key_default,H_H_EL, 1 ;1
   key_default,H_P_EL, 1 ;1
   key_default,H_H2_EL, 1 ;1
   key_default,H_P_CX, 1 ;1
   key_default,ni_correct,1 ;1
   key_default, No_Recomb, 0 ;0
   key_default, No_Johnson_Hinnov, 0 ;0
   key_default,mu,2
   key_default,Truncate,1.0e-4
   print,'H_H_EL:'+sval(H_H_EL)
   print,'H_P_EL:'+sval(H_P_EL)
   print,'H_H2_EL:'+sval(H_H2_EL)
   print,'H_P_CX:'+sval(H_P_CX)
   print,'ni_correct:'+sval(ni_correct)
   print,'Simple_CX:'+sval(Simple_CX)
   print,'mu:'+sval(mu)
   print,'Truncate:'+sval(truncate)
   compute_errors=1

   test = 1

;
; Test #1 - 3 eV neutrals
;
   if test eq 1 then begin
      nx=8
      xw=0.0
      xlim=0.2
      xa=xlim-.02
      xb=0.28
      x1=xw+(xa-xw)*findgen(nx)/nx

      nx=100
      x2=xa+(xb-xa)*findgen(nx)/(nx-1)
      x=[x1,x2]
      nx=n_elements(x)

      Ti=10.0*exp((x-xlim)/.025)
      Ti=Ti > 1.0
      Te=Ti
      n=1.0e19*exp((x-xlim)/.03)
      n=n > 1.0e15
      n=n < 2.0e20
 
      nv=10
      fctr=1.0
      PipeDia=dblarr(nx)
      PipeDia(*)=0.5
      Create_Kinetic_H_Mesh,nv,mu,x,Ti,Te,n,PipeDia,xH,TiH,TeH,neH,PipeDiaH,vx,vr,Tnorm,fctr=fctr


      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(xH)
      vxi=fltarr(nx)

      GammaxHBC=1.0e23
      fHBC=fltarr(nvr,nvx)
      Tneut=3.0
      for i=ip(0),nvx-1 do begin
         arg=-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm)
         fHBC(*,i)=exp(arg > (-80))
      endfor
  endif
  if test eq 2 then begin
;
; Test #2 no ionization, Ti=T0
;
      max_gen=100
      nx=70
      xa=0.0
      xb=0.05
      x=xa+(xb-xa)*findgen(nx)/(nx-1)

      Ti=replicate(10.0,nx)

      nv=20
      Create_VrVxMesh,nv,Ti,vx,vr,Tnorm

      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(x)
      vxi = fltarr(nx)

      GammaxHBC=1.0e23
;      n=1.0e19*exp(x/.030)
      n=replicate(5.0e19,nx)
      Te=replicate(0.1,nx)
      PipeDiaH=dblarr(nx)
      fHBC=fltarr(nvr,nvx)
      Tneut=Ti(0)
      for i=ip(0),nvx-1 do begin
         arg=-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm)
         fHBC(*,i)=exp(arg > (-80))
      endfor

      xH = x
      TiH = Ti
      TeH = Te
      neH = n
   endif
  if test eq 3 then begin
;
; Test #3 with large ionization/charge exchange fraction, Ti=T0 
;
      max_gen=100
      nx=70
      xa=0.0
      xb=0.05
      x=xa+(xb-xa)*findgen(nx)/(nx-1)

      Ti=replicate(10.0,nx)

      nv=40
      Create_VrVxMesh,nv,Ti,nv,vr,Tnorm


      ip=where(vx gt 0)
      nvr=n_elements(vr)
      nvx=n_elements(vx)
      nx=n_elements(x)

      GammaxHBC=1.0e22
      n=replicate(5.0e19,nx)
      Te=replicate(30.0,nx)
      fHBC=fltarr(nvr,nvx)
      Tneut=Ti(0)
      for i=ip(0),nvx-1 do begin
         arg=-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm)
         fHBC(*,i)=exp(arg > (-80))
      endfor

      xH = x
      TiH = Ti
      TeH = Te
      neH = n
      STOP
   endif
;
   print,'Test ='+sval(test)

   multiplot,x,n,x,te,x,ti,title='Inputted profiles',xtitle='x',$
             ytitle=['n','Te','Ti'],color=[2,3,4],$
             varlabel=['Density (m^-3)','Te (ev)','Ti (eV)']
   press_return
   plot,vx,fHBC(0,*),/nodata,yrange=[0,max(fHBC)],title='Inputted fHBC'
   for i=0,nvr-1 do oplot,vx,fHBC(i,*),color=(i mod 8)+2
   press_return

   ; print, fH[*,0,0]
   ; press_return
   

   Kinetic_H,vx,vr,xH,Tnorm,mu,TiH,TeH,neH,vxi,fHBC,GammaxHBC,PipeDiaH,fH2,fSH,nHP,THP,$
       fH, nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,WallH,$
       truncate=truncate,Simple_CX=Simple_CX,Max_Gen=Max_Gen,$
       No_Johnson_Hinnov=No_Johnson_Hinnov,No_Recomb=No_Recomb,$
       H_H_EL=H_H_EL,H_P_EL=H_P_EL,H_H2_EL=_H_H2_EL,H_P_CX=H_P_CX,ni_correct=ni_correct,$
       error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause

   common Kinetic_H_Output,piH_xx,piH_yy,piH_zz,RxHCX,RxH2_H,RxP_H,RxW_H,EHCX,EH2_H,EP_H,EW_H,Epara_PerpH_H,SourceH,SRecomb

   common Kinetic_H_Errors,Max_dx,vbar_error,mesh_error,moment_error,C_Error,CX_Error,H_H_error,$
                           qxH_total_error,QH_total_error

;   save,file=file,vx,vr,x,Tnorm,mu,Ti,Te,n,vxi,fHBC,GammaxHBC,fH2,fSH,nHP,THP,$
;       fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,$
;       truncate,Simple_CX,Max_Gen,$
;       No_Johnson_Hinnov,No_Recomb,$
;       H_H_EL,H_P_EL,H_H2_EL,H_P_CX,$
;       error,compute_errors,$
;       plot,debug,debrief,pause,$
;       piH_xx,piH_yy,piH_zz,RxH2_H,RxP_H,EH2_H,EP_H,$
;       Max_dx,vbar_error,mesh_error,moment_error,C_Error,H_H_error,$
;       qxH_total_error,QH_total_error

   end
