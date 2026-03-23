;
; Test_Kinetic_H2.pro
;
   tek
   key_default,plot,1
   key_default,debug,0
   key_default,debrief,1
   key_default,pause,0
   key_default,H2_H2_EL,1
   key_default,H2_P_EL,1
   key_default,H2_H_EL,1
   key_default,H2_HP_CX,1
   key_default,ni_correct,1
   key_default,Simple_CX,1
   key_default,Truncate,1.0e-4
   key_default,Compute_H_Source,1
   key_default,mu,2
   print,'H2_H2_EL:'+sval(H2_H2_EL)
   print,'H2_P_EL:'+sval(H2_P_EL)
   print,'H2_H_EL:'+sval(H2_H_EL)
   print,'H2_HP_CX:'+sval(H2_HP_CX)
   print,'ni_correct:'+sval(ni_correct)
   print,'Simple_CX:'+sval(Simple_CX)
   print,'mu:'+sval(mu)
   print,'Truncate:'+sval(truncate)

   compute_errors=1
   nx=8
   xw=0.0
   xlim=0.2
   xa=xlim-.02
   xb=0.25
   x1=xw+(xa-xw)*findgen(nx)/nx

   nx=100
   x2=xa+(xb-xa)*findgen(nx)/(nx-1)
   x=[x1,x2]
   nx=n_elements(x)

   Ti=10.0*exp((x-xlim)/.025)
   Ti=Ti > 1.0
   Te=Ti
   n=1.0e19*exp((x-xlim)/.030)
   n=n > 1.0e17
   PipeDia=fltarr(nx)
   PipeDia(*)=0.5

   nv=6
   Eneut=[0.003,0.01,0.03,0.1,0.3,1.0,3.0]
   fctr=0.8
   Create_Kinetic_H2_Mesh,nv,mu,x,Ti,Te,n,PipeDia,xH2,TiH2,TeH2,neH2,PipeDiaH2,vx,vr,Tnorm,E0=Eneut,ixE0=ixE0,irE0=irE0,fctr=fctr

   ip=where(vx gt 0)
   nvr=n_elements(vr)
   nvx=n_elements(vx)
   nx=n_elements(xH2)
   vxi=fltarr(nx)
   Nuloss=fltarr(nx)

   fH2BC=fltarr(nvr,nvx)
   Tneut=1.0/40
   for i=ip(0),nvx-1 do begin
      arg=-(vr(*)^2+vx(i)^2)/(Tneut/Tnorm/2)
      fH2BC(*,i)=exp(arg > (-80))
; divide Tneut by 2 to account for velocity of diatomic molecule = sqrt(2) of atom
   endfor
   GammaxH2BC=1.0e23

   multiplot,xh2,neh2,xH2,teH2,xH2,tiH2,title='Inputted profiles',xtitle='x',$
             ytitle=['n','Te','Ti'],color=[2,3,4],$
             varlabel=['Density (m^-3)','Te (ev)','Ti (eV)']
   press_return

   Kinetic_H2,vx,vr,xH2,Tnorm,mu,TiH2,TeH2,neH2,vxi,fH2BC,GammaxH2BC,NuLoss,PipeDiaH2,fH,SH2,$
       fH2,nH2,GammaxH2,VxH2,pH2,TH2,qxH2,qxH2_total,Sloss,QH2,RxH2,QH2_total,AlbedoH2,WallH2,$
       truncate=truncate,$
       nHP,THP,fSH,SH,SP,SHP,NuE,NuDis,$
       Simple_CX=Simple_CX,Max_Gen=Max_Gen,Compute_H_Source=Compute_H_Source,$
       No_Sawada=No_Sawada,H2_H2_EL=H2_H2_EL,H2_P_EL=H2_P_EL,H2_H_EL=_H2_H_EL,H2_HP_CX=H2_HP_CX,ni_correct=ni_correct,$
       ESH=ESH,Eaxis=Eaxis,error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause

   common Kinetic_H2_Output,piH2_xx,piH2_yy,piH2_zz,RxH2CX,RxH_H2,RxP_H2,RxW_H2,EH2CX,EH_H2,EP_H2,EW_H2,Epara_PerpH2_H2

   common Kinetic_H2_Errors,Max_dx,vbar_error,mesh_error,moment_error,C_Error,CX_Error,Swall_Error,H2_H2_error,Source_Error,$
                            qxH2_total_error,QH2_total_error

   end
