;
; Test_KN1d_no_pipe.pro
;
   if !d.name ne 'PS' then tek
   newfile=1
;   plot=1
;   debug=1
   refine=1
;   pause=1
   compute_errors=1
   debrief=1
;   Hplot=1
   Hdebrief=1
;   H2plot=1
   H2debrief=1
   mu=2
   nx=8
   xw=0.0
   xlim=0.2
   xa=xlim-.02
   xb=0.265
   x1=xw+(xa-xw)*findgen(nx)/nx

   nx=100
   x2=xa+(xb-xa)*findgen(nx)/(nx-1)
   x=[x1,x2]
   nx=n_elements(x)

   Ti=10.0*exp((x-xlim)/.025)
   Ti=Ti > 1.5
   Te=Ti
   n=1.0e19*exp((x-xlim)/.03)
   n=n > 1.0e15
   n=n < 2.0e20

   nx=n_elements(x)
   LC=fltarr(nx)
   LC(where(x le xlim))=1.1
   vxi=fltarr(nx)
   xlimiter=0.2
   xsep=.25
;
   key_default,GaugeH2,1.0
   key_default,file,'Test_1mtorr_no_pipe'
   PipeDia=fltarr(nx)
;   PipeDia(where(x le xlimiter))=1.0
;   PipeDia(where(x le 0.05))=0.2

;   common KN1D_collisions,H2_H2_EL,H2_P_EL,H2_H_EL,H2_HP_CX,H_H_EL,H_P_EL,H_P_CX,Simple_CX
;   H2_H2_EL=0
;   H_H_EL=0

   KN1D,x,xlimiter,xsep,GaugeH2,mu,Ti,Te,n,vxi,LC,PipeDia,$
       xH2,nH2,GammaxH2,TH2,qxH2_total,nHP,THP,SH,SP,$
       xH,nH,GammaxH,TH,qxH_total,NetHSource,Sion,QH_total,SideWallH,Lyman,Balmer,$
       GammaHLim,$
       truncate=truncate,refine=refine,File=File,NewFile=NewFile,ReadInput=ReadInput,$
       error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause,$
       Hplot=Hplot,Hdebug=Hdebug,Hdebrief=Hdebrief,Hpause=Hpause,$
       H2plot=H2plot,H2debug=H2debug,H2debrief=H2debrief,H2pause=H2pause
end
