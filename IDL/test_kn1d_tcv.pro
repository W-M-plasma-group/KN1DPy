;
; Test_KN1d_tcv.pro
;
   if !d.name ne 'PS' then tek
   newfile=0
   plot=1
   pause=1
   debug=1
   refine=0
   key_default,pause,0
   compute_errors=1
   debrief=1
   Hdebrief=1
   H2debrief=1
   mu=2
   nx=8
   xw=0.0
   xlim=0.01
   xa=xlim-.002
   xb=0.4
   x1=xw+(xa-xw)*findgen(nx)/nx

;   nx=100
   nx=500
   x2=xa+(xb-xa)*findgen(nx)/(nx-1)
   x=[x1,x2]
   nx=n_elements(x)

   Ti=400*(1-((max(x)-x)/max(x))^2)^1.5+1.
; this is OK !!! (AK)
;   Te=800*(1-((max(x)-x)/max(x))^2)^1.5+1.
; !!!!! NOK !!!!! (AK)
   Te=900*(1-((max(x)-x)/max(x))^2)^1.5+1.
;   Te=810*(1-((max(x)-x)/max(x))^2)^1.5+1.

   Ti=Ti > 5.
   Te=Te > 5.
;   Te=Te < 800.
;   Te=Te < 900.
;   Te=Ti
   n=2e19*(1-((max(x)-x)/max(x))^2)^1.+1e17
   n=n > 1.0e15
   n=n < 2.0e20


   nx=n_elements(x)
;   GaugeH2=.01
   GaugeH2=.1
   LC=fltarr(nx)
   LC(where(x le xlim))=1.1
   vxi=fltarr(nx)
   xlimiter=0.02
   xsep=.02
;
   file=''
   PipeDia=fltarr(nx)
;
;   PipeDia(where(x le xlimiter))=1.0
;   PipeDia(where(x le 0.05))=0.1
;
; Turn off self collisions
;
    common KN1D_collisions,H2_H2_EL,H2_P_EL,H2_H_EL,H2_HP_CX,H_H_EL,H_P_EL,H_P_CX,Simple_CX
   H_H_EL=0
   H2_H2_EL=0


;       H2_H2_EL- if set, then include H2 -> H2 elastic self collisions
;        H2_P_EL- if set, then include H2 -> H(+) elastic collisions 
;        H2_H_EL- if set, then include H2 <-> H elastic collisions 
;       H2_HP_CX- if set, then include H2 -> H2(+) charge exchange collisions
;  H_H_EL- if set, then include H -> H elastic self collisions
;  H_P_CX- if set, then include H -> H(+) charge exchange collisions 
;         H_P_EL- if set, then include H -> H(+) elastic collisions 
;             Simple_CX- if set, then use CX source option (B): Neutrals are born
;                         in velocity with a distribution proportional
;                         to the local
;                         ion distribution function. Simple_CX=1 is default.

   KN1D,x,xlimiter,xsep,GaugeH2,mu,Ti,Te,n,vxi,LC,PipeDia,$
       xH2,nH2,GammaxH2,TH2,qxH2_total,nHP,THP,SH,SP,$

xH,nH,GammaxH,TH,qxH_total,NetHSource,Sion,QH_total,SideWallH,Lyman,Balmer,$
       GammaHLim,$

truncate=truncate,refine=refine,File=File,NewFile=NewFile,ReadInput=ReadInput,$
       error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause,$
       Hplot=Hplot,Hdebug=Hdebug,Hdebrief=Hdebrief,Hpause=Hpause,$
       H2plot=H2plot,H2debug=H2debug,H2debrief=H2debrief,H2pause=H2pause

;
saveoutput,x,xlimiter,xsep,gaugeh2,mu,lc,n,te,ti,xh2,nh2,th2,GammaxH2,xh,nh,th,gammaxh,Lyman,Balmer,qxH2_total,SP,SH,THP,nHP,qxH_total,Sion

end
