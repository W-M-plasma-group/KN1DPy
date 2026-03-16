;
; plot_sigmav_EL_H_P.pro
;
Emin=0.1
Emax=2e4
Ti=10.0^(alog10(Emin)+(alog10(Emax)-alog10(Emin))*findgen(101)/100)
for i=0,5 do begin
   e0=10.0^(i-1)
   E=replicate(e0,n_elements(Ti))
   sigv=sigmav_EL_H_P(Ti,E)
   if i eq 0 then begin
      plot,/nodata,Ti,sigv*1e6,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
           yrange=[1.0e-9,2.e-6],ystyle=1,$
           title='p + H -> p + H elastic',xtitle='Ti (eV)',$
           ytitle='<sigma v> (cm!U3!N s!U-1!N)',xticklen=1,yticklen=1
   endif
   oplot,Ti,sigv*1e6,color=i+2,thick=3.
   xyouts,.15,.7+i*.04,/normal,'E0 = '+sval(e0)+' eV',color=i+2
endfor
   yleg=.25
   dy=.03
   cs=0.8
   xyouts,.15,yleg,/normal,'output from idl function SIGMAV_EL_P_H.PRO which evaluates the reaction rates using',charsize=0.8
   xyouts,.15,yleg-dy,/normal,'data from "Atomic and Molecular Processes in Fusion Edge Plasmas", Edited by R.K. Janev,',charsize=cs
   xyouts,.15,yleg-2*dy,/normal,'Chapter 1: Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen',charsize=cs
   xyouts,.15,yleg-3*dy,/normal,'and Helium Ions, Neutrals, and Isotopes  by D.R. Schultz, S.Yu. Ovchinnikov, and ',charsize=cs
   xyouts,.15,yleg-4*dy,/normal,'S.V. Passovets, pages 298. (Plenum Press, New York 1995)',charsize=cs
end
