;
; test_sigmav_cx_hh.pro
;
Ti=10.0^(-1+(5+alog10(2))*findgen(101)/100)
for i=0,5 do begin
   e0=10.0^(i-1)
   E=replicate(e0,n_elements(Ti))
   sigv=sigmav_cx_hh(Ti,E)
   if i eq 0 then begin
      plot,/nodata,Ti,sigv*1e6,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
           yrange=[1.0e-10,2.e-7],ystyle=1,$
           title='H!D2!N!U+!N + H!D2!N -> H!D2!N + H!D2!N!U+!N',xtitle='T (eV)',$
           ytitle='<sigma v> (cm!U3!N s!U-1!N)',xticklen=1,yticklen=1
   endif
   oplot,Ti,sigv*1e6,color=i+2,thick=3.
   xyouts,.15,.7+i*.04,/normal,'E = '+sval(e0)+' eV',color=i+2
   xyouts,.15,.18,/normal,'output from idl function SIGMAV_CX_HH.PRO which evaluates the reaction rates using',charsize=0.8
   xyouts,.15,.15,/normal,'data from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 292.',charsize=0.8
endfor
end
