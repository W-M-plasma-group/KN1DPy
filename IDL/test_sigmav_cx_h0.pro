;
; test_sigmav_cx_h0.pro
;
Ti=10.0^(-1+(5+alog10(2))*findgen(101)/100)
for i=0,5 do begin
   e0=10.0^(i-1)
   E=replicate(e0,n_elements(Ti))
   sigv=sigmav_cx_h0(Ti,E)
   if i eq 0 then begin
      plot,/nodata,Ti,sigv*1e6,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
           yrange=[1.0e-9,2.e-6],ystyle=1,$
           title='p + H(1s) -> H(1s) + p',xtitle='Ti (eV)',$
           ytitle='<sigma v> (cm^3/s)',xticklen=1,yticklen=1
   endif
   oplot,Ti,sigv*1e6,color=i+2,thick=3.
   xyouts,.15,.7+i*.04,/normal,'E0 = '+sval(e0)+' eV',color=i+2
   xyouts,.15,.18,/normal,'output from idl function SIGMAV_CX_H0.PRO which evaluates the reaction rates using',charsize=0.8
   xyouts,.15,.15,/normal,'data from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 272.',charsize=0.8
endfor
end
