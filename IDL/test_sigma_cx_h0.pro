;
; test_sigma_cx_h0.pro
;
   key_default,both,0
   if both then begin
      freeman=1
      E=10.0^(-1+(5+alog10(2))*findgen(101)/100)
      sig=sigma_cx_h0(E,freeman=freeman)
         plot,/nodata,E,sig*1e4,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
              yrange=[1.0e-16,1.e-14],ystyle=1,$
              title='p + H(1s) -> H(1s) + p',xtitle='E (eV)',$
              ytitle='sigma (cm^2)',xticklen=1,yticklen=1
      oplot,E,sig*1e4,color=2,thick=3.
      freeman=0
      sig1=sigma_cx_h0(E,freeman=freeman)
      oplot,E,sig1*1e4,color=4,thick=3.
      xyouts,.15,.20,/normal,'output from idl function SIGMA_CX_H0.PRO which evaluates the CX Cross-Section using',charsize=0.8
      xyouts,.15,.17,/normal,'Freeman and Jones,"Atomic Collision Processes in Plasma Physics Experiments", table 2.',charsize=0.8,$
              color=2
      xyouts,.15,.14,/normal,'polynomial fit from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 250.',charsize=0.8,$
              color=4
   endif else begin
      freeman=0
      E=10.0^(-1+(5+alog10(2))*findgen(101)/100)
      sig=sigma_cx_h0(E,freeman=freeman)
         plot,/nodata,E,sig*1e4,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
              yrange=[1.0e-16,1.e-14],ystyle=1,$
              title='p + H(1s) -> H(1s) + p',xtitle='E (eV)',$
              ytitle='sigma (cm^2)',xticklen=1,yticklen=1
      oplot,E,sig*1e4,color=4,thick=3.
      xyouts,.15,.20,/normal,'output from idl function SIGMA_CX_H0.PRO which evaluates the CX Cross-Section using',charsize=0.8
      xyouts,.15,.17,/normal,'polynomial fit from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 250.',charsize=0.8,$
              color=4
   endelse
end
