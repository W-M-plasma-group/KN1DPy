;
; test_sigma_cx_hh.pro
;
      E=10.0^(-1+(5+alog10(2))*findgen(101)/100)
      sig=sigma_cx_hh(E)
         plot,/nodata,E,sig*1e4,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
              yrange=[1.0e-16,1.e-14],ystyle=1,$
              title='H!D2!N!U+!N + H!D2!N -> H!D2!N + H!D2!N!U+!N',xtitle='E (eV)',$
              ytitle='Sigma (cm!U2!N)',xticklen=1,yticklen=1
      oplot,E,sig*1e4,color=4,thick=3.
      xyouts,.15,.20,/normal,'output from idl function SIGMA_CX_HH.PRO which evaluates the CX Cross-Section using',charsize=0.8
      xyouts,.15,.17,/normal,'polynomial fit from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 253.',charsize=0.8,$
              color=4
end
