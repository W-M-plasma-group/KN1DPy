;
; test_sigmav_rec_H1s.pro
;
   Te=10.0^(-1+(5+alog10(2))*findgen(101)/100)
   sigv=sigmav_rec_H1s(Te)
      plot,/nodata,Te,sigv*1e6,/xlog,/ylog,xrange=[0.1,2e4],xstyle=1,$
           yrange=[1.0e-11,1.e-7],ystyle=1,$
           title='e + H!U+!N -> H(1s) + hv',xtitle='Te (eV)',$
           ytitle='<sigma v> (cm^3/s)',xticklen=1,yticklen=1
   oplot,Te,sigv*1e6,color=2,thick=3.0
   xyouts,.15,.90,/normal,'output from idl function SIGMAV_REC_H1s.PRO which evaluates the reaction rates using',charsize=0.8
   xyouts,.15,.87,/normal,'data from JANEV et al.,"Elementary Processes in Hydrogen-Helium Plasmas", p 32.',charsize=0.8
   end
