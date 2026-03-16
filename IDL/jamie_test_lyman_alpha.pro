;
; Jamie_Test_Lyman_Alpha.pro
;
tek
; 
; density scan test
;
density=10^(16.0+(21.-16)*findgen(101)/100)
Te=replicate(20,101)
n0=replicate(3.3e16,101)
ly=lyman_alpha(density,te,n0)
dbplot,density,ly,/logx,/logy,title='Lyman-alpha',xtitle='Density (m^-3)',ytitle='watts m^-3'
press_return

;
; Te scan test
;
Te=10^(alog10(0.35)+(alog10(700)-alog10(0.35))*findgen(101)/100)
density=replicate(1.0e19,101)
n0=replicate(3.3e16,101)
ly=lyman_alpha(density,te,n0)
dbplot,Te,ly,/logx,/logy,title='Lyman-alpha',xtitle='Te (eV)',ytitle='watts m^-3'
press_return


end