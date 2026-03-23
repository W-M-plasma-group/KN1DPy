;
; Jamie_Test_Balmer_Alpha.pro
;
tek
; 
; density scan test
;
density = 10^(16.0 + (21.0 - 16.0) * findgen(101) / 100.0)
Te      = replicate(20.0, 101)
n0      = replicate(3.3e16, 101)
bal     = Balmer_Alpha(density, Te, n0, /no_null)
dbplot, density, bal, /logx, /logy, $
       title = 'Balmer-alpha', $
       xtitle = 'Density (m^-3)', $
       ytitle = 'watts m^-3'
press_return

;
; Te scan test
;
Te      = 10^(alog10(0.35) + (alog10(700.0) - alog10(0.35)) * findgen(101) / 100.0)
density = replicate(1.0e19, 101)
n0      = replicate(3.3e16, 101)
bal     = Balmer_Alpha(density, Te, n0, /no_null)
dbplot, Te, bal, /logx, /logy, $
       title = 'Balmer-alpha', $
       xtitle = 'Te (eV)', $
       ytitle = 'watts m^-3'
press_return

end
