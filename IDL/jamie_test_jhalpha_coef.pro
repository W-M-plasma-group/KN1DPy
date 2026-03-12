;+
; jamie_test_jhalpha_coef.pro
;
; Test and plot Johnson–Hinnov recombination (α) coefficients by calling JHAlpha_Coef.
;+
pro jamie_test_jhalpha_coef
  common JH_Coef,DKnot,TKnot,order,                          $
         LogR_BSCoef,LogS_BSCoef,LogAlpha_BSCoef,            $
         A_Lyman,A_Balmer

  ;— 1) Define densities & temperature grid —
  nt = 1000
  TEMP = findgen(nt) / (nt-1) * (700.0 - 1.0) + 1.0              ; 1→700 eV
  dens_m3 = [1e18, 1e19, 1e20]                                  ; m^-3 = [1e12,1e13,1e14 cm^-3]

  ;— 2) Loop over densities: compute α, print & stash for plotting —
  for i=0,2 do begin
    d = dens_m3[i]
    density = replicate(d, nt)
    alpha = JHAlpha_Coef(density, TEMP, no_null=1)               ; m^3/s

    print, 'density:', d/1e6, 'cm^-3'
    print, 'Alpha_vals (m^3/s):', alpha

    case i of
      0: alpha1 = alpha
      1: alpha2 = alpha
      2: alpha3 = alpha
    endcase
  endfor

  ;— 3) Draw black axes & plot the three α‐curves (converted to cm^3/s) —
  !X.TICKLEN = 1
  !Y.TICKLEN = 1

  dbplot, /LOGX, /LOGY, TEMP, alpha3*1e6, /NODATA,    $
         XRANGE=[0.1,2000.0], YRANGE=[1e-16,1e-11],   $ ; adjust yrange as needed
         TITLE='e^- + H^+ → H + 2 e^- Recombination',  $
         XTITLE='Te (eV)',                            $
         YTITLE='<α>  cm^3 s^-1'

  dboplot, TEMP, alpha1*1e6, COLOR=2, THICK=2           ; red:   1e12 cm^-3
  dboplot, TEMP, alpha2*1e6, COLOR=3, THICK=2           ; blue:  1e13 cm^-3
  dboplot, TEMP, alpha3*1e6, COLOR=4, THICK=2           ; green: 1e14 cm^-3

  xyouts, /NORMAL, .15, .85, 'n!De!N = 10!U12!N cm!U-3!N', COLOR=2
  xyouts, /NORMAL, .15, .80, 'n!De!N = 10!U13!N cm!U-3!N', COLOR=3
  xyouts, /NORMAL, .15, .75, 'n!De!N = 10!U14!N cm!U-3!N', COLOR=4

  ;— 4) Citation in black —
  yleg = 0.25
  dy   = 0.03
  cs   = 0.8
  xyouts, .42, yleg,      /NORMAL, 'data from collisional-radiative model',                             CHARSIZE=cs
  xyouts, .42, yleg-dy,   /NORMAL, 'of L.C. Johnson and E. Hinnov, J. Quant. Spectrosc. Radiat. Transfer', CHARSIZE=cs
  xyouts, .42, yleg-2*dy, /NORMAL, 'vol.13 pp.333–358',                                                      CHARSIZE=cs
end
