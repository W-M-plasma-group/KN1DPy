pro jamie_test_jhs_coef
  common JH_Coef,DKnot,TKnot,order,LogR_BSCoef,LogS_BSCoef,LogAlpha_BSCoef,A_Lyman,A_Balmer

  ;— 1 Define densities & temperatures —
  dens_m3 = [1e18, 1e19, 1e20]              ; [m^-3] = [10^12,10^13,10^14 cm^-3]
  nt      = 1000
  TEMP    = findgen(nt) / (nt-1) * (700.0-1.0) + 1.0  ; 1→700 eV

  ;— 2 Loop: compute S, print density, TEMP, S_vals —
  for i=0,2 do begin
    d = dens_m3[i]
    density = replicate(d, nt)
    S       = JHS_Coef(density, TEMP, no_null=1) * 1e6  ; cm^3/s

    ; Print exactly as in Python:
    print, 'density:', d/1e6, ' cm^-3'      ; convert back to cm^-3 for readability
    print, 'TEMP array (eV):', TEMP
    print, 'S_vals (cm^3/s):', S

    ; Save into separate vectors for plotting
    case i of
      0: sv1 = S
      1: sv2 = S
      2: sv3 = S
    endcase
  endfor

  ;— 3 Now do your black‐axes + coloured curves plot —
  !X.TICKLEN = 1
  !Y.TICKLEN = 1

  dbplot, /LOGX, /LOGY, TEMP, sv3, /NODATA, $
         XRANGE=[0.1,2000.0], YRANGE=[1e-11,1e-7], XSTYLE=1, $
         TITLE='e^- + H → H^+ + 2 e^-', $
         XTITLE='Te (eV)', $
         YTITLE='<σv>  cm^3 s^-1'

  dboplot, TEMP, sv1, COLOR=2, THICK=2    ; red   (10^12)
  dboplot, TEMP, sv2, COLOR=3, THICK=2    ; blue  (10^13)
  dboplot, TEMP, sv3, COLOR=4, THICK=2    ; green (10^14)

  xyouts, /NORMAL, .15, .85, 'n!De!N = 10!U12!N cm!U-3!N', COLOR=2
  xyouts, /NORMAL, .15, .80, 'n!De!N = 10!U13!N cm!U-3!N', COLOR=3
  xyouts, /NORMAL, .15, .75, 'n!De!N = 10!U14!N cm!U-3!N', COLOR=4

  yleg = 0.25
  dy   = 0.03
  cs   = 0.8
  xyouts, .42, yleg,      /NORMAL, 'data from collisional-radiative model',                             CHARSIZE=cs
  xyouts, .42, yleg-dy,   /NORMAL, 'of L.C. Johnson and E. Hinnov, J. Quant. Spectrosc. Radiat. Transfer', CHARSIZE=cs
  xyouts, .42, yleg-2*dy, /NORMAL, 'vol.13 pp.333–358',                                                      CHARSIZE=cs
end




