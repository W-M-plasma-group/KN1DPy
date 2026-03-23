;+
; run_mastu_test_idl.pro
;
; MAST-U KN1D example - IDL
; ==========================
; Run from the KN1DPy root directory:
;
;   .run examples/MAST-U/run_mastu_test_idl.pro
;
; Inputs are loaded from the shared .sav file (same source as the Python script).
; To generate mastu_test_in.sav from the raw .dat files, run create_mastu_sav.pro
; in IDL once (see examples/MAST-U/create_mastu_sav.pro).
;-

pro run_mastu_test_idl

  ; ---------------------------------------------------------------- ;
  ;  Input file
  ; ---------------------------------------------------------------- ;
  sav_file = 'examples/MAST-U/input/mastu_test_in.sav'

  print, 'Loading inputs: ' + sav_file
  restore, sav_file
  ; Variables now in memory: x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia

  ; Unit conversions (sav file stores keV and 1e20 m^-3)
  Ti = Ti * 1e3    ; keV -> eV
  Te = Te * 1e3    ; keV -> eV
  n  = n  * 1e20   ; 1e20 m^-3 -> m^-3

  ; ---------------------------------------------------------------- ;
  ;  Run options
  ; ---------------------------------------------------------------- ;
  File           = 'examples/MAST-U/IDL_output/mastu_test'
  refine         = 0
  ReadInput      = 0
  NewFile        = 1
  compute_errors = 1
  Hdebrief       = 1
  H2debrief      = 1
  debrief        = 1

  ; ---------------------------------------------------------------- ;
  ;  Call KN1D
  ; ---------------------------------------------------------------- ;
  print, 'Starting KN1D ...'
  t0 = systime(1)

  KN1D, x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, $
        File=File, refine=refine, ReadInput=ReadInput, NewFile=NewFile, $
        compute_errors=compute_errors, $
        debrief=debrief, Hdebrief=Hdebrief, H2debrief=H2debrief

  print, string(systime(1)-t0, format='("Elapsed time: ", F6.1, " s")')
  print, 'Results written to ' + File + '.KN1D_H etc.'

end
