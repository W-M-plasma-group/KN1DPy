;
; jamie_test_path_interp_2d.pro
;
; IDL example: use a handful of distinct (n,Te) pairs
;
pro jamie_test_path_interp_2d

  ;— 1 Lookup tables (log‐space) —
  Te_table = alog([5, 20, 100])
  Ne_table = alog([1e14,1e17,1e18,1e19,1e20,1e21,1e22])

  fctr_Table = fltarr(7,3)
  fctr_Table(*,0) = [2.2,2.2,2.1,1.9,1.2,1.1,1.05]/5.3
  fctr_Table(*,1) = [5.1,5.1,4.3,3.1,1.5,1.25,1.25]/10.05
  fctr_Table(*,2) = [1.3,1.3,1.1,0.8,0.38,0.24,0.22]/2.1

  ;— 2 Pick a few distinct test points —
  Te = [10.0, 20.0, 50.0, 100.0]      ; eV
  n  = [1e15, 1e17, 1e19, 1e21]       ; m^-3

  ;— 4 Interpolate along the path —
  fctr = Path_Interp_2D(fctr_Table, Ne_table, Te_table, alog(n), alog(Te))

  ;— 5 Print with labels so you can compare directly —
  print, '   n(m^-3) :', n
  print, '   Te (eV) :', Te
  print, 'fctr (IDL):', fctr
end


