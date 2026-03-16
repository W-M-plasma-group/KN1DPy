;
; jamie_test_locate.pro
;
; Test the LOCATE function for ascending and descending tables
;
pro jamie_test_locate

  ;— 1 Ascending table tests —
  table1 = [1, 2, 3, 4]
  values = [3.5,  3.0,  0.5,  4.5]
  idx1 = LOCATE(table1, values)
  print, 'Ascending table:', table1
  print, 'Values:           ', values
  print, 'LOCATE (IDL):     ', idx1
  ; expected: [2, 2, -1, 3]

  ;— 2 Descending table tests —
  table2 = [4, 3, 2, 1]
  idx2 = LOCATE(table2, values)
  print, 'Descending table:', table2
  print, 'Values:            ', values
  print, 'LOCATE (IDL):      ', idx2
  ; expected: [3, 3, 0, -1]

  ;— 3 Scalar lookup tests —
  scalar_val = 2.7
  idx_s1 = LOCATE(table1, scalar_val)
  idx_s2 = LOCATE(table2, scalar_val)
  print, 'Scalar ', scalar_val, ' in ascending =>', idx_s1  ; expect 1
  print, 'Scalar ', scalar_val, ' in descending=>', idx_s2  ; expect 2

end
