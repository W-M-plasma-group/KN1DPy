;
; jamie_test_normal2data.pro
;
; Test the NORMAL2DATA function for both linear and log axes
;
pro jamie_test_normal2data

  ;— 1 Linear‐axis test —
  norm_lin = [2.0, 4.0, 6.0]
  ; set up X‐axis metadata: offset=2.0, scale=2.0, linear type
  !X.S     = [2.0, 2.0]
  !X.type  = 0
  data_lin = normal2data(norm_lin, x=0)  
  print, 'Linear axis: norm=', norm_lin, '=> data=', data_lin
  ; expected: [(2–2)/2, (4–2)/2, (6–2)/2] = [0,1,2]

  ;— 2 Log‐axis test —
  norm_log = [2.0, 3.0, 4.0]
  ; set up Y‐axis metadata: offset=1.0, scale=1.0, log type
  !Y.S     = [1.0, 1.0]
  !Y.type  = 1
  data_log = normal2data(norm_log, y=1)
  print, 'Log axis:    norm=', norm_log, '=> data=', data_log
  ; expected: 10^((norm−1)/1) = [10, 100, 1000]

end
