;+
; jamie_test_create_kinetic_h2_mesh.pro
;   Stand‐alone test of Create_Kinetic_H2_Mesh.pro
;_______________________________________________________________________________
pro jamie_test_create_kinetic_h2_mesh

  ;— 0 Physical constants (just for reference) —
  mH = 1.6726231D-27
  q  = 1.602177D-19
  kB = 1.380658D-23

  ;— 1 Fake plasma profile inputs —
  nv       = 12
  mu       = 2.0d0
  nX       = 8
  x        = findgen(nX)*0.15d0                    ; 0..1.05
  Ti       = 5.0d0 + 0.2d0*findgen(nX)
  Te       = 6.0d0 + 0.3d0*findgen(nX)
  n        = 1.0e19*(1 + 0.05d0*findgen(nX))
  PipeDia  = 0.2d0 + 0.01d0*findgen(nX)

  ;— 2 Call the IDL routine under test —
  Create_Kinetic_H2_Mesh, nv, mu, x, Ti, Te, n, PipeDia, fctr=1.0, $
        xH2, TiH2, TeH2, neH2, PipeDiaH2, vx, vr, Tnorm, $
        E0=E0, ixE0=ixE0, irE0=irE0

  ;— 3 Print out sizes and first/last entries —
  print, 'IDL → Create_Kinetic_H2_Mesh results:'
  print, ' xH2  (npts)=', n_elements(xH2), xH2(0), xH2(n_elements(xH2)-1)
  print, ' TiH2 (npts)=', n_elements(TiH2), TiH2(0), TiH2(n_elements(TiH2)-1)
  print, ' TeH2        =', TeH2(0), TeH2(n_elements(TeH2)-1)
  print, ' neH2        =', neH2(0), neH2(n_elements(neH2)-1)
  print, ' PipeDiaH2   =', PipeDiaH2(0), PipeDiaH2(n_elements(PipeDiaH2)-1)
  print, ' vx (len)=', n_elements(vx), '  vr (len)=', n_elements(vr)
  print, ' Tnorm      =', Tnorm
  print, ' E0         =', E0, ' ixE0=', ixE0, ' irE0=', irE0

end
