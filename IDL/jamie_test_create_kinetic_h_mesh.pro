;+
; jamie_test_create_kinetic_h_mesh.pro
;    Stand‐alone test of Create_Kinetic_H_Mesh.pro
;_______________________________________________________________________________
pro jamie_test_create_kinetic_h_mesh

  ;— 0) Choose fake inputs —
  nv      = 12
  mu      = 2.0d0
  nX      = 8
  x       = findgen(nX) * 0.15d0            ; 0 .. 1.05
  Ti      = 5.0d0 + 0.2d0 * findgen(nX)
  Te      = 6.0d0 + 0.3d0 * findgen(nX)
  n       = 1.0e19 * (1 + 0.05d0 * findgen(nX))
  PipeDia = 0.2d0 + 0.01d0 * findgen(nX)

  ;— 1) Call the routine under test —
  Create_Kinetic_H_Mesh, nv, mu, x, Ti, Te, n, PipeDia, $
        xH, TiH, TeH, neH, PipeDiaH, vx, vr, Tnorm, $
        E0=E0, ixE0=ixE0, irE0=irE0, fctr=1.0

  ;— 2) Print out summary —
  print, 'IDL → Create_Kinetic_H_Mesh results:'
  print, ' xH  (npts)=', n_elements(xH), xH(0), xH(n_elements(xH)-1)
  print, ' TiH (npts)=', n_elements(TiH), TiH(0), TiH(n_elements(TiH)-1)
  print, ' TeH        =', TeH(0), TeH(n_elements(TeH)-1)
  print, ' neH        =', neH(0), neH(n_elements(neH)-1)
  print, ' PipeDiaH   =', PipeDiaH(0), PipeDiaH(n_elements(PipeDiaH)-1)
  print, ' vx (len)=', n_elements(vx), '  vr (len)=', n_elements(vr)
  print, ' Tnorm      =', Tnorm
  print, ' E0         =', E0, ' ixE0=', ixE0, ' irE0=', irE0

end
