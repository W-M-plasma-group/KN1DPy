;+
; jamie_test_create_vrvxmesh.pro
;
; Test Create_VrVxMesh in three scenarios
;+
pro jamie_test_create_vrvxmesh

  ;-- common parameters --
  nv = 5
  Ti1 = [10.0, 20.0, 30.0]    ; a simple Ti profile

  ;-- 1 baseline: no E0, no Tmax --
  Create_VrVxMesh, nv, Ti1, vx1, vr1, Tnorm1
  print, '--- Test 1: no E0, no Tmax ---'
  print, 'nv=', nv, 'Ti=', Ti1
  print, 'vx1 =', vx1
  print, 'vr1 =', vr1
  print, 'Tnorm1 =', Tnorm1
  print

  ;-- 2 with E0 = 5 eV --
  E0 = [5.0]
  Create_VrVxMesh, nv, Ti1, vx2, vr2, Tnorm2, E0=E0, ixE0=ix2, irE0=ir2
  print, '--- Test 2: E0 =', E0, '---'
  print, 'vx2 =', vx2
  print, 'vr2 =', vr2
  print, 'ixE0 =', ix2, ' irE0 =', ir2
  print

  ;-- 3 with Tmax = 25 eV (drop Ti above 25) --
  Tmax = 25.0
  Create_VrVxMesh, nv, Ti1, vx3, vr3, Tnorm3, Tmax=Tmax
  print, '--- Test 3: Tmax =', Tmax, '---'
  print, 'Ti used =', Ti1[where(Ti1 lt Tmax)]
  print, 'vx3 =', vx3
  print, 'vr3 =', vr3
  print, 'Tnorm3 =', Tnorm3

end
