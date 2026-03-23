;+
; jamie_test_create_shifted_maxwellian.pro
;   Standalone test of CREATE_SHIFTED_MAXWELLIAN + NuLoss + extra moments
;________________________________________________________________________________
pro jamie_test_create_shifted_maxwellian

  ;— physical constants —
  mH = 1.6726231D-27
  q  = 1.602177D-19

  ;— 1) Define some simple test meshes —
  nvrM = 20
  nvxM = 31
  vrM   = findgen(nvrM) * 5.0d0 / (nvrM - 1)         ; 0 → 5 m/s
  vxM   = findgen(nvxM) * 10.0d0 / (nvxM - 1) - 5.0d0 ; −5 → +5 m/s

  ;— 2) Positive‐vx indices —
  ipM = where(vxM gt 0.0d0, n_ip)

  ;— 3) Gas parameters for H₂ beam —
  GaugeH2     = 0.1d0
  v0_bar      = 1.0d5    ; m/s
  DensM       = 3.537D19 * GaugeH2
  GammaxH2BC  = 0.25d0 * DensM * v0_bar

  ;— 4) Shifted Maxwellian inputs —
  mu        = 2.0d0     ; reduced‐mass
  mol       = 2         ; diatomic
  Twall     = 2.0d3     ; eV
  Tmaxwell  = [Twall]
  vx_shift  = [0.0d0]
  TnormM    = total(Tmaxwell) / n_elements(Tmaxwell)

  ;— 5) Build and call —
  create_shifted_maxwellian, vrM, vxM, Tmaxwell, vx_shift, mu, mol, TnormM, Maxwell

  ;— 6) Beam slice —
  fH2BC = fltarr(nvrM, nvxM)
  fH2BC[*, ipM] = Maxwell[*, ipM, 0]

  ;— 7) Toy plasma→NuLoss (unchanged) —
  nX   = 50
  x    = findgen(nX) * 0.1d0
  Ti   = 1.0d0 + 0.1d0 * findgen(nX)
  Te   = 2.0d0 + 0.1d0 * findgen(nX)
  LC   = 0.2d0 + 0.01d0 * findgen(nX)
  xH2  = x[where(x le 2.0d0, cnt)]
  Cs_LC = fltarr(nX)
  idx = where(LC gt 0.0d0, n_gt)
  if n_gt gt 0 then $
    Cs_LC[idx] = sqrt(q*(Ti[idx]+Te[idx])/(mu*mH))/LC[idx]
  NuLoss = interpol(Cs_LC, x, xH2)

  ;— 8) Now compute extra moments of the shifted‐Maxwellian —
  ;   get the velocity‐space weights
  Make_dVr_dVx, vrM, vxM, VR2pidVr, VrVr4pidVr, dVx, $
              vrL=vrL, vrR=vrR, vxL=vxL, vxR=vxR, $
              Vol=Vol, Vth_DeltaVx=vth_DVx, Vx_DeltaVx=vx_DVx, $
              Vr_DeltaVr=Vr_DVr, Vr2Vx2=vr2vx2_2d


  Vth = sqrt(2*q*TnormM/(mu*mH))

  ; density
  dens  = total( VR2pidVr * ( Maxwell(*,*,0) # dVx ) )
  ; mean vx
  ux_out = Vth * total( VR2pidVr * ( Maxwell(*,*,0) # (vxM * dVx) ) )
  ; temperature (include molecular weight)
  T_out  = (mol*mu*mH)*Vth^2 * total( VR2pidVr * ((vr2vx2_2d * Maxwell(*,*,0)) # dVx) ) / (3*q)

  ;— 9) Print everything —
  print, '--- Test Results ---'
  print, 'DensM            =', DensM
  print, 'GammaxH2BC       =', GammaxH2BC
  print, 'NuLoss(', n_elements(xH2), 'points)=', NuLoss
  print, '--- shifted‐Maxwellian moments:'
  print, ' target  vx shift =', vx_shift(0), '   actual ux_out =', ux_out
  print, ' target Tmaxwell   =', Tmaxwell(0), '   actual T_out  =', T_out
  print, ' density norm      =', dens

end

