;+
; jamie_test_make_sigv.pro
;
; Test Make_SigmaV with Sigma_Function='SQRT' (σ(E)=√E)
;+
PRO jamie_test_make_sigv

  ;— 1 Define the same test grid —
  mE = 3
  Emin = 1.0  & Emax = 100.0
  E_Particle = 10.0^(alog10(Emin) + (alog10(Emax)-alog10(Emin)) * findgen(mE)/(mE-1))

  nT = 2
  Tmin = 1.0  & Tmax = 100.0
  T_Target = 10.0^(alog10(Tmin) + (alog10(Tmax)-alog10(Tmin)) * findgen(nT)/(nT-1))

  mu_particle = 1.0
  mu_target   = 1.0

  PRINT, 'E_Particle (eV):', E_Particle

  ;— 2 Loop over T_Target calling Make_SigmaV with SQRT —
  FOR i=0, nT-1 DO BEGIN
    PRINT, ' T_Target (eV):', T_Target[i]
    ; pass 'SQRT' so that σ(E_rel)=SQRT(E_rel)
    SigmaV = Make_SigmaV(E_Particle, mu_particle, T_Target[i], mu_target, 'SQRT')
    PRINT, ' SigmaV (IDL) [m^2/s]:', SigmaV
  ENDFOR

END
