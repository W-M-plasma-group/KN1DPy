;+
; SigmaV_H1s_Hn3_HH.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen resulting in one H atom in
; the 1s state and one H atom in the n=3 state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.259.
;
; Also returns minimum, maximum, and average energy of the resultant H(1s), H(n=3) atoms.
;
;________________________________________________________________________________
   Function SigmaV_H1s_Hn3_HH,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;  Output Keywords:
;	E0_ave	- float, average energy of H(2p), H(2s) atoms (eV).
;	E0_max	- float, maximum energy of H(2p), H(2s) atoms (eV).
;	E0_min	- float, minimum energy of H(2p), H(2s) atoms (eV).
;________________________________________________________________________________
;-
   E0_ave=2.5
   E0_max=1.25
   E0_min=3.75
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-3.884976142596e+1, $
         1.520368281111e+1, $
        -6.078494762845e+0, $
         1.535455119900e+0, $
        -2.628667482712e-1, $
         2.994456451213e-2, $
        -2.156175515382e-3, $
         8.826547202670e-5, $
        -1.558890013181e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
