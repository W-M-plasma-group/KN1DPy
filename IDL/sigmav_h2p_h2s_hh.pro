;+
; SigmaV_H2p_H2s_HH.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen resulting in one H atom in
; the 2p state and one H atom in the 2s state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.259.
;
; Also returns minimum, maximum, and average energy of the resultant H(2p), H(2s) atoms.
;
;________________________________________________________________________________
   Function SigmaV_H2p_H2s_HH,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
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
   E0_ave=4.85
   E0_max=5.85
   E0_min=2.85
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-4.794288960529e+1, $
         2.629649351119e+1, $
        -1.151117702256e+1, $
         2.991954880790e+0, $
        -4.949305181578e-1, $
         5.236320848415e-2, $
        -3.433774290547e-3, $
         1.272097387363e-4, $
        -2.036079507592e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
