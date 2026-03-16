;+
; SigmaV_H1s_H2s_HH.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen resulting in one H atom in
; the 1s state and one H atom in the 2s state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.259.
;
; Also returns minimum, maximum, and average energy of the resultant H(1s), H(2s) atoms.
;
;________________________________________________________________________________
   Function SigmaV_H1s_H2s_HH,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;  Output Keywords:
;	E0_ave	- float, average energy of H(1s), H(2s) atoms (eV).
;	E0_max	- float, maximum energy of H(1s), H(2s) atoms (eV).
;	E0_min	- float, minimum energy of H(1s), H(2s) atoms (eV).
;________________________________________________________________________________
;-
   E0_ave=0.3
   E0_max=0.55
   E0_min=0.0
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-3.454175591367e+1, $
         1.412655911280e+1, $
        -6.004466156761e+0, $
         1.589476697488e+0, $
        -2.775796909649e-1, $
         3.152736888124e-2, $
        -2.229578042005e-3, $
         8.890114963166e-5, $
        -1.523912962346e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
