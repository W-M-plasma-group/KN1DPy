;+
; SigmaV_P_Hn2_HP.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen ions resulting in 
; one proton and one H atom in the n=2 state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.260.
;
; Also returns minimum, maximum, and average energy of the resultant proton and H(n=2) atom.
;
;________________________________________________________________________________
   Function SigmaV_P_Hn2_HP,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;  Output Keywords:
;	E0_ave	- float, average energy of P, H(n=2) atom (eV).
;	E0_max	- float, maximum energy of P, H(n=2) atom (eV).
;	E0_min	- float, minimum energy of P, H(n=2) atom (eV).
;________________________________________________________________________________
;-
   E0_ave=1.5
   E0_max=1.5
   E0_min=1.5
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-3.408905929046e+1, $
         1.573560727511e+1, $
        -6.992177456733e+0, $
         1.852216261706e+0, $
        -3.130312806531e-1, $
         3.383704123189e-2, $
        -2.265770525273e-3, $
         8.565603779673e-5, $
        -1.398131377085e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
