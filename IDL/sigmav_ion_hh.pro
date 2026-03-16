;+
; SigmaV_Ion_HH.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; ionization of molecular hydrogen. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.259.
;
;________________________________________________________________________________
   Function SigmaV_Ion_HH,Te
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;________________________________________________________________________________
;-
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-3.568640293666e+1, $
         1.733468989961e+1, $
        -7.767469363538e+0, $
         2.211579405415e+0, $
        -4.169840174384e-1, $
         5.088289820867e-2, $
        -3.832737518325e-3, $
         1.612863120371e-4, $
        -2.893391904431e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
