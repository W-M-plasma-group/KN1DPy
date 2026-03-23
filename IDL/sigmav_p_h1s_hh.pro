;+
; SigmaV_P_H1s_HH.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; ionization and disociation of Molecular hydrogen resulting in one 
; proton and one H atom in the 1s state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.260.
;
;________________________________________________________________________________
   Function SigmaV_P_H1s_HH,Te
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
   b = [-3.834597006782e+1, $
         1.426322356722e+1, $
        -5.826468569506e+0, $
         1.727940947913e+0, $
        -3.598120866343e-1, $
         4.822199350494e-2, $
        -3.909402993006e-3, $
         1.738776657690e-4, $
        -3.252844486351e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
