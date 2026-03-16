;+
; SigmaV_P_H1s_HP.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen ions resulting in 
; one proton and one H atom in the 1s state. Coefficients are taken 
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.260.
;
; Also returns minimum, maximum, and average energy of the resultant proton and H(1s) atom.
;
;________________________________________________________________________________
   Function SigmaV_P_H1s_HP,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;  Output Keywords:
;	E0_ave	- float, average energy of P, H(1s) atom (eV).
;	E0_max	- float, maximum energy of P, H(1s) atom (eV).
;	E0_min	- float, minimum energy of P, H(1s) atom (eV).
;________________________________________________________________________________
;-
   E0_ave=4.3
   E0_max=4.3
   E0_min=4.3
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-1.781416067709e+1, $
         2.277799785711e+0, $
        -1.266868411626e+0, $
         4.296170447419e-1, $
        -9.609908013189e-2, $
         1.387958040699e-2, $
        -1.231349039470e-3, $
         6.042383126281e-5, $
        -1.247521040900e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
