;+
; SigmaV_P_P_HP.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociation of molecular hydrogen ions resulting in 
; two protons. Coefficients are taken from Janev, 
; "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.260.
;
;________________________________________________________________________________
   Function SigmaV_P_P_HP,Te
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
   b = [-3.746192301092e+1, $
         1.559355031108e+1, $
        -6.693238367093e+0, $
         1.981700292134e+0, $
        -4.044820889297e-1, $
         5.352391623039e-2, $
        -4.317451841436e-3, $
         1.918499873454e-4, $
        -3.591779705419e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
