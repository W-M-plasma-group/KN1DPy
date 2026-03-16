;+
; SigmaV_H1s_Hn_HP.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; dissociative recombination of molecular hydrogen ions resulting in 
; one H atom in the 1s state and one H atom in state n > or = 2. Coefficients 
; are taken from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.260.
;
;________________________________________________________________________________
   Function SigmaV_H1s_Hn_HP,Te
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;________________________________________________________________________________
;-
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-1.670435653561e+1, $
        -6.035644995682e-1, $
        -1.942745783445e-8, $
        -2.005952284492e-7, $
         2.962996104431e-8, $
         2.134293274971e-8, $
        -6.353973401838e-9, $
         6.152557460831e-10, $
        -2.025361858319e-11]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
