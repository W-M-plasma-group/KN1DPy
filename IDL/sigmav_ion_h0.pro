;+
; SigmaV_ion_H0.pro
;
; Returns maxwellian averaged <sigma V) for electron impact
; ionization of atomic hydrogen. Coefficients are taken
; from Janev, "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.258.
;
;________________________________________________________________________________
   Function SigmaV_Ion_H0,Te
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
   b = [-3.271396786375e+1, $
         1.353655609057e+1, $
        -5.739328757388e+0, $
         1.563154982022e+0, $
        -2.877056004391e-1, $
         3.482559773737e-2, $
        -2.631976175590e-3, $
         1.119543953861e-4, $
        -2.039149852002e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
