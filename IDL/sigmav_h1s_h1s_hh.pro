;+
; SigmaV_H1s_H1s_HH.pro
;
; Returns maxwellian averaged <sigma V> for electron impact
; dissociation of molecular hydrogen resulting in two H atoms in
; the 1s state. Coefficients are taken from Janev, 
; "Elementary Processes in Hydrogen-Helium Plasmas",
; Springer-Verlag, 1987, p.259.
;
;  NOTE: The coefficients in this edition of Janev's book for Reaction A2.5 are incorrectly listed on page 259.
;        Darren Stotler and Doug Post have since provided the correct coefficients (March 3, 1995).
;        The corrected peak reaction rate is now 9.81e-9 instead of the value 3.82e-09 cm^3/s.
;        The plot on page 45 also reflects the incorrect coefficients.
;
; This routine uses Darren Stotler and Doug Post's corrected coefficients.
;
; The routine also returns minimum, maximum, and average energy of the resultant H(1s) atoms.
;
;________________________________________________________________________________
   Function SigmaV_H1s_H1s_HH,Te,E0_ave=E0_ave,E0_min=E0_min,E0_max=E0_max
;________________________________________________________________________________
;  Input:
;	Te	- fltarr(*) or float, electron temperature (eV)
;
;  Output:
;	returns <sigma V> for 0.1 < Te < 2e4.
;	units: m^3/s
;
;  Output Keywords:
;	E0_ave	- float, average energy of H(1s) atoms (eV).
;	E0_max	- float, maximum energy of H(1s) atoms (eV).
;	E0_min	- float, minimum energy of H(1s) atoms (eV).
;________________________________________________________________________________
;-
   E0_ave=3.0
   E0_max=4.25
   E0_min=2.0
   t=type_of(Te,nDim=nDim)
   _Te=[Te]
   b = [-2.787217511174e+1, $
         1.052252660075e+1, $
        -4.973212347860e+0, $
         1.451198183114e+0, $
        -3.062790554644e-1, $
         4.433379509258e-2, $
        -4.096344172875e-3, $
         2.159670289222e-4, $
        -4.928545325189e-6]
   _Te =_Te > 0.1
   _Te =_Te < 2.01e4
   result=EXP(poly(ALOG(_Te),b))*1e-6
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
