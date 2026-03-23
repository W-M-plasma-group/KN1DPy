;+
; Sigma_cx_HH.pro
;
; Returns charge exchange cross section for molecular hydrogen. Data are taken
; the polynomial fit in
;
;     Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.253.
;
;________________________________________________________________________________
   Function Sigma_CX_HH,E
;________________________________________________________________________________
;  Input:
;	E	- fltarr(*) or float, energy of molecule corresponding to the
;                 relative velocity between molecule and molecular ion. (eV)
;  Output:
;	returns sigma_CX for 0.1 < E < 2e4
;	units: m^-2
;________________________________________________________________________________
;-
   ty=type_of(E,nDim=nDim)
   _E=double([E])
   _E = _E > 0.1d0
   _E = _E < 2.01D4
   alpha=[-3.427958758517D+01, -7.121484125189D-02, 4.690466187943D-02,$
          -8.033946660540D-03, -2.265090924593D-03,-2.102414848737D-04,$
           1.948869487515D-04, -2.208124950005D-05, 7.262446915488D-07]
   result=EXP(poly(ALOG(_E),alpha))*1D-4
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
