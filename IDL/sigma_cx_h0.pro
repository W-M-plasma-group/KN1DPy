;+
; Sigma_cx_H0.pro
;
; Returns charge exchange cross section for atomic hydrogen. Data are taken
; either from the polynomial fit in
;
;     Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987, p.250.
;
; from Freeman and Jone's analytic fit tabulated in 
;
;     Freeman, E.L., Jones, E.M., "Atomic Collision Processes in Plasma Physics
;     Experiments", UKAEA Report No. CLM-R137 (Culham Laboratory, Abington, England 1974)
;
;________________________________________________________________________________
   Function Sigma_CX_H0,E,freeman=freeman
;________________________________________________________________________________
;  Input:
;	E	- fltarr(*) or float, energy of proton corresponding to the
;                 relative velocity between proton and hydrogen atom. (eV)
;
;    Keywords:
;
;	Freeman - if set, then return CX based on Freeman and Jones' analytic fit in
;		  Freeman, E.L., Jones, E.M., "Atomic Collision Processes in Plasma Physics
;                 Experiments", UKAEA Report No. CLM-R137 (Culham Laboratory, Abington, England 1974)
;
;                 otherwise return CX based on polynomial fit in
;		  Janev, "Elementary Processes in Hydrogen-Helium Plasmas", 
;	          Springer-Verlag, 1987, p.250, other
;                 
;                 Default is Freeman=0
;
;  Output:
;	returns sigma_CX for 0.1 < E < 2e4
;	units: m^-2
;________________________________________________________________________________
;-
   key_default,freeman,0
   ty=type_of(E,nDim=nDim)
   if freeman then begin
      _E=double([E])
      _E = _E > 0.1D0
      _E = _E < 1.0D5
      result=1.0D-4 * 0.6937D-14*(1.0D0-0.155D0*ALOG10(_E))^2/(1.0D0+0.1112D-14*_E^3.3)
   endif else begin
      _E=double([E])
      _E = _E > 0.1d0
      _E = _E < 2.01D4
      alpha=[-3.274123792568D+01, -8.916456579806D-02, -3.016990732025D-02,$
              9.205482406462D-03,  2.400266568315D-03, -1.927122311323D-03,$
              3.654750340106D-04, -2.788866460622D-05,  7.422296363524D-07]
      result=EXP(poly(ALOG(_E),alpha))*1D-4
   endelse
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
