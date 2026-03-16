;+
; Sigma_EL_P_H.pro
;
;
; Returns momentum transfer cross section for elastic collisions of H+ onto H 
; for specified energy of H+. Data are taken from 
; 
; Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
; Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
; Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
; Passovets, page 298.
;
;________________________________________________________________________________
   Function Sigma_EL_P_H,E
;________________________________________________________________________________
;  Input:
;	E	- fltarr(*) or float, energy of H+ ion (target H atom is at rest)
;
;  Output:
;	Returns Sigma for 0.001 < E < 1e5. For E outside this range, 
;       the value of Sigma at the 0.001 or 1e5 eV boundary is returned.
;
;	units: m^-2
;________________________________________________________________________________
;-
   ty=type_of(E,nDim=nDim)
   _E=float([E])
   _E = _E > 0.001e0
   _E = _E < 1.01e5
   result=_E & result(*)=0.0
   ilow=where(_E le 10.0,count)
   if count gt 0 then begin
      a=[-3.233966e1, -1.126918e-1, 5.287706e-3, -2.445017e-3, -1.044156e-3, 8.419691e-5, 3.824773e-5]
      result(ilow)=EXP(poly(ALOG(_E(ilow)),a))*1e-4
   endif
   ihigh=where(_E gt 10.0,count)
   if count gt 0 then begin
      a=[-3.231141e1, -1.386002e-1]
      result(ihigh)=EXP(poly(ALOG(_E(ihigh)),a))*1e-4
   endif
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
