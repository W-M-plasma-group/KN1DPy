;+
; Sigma_EL_H_HH.pro
;
;
; Returns momentum transfer cross section for elastic collisions of H onto H2 
; for specified energy of H. Data are taken from 
; 
; Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
; Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
; Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
; Passovets, page 305.
;
;________________________________________________________________________________
   Function Sigma_EL_H_HH,E
;________________________________________________________________________________
;  Input:
;	E	- fltarr(*) or float, energy of H atom (target H2 molecule is at rest)
;
;  Output:
;	Returns Sigma for 0.03 < E < 1e4. For E outside this range, 
;       the value of Sigma at the 0.03 or 1e4 eV boundary is returned.
;
;	units: m^-2
;________________________________________________________________________________
;-
   ty=type_of(E,nDim=nDim)
   _E=float([E])
   _E = _E > 0.03e0
   _E = _E < 1.01e4
   a=[-3.495671e1, -4.062257e-1, -3.820531e-2, -9.404486e-3, 3.963723e-4]
   result=EXP(poly(ALOG(_E),a))*1e-4
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
