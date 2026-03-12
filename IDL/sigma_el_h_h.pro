;+
; Sigma_EL_H_H.pro
;
;
; Returns momentum transfer cross section for elastic collisions of H onto H 
; for specified energy of H. Data are taken from 
; 
; Janev, "Atomic and Molecular Processes in Fusion Edge Plasmas", Chapter 11 - 
; Elastic and Related Cross Sections for Low-Energy Collisions among Hydrogen and 
; Helium Ions, Neutrals, and Isotopes  by D.R. Sdhultz, S. Yu. Ovchinnikov, and S.V.
; Passovets, page 305.
;
;________________________________________________________________________________
   Function Sigma_EL_H_H,E,vis=vis
;________________________________________________________________________________
;  Input:
;	E	- fltarr(*) or float, energy of H atom (target H atom is at rest)
;
;  KEYWORD input:
;	VIS	- if set, then return viscosity cross section instead of momentum
;		  transfer cross section
;
;  Output:
;	Returns Sigma for 0.03 < E < 1e4. For E outside this range, 
;       the value of Sigma at the 0.03 or 1e4 eV boundary is returned.
;
;	units: m^-2
;________________________________________________________________________________
;-
   key_default,vis,0
   ty=type_of(E,nDim=nDim)
   _E=float([E])
   _E = _E > 0.03e0
   _E = _E < 1.01e4
   if vis then begin
      a=[-3.344860e1, -4.238982e-1, -7.477873e-2, -7.915053e-3, -2.686129e-4]
      result=EXP(poly(ALOG(_E),a))*1e-4
   endif else begin
      a=[-3.330843e1, -5.738374e-1, -1.028610e-1, -3.920980e-3, 5.964135e-4]
      result=EXP(poly(ALOG(_E),a))*1e-4
   endelse
   if nDim eq 0 then result=result(0)
   RETURN,result
END 
