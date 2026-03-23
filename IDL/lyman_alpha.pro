;+
; Lyman_Alpha.pro
;
;     Computes Lyman-alpha emissivity (watts m^-3) given the local
;  electron density, electron temperature, and ground-state
;  neutral density.
;
;  Method:
;     (1) Compute the local n=2 population density using the Johnson-Hinnov
; rate equations and coefficients [L.C.Johnson and E. Hinnov, J. Quant. 
; Spectrosc. Radiat. Transfer. vol. 13 pp.333-358]
;     (2) Multiply by the n=2->1 spontaneous emission coefficient
;     (3) Convert to watts/m^3
;
Function Lyman_Alpha,Density,Te,N0,photons=photons,create=create,no_null=no_null
;________________________________________________________________________________
; Input:
;  	Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
;  	Te	- fltarr, electron temperature (eV
;  	N0	- fltarr, ground state neutral density (m^-3)
;
; Keywords:
;	photons - returns emissivity in number of photons m^-3 s^-1
;	create	- if set, then create bi-cubic spline coefficients for
;		  interpolation of r0(p) r1(p) and save them in the
;		  default save set. 
;	No_Null	- if set, then rather than generate a NULL value when Density and Te
;                 are not null but still outside the data range, compute the rate based on the min or max
;		  data range values.
;________________________________________________________________________________
; History:
;    Coding by B. LaBombard  6/29/99
;    Coefficients from J. Terry's idl code JH_RATES.PRO
;   
;-
   common JH_Coef,DKnot,TKnot,order,LogR_BSCoef,LogS_BSCoef,LogAlpha_BSCoef,A_Lyman,A_Balmer
   if keyword_set(create) then create_JH_BSCoef
   if type_of(LogR_BSCoef) eq 0 then begin
      restore,'jh_bscoef.dat'
   endif
;
;  From Johnson and Hinnov, eq (11):
;
;   n(2) =  ( r0(2) + r1(2)*n(1)/NHsaha(1) )*NHsaha(2)
;
   if n_elements(Density) ne n_elements(Te) then message,'Number of elements of Density and Te are different!'
   if n_elements(Density) ne n_elements(N0) then message,'Number of elements of Density and N0 are different!'
   result=Density & result(*)=1.0e32
   photons=result
   r02=JHR_Coef(Density,Te,0,2,no_null=no_null)
   print, 'r02'
   print, r02
   r12=JHR_Coef(Density,Te,1,2,no_null=no_null)
   NHSaha1=NHSaha(Density,Te,1)
   NHSaha2=NHSaha(Density,Te,2)
   ok=where(N0 lt 1e32 and N0 gt 0 and r02 lt 1.0e32 and r12 lt 1.0e32 and NHSaha1 lt 1.0e32 and NHSaha2 lt 1.0e32,count)
   if count gt 0 then begin
      photons(ok)=A_Lyman(0)*(r02(ok)+r12(ok)*N0(ok)/NHSaha1(ok))*NHSaha2(ok)
      result(ok)=13.6057*0.75*photons(ok)*1.6e-19
   endif
   return,result
   end
