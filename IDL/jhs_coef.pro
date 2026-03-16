;+
; JHS_Coef.pro
;
;  Evaluates the ionization rate coefficient, S (m^-3 s^-1), from  Johnson-Hinnov table 2 (MKS units).
;
Function JHS_Coef,Density,Te,create=create,no_null=no_null
;
;________________________________________________________________________________
; Input:
;  	Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
;  	Te	- fltarr, electron temperature (eV
;
; Keywords:
;	create	- if set, then create bi-cubic spline coefficients for
;		  interpolation of S and save them in the default save set. 
;	No_Null	- if set, then rather than generate a NULL value when Density and Te
;                 are outside the data range, compute the rate based on the min or max
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
;      restore,'jh_bscoef.dat'
      restore,'/home/labombard/edge/jh/jh_bscoef.dat'
   endif
;
; Evaluate S coefficients
;
   if n_elements(Density) ne n_elements(Te) then message,'Number of elements of Density and Te are different!'
   Result=Density & Result(*)=1.0e32
   LDensity=Alog(Density)
   LTe=Alog(Te)
   if keyword_set(No_Null) then begin
      LDensity=LDensity > (min(Dknot)+.001)
      LDensity=LDensity < (max(Dknot)-.001)
      LTe=LTe > (min(Tknot)+.001)
      LTe=LTe < (max(Tknot)-.001)
      count=n_elements(LDensity)
      ok=lindgen(count)
   endif else begin
      ok=where(LDensity le max(Dknot) and LDensity ge min(Dknot) and LTe le max(Tknot) and LTe ge min(Tknot),count)
   endelse
   if count gt 0 then Result(ok)=exp(BS2DR(0,0,LDensity(ok),LTe(ok),order,order,Dknot,Tknot,LogS_BSCoef))
   return,result
   end
