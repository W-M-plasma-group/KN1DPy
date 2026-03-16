;+
; NHSaha.pro
;
; Evaluates the Saha equilibrium population density (m^-3)
; for atomic hydrogren level p. 
;
Function NHSaha,Density,Te,p
;
;________________________________________________________________________________
;  Density	- fltarr, electron density (=hydrogen ion density) (m^-3)
;  Te		- fltarr, electron temperature (eV
;  p		- intarr, hydrogen energy level, p=1 is ground state
;________________________________________________________________________________
;  Coding by B. LaBombard   6/29/99
;-
;
;  a=h^3/(2*!pi*m*k)^1.5
;  h=6.626e-34 !j-s
;  k=1.603e-19 !j/eV
;  m=0.91094e-30 kg
;  a=(h/(sqrt(2*!pi*m)*sqrt(k)))^3
;  a = {j^3 s^3 kg^-1.5 j^-1.5 eV^1.5)
;  a = {j^1.5 s^3 kg^-1.5 eV^1.5)
;  j=kg m^2 s^-2
;  a = {kg^1.5 m^3 s^-3 s^3 kg^-1.5 eV^1.5)
;  a = {m^3 eV^1.5)
;  a=3.310e-28
;
   if n_elements(Density) ne n_elements(Te) then message,'Number of elements of Density and Te are different!'
   if n_elements(p) ne 1 then message,'"p" must be a scalar'
   if p lt 0 then message,'"p" must greater than 0'
   Result=Density & Result(*)=1.0e32
   ok=where(Density lt 1.0e32 and Density gt 0.0 and Te lt 1.0e32 and Te gt 0.0,count)
   if count gt 0 then begin
      result(ok)=Density(ok)*(3.310E-28*Density(ok))*p*p*exp(13.6057/(p*p*Te(ok)))/Te(ok)^1.5
   endif
   return,result
   end
