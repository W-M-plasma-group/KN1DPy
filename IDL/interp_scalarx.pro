;+
; Interp_ScalarX.pro
;
;   Interpolates 'density' profiles used by Kinetic_Neutrals.pro,
; Kinetic_H2.pro, Kinetic_H2.pro, and other related procedures.
;
pro Interp_ScalarX,fa,Xa,fb,Xb,warn=warn,debug=debug
;________________________________________________________________________________
;  Input:
;     Input function 'a'
;	fa	- dblarr(nXa) density 
;       Xa	- fltarr(nXa)  - spatial coordinate
;
;    Desired space coordinates of Output function 'b'
;       Xb	- fltarr(nXb)  - spatial coordinate
;
;  Output:
;     Interpolated function 'b'
;	fb	- dblarr(nVrb,nVxb,nXb)
;
;  Keywords:
;     Input:
;	Warn	- float, acceptable truncation level.
;		  For interpolations outside the space set by (Xa), the values of fb are set to zero.
;		  This may not be acceptable. A test is performed on
;		  fb at the boundaries. If fb at the boundaries is greater
;		  than Warn times the maximum value of fb,
;		  a warning message is generated.
;________________________________________________________________________________
;-
;
   key_default,debug,0
   nxa=n_elements(xa)

   if n_elements(fa) ne nxa then begin 
      message,'Number of elements in fa and Xa do not agree!' 
   endif

   okk=where(Xb le max(Xa) and Xb ge min(Xa),nk)
   if nk lt 1 then message,'No values of Xb are within range of Xa'
   k0=okk(0) & k1=okk(nk-1)

   nxb=n_elements(xb)
   fb=dblarr(nXb)
;
; Call interpol
;
   fb(k0:k1)=interpol(fa,xa,xb(okk))
;
   if keyword_set(warn) then begin
;
; Test Boundaries
;
; k0 & k1
;
      big=max(abs(fb))
      if k0 gt 0 or k1 lt nXb-1 then begin
         if (k0 gt 0) and (abs(fb(k0)) gt warn*big) then begin
            message,'Non-zero value of fb detected at min(Xa) boundary',/info
         endif
         if (k1 lt nxb-1) and (abs(fb(k1)) gt warn*big) then begin
            message,'Non-zero value of fb detected at max(Xa) boundary',/info
         endif
      endif
   endif

   return
   end
