;-------------------------------------------------------------
;+
; NAME:
;       SIGN
; PURPOSE:
;       Return the mathematical sign of the argument.
; CATEGORY:
; CALLING SEQUENCE:
;       s = sign(x)
; INPUTS:
;       x = value of array of values.     in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       s = sign of value(s).             out
; COMMON BLOCKS:
; NOTES:
;       Note:
;         s = -1 for x < 0
;         s =  0 for x = 0
;         s =  1 for x > 0
; MODIFICATION HISTORY:
;       R. Sterner, 7 May, 1986.
;       Johns Hopkins University Applied Physics Laboratory.
;       RES 15 Sep, 1989 --- converted to SUN.
;	RES 23 Sep, 1991 --- rewrote, reducing 11 lines of code to 1.
;
; Copyright (C) 1986, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	FUNCTION SIGN, X, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Return the mathematical sign of the argument.'
	  print,' s = sign(x)'
	  print,'   x = value of array of values.     in'
	  print,'   s = sign of value(s).             out'
	  print,' Note:'
	  print,'   s = -1 for x < 0'
	  print,'   s =  0 for x = 0'
	  print,'   s =  1 for x > 0'
	  return, -1
	endif
 
	return, fix(x gt 0.) - fix(x lt 0.)
	end
