;+
;________________________________________________________________________________
; Type_of.pro
;________________________________________________________________________________
; Tests argument to see if variable has been defined.
; Retuns the following integer based on data type:
;	0 -  undefined
;	1 - Byte
;	2 - Integer
;	3 - Longword Integer
;	4 - Floating point
;	5 - Double Precision
;	6 - Complex Floating
;	7 - String
;	8 - Structure
Function Type_of,arg,nDim=nDim,D1=D1,D2=D2,D3=D3,D4=D4,D5=D5,D6=D6,D7=D7,type=type
;-
_type=['UNDEFINED','BYTE','INTEGER','LONG','FLOAT','DOUBLE','COMPLEX','STRING','STRUCTURE']
s=size(arg)
n=n_elements(s)
id=s(n-2)
nDim=s(0)
D1=0 & D2=0 & D3=0 & D4=0 & D5=0 & D6=0 & D7=0
if nDim gt 0 then D1=S(1)
if nDim gt 1 then D2=S(2)
if nDim gt 2 then D3=S(3)
if nDim gt 3 then D4=S(4)
if nDim gt 4 then D5=S(5)
if nDim gt 5 then D6=S(6)
if nDim gt 6 then D7=S(7)
if id gt -1 and id lt 9 then begin
   type=_type(id)
endif else begin
   type='type not decoded'
endelse
return,id
end
