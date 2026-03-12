;
; sval.pro
;
;+
; function sval,val,len=len
;-
function sval,val,len=len,rfill0=rfill0,rfill32=rfill32
s=strtrim(string(val),2)
if keyword_set(len) then begin
   s=STRMID(s,0,len)
   return,s
endif
if keyword_set(rfill0) then begin
   s=rightfill(s,rfill0,'0')
   return,s
endif
if keyword_set(rfill32) then begin
   s=rightfill(s,rfill32,' ')
   return,s
endif
return,s
end
;
