;+
; Key_Default.pro
;
pro key_default,key,default
on_error,2
if type_of(default) eq 0 then message,'Default key value not defined'
if type_of(key) eq 0 then begin
   key=default
   return
endif
if (type_of(key) gt 0 and type_of(key) lt 6) and (type_of(default) gt 0 and type_of(default) lt 6) then return
if type_of(key) ne type_of(default) then key=default
return
end
;-
