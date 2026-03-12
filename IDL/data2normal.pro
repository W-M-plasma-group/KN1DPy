function data2normal,data,x=x,y=y
   if keyword_set(y) then begin
      scale=!Y.S(1)
      offset=!Y.S(0)
      type=!y.type
   endif else begin
      scale=!x.S(1)
      offset=!x.S(0)
      type=!x.type
   endelse
   d=data
   if type ne 0 then d=alog10(data)
   norm=scale*d+offset
   return,norm
   end
