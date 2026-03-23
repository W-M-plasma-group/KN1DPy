;
; File_Present.pro
;
function File_Present,file,errmsg=errmsg
   get_lun,unit
   openr,unit,file,error=err
   if (err ne 0) then begin
      id=0
      errmsg=!err_string
   endif else begin
      id=1
      errmsg=''
   endelse
   close,unit
   free_lun,unit
   return,id
   end
