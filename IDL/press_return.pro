;
; press_return
;
; B. LaBombard 1994
;
pro press_return,prompt=prompt,again=again,psprompt=psprompt,landscape=landscape,psfont=psfont,full=full
common press_return,file
On_Error,2
if not keyword_set(prompt) then Prompt='(Press any key to continue)'
if keyword_set(again) then again=1 else again=0
if keyword_set(psprompt) then psprompt=1 else psprompt=0
if psprompt then Prompt='(Press any key to continue, P=enter postscript filename)'
if psprompt and again then begin
   print,'Plotted to file: '+file
   device,/close
   tek
endif
print,prompt
a=get_kbrd(1)
if strupcase(a) eq 'S' then begin
   msg='Stopping program by request'
   message,msg,level=-1
endif
if strupcase(a) eq 'P' and psprompt then begin
   if type_of(file) ne 0 then print,'Previous postscript file was: '+file
   print,format='("Enter postscript file name",$)'
   file=''
   on_ioerror,abort
   read,file
   on_ioerror,null
   ps,landscape=landscape,/color,file=file,psfont=psfont,full=full
   again=1
   return
endif
print,' '
Abort:
   again=0
   on_ioerror,null
return
end
