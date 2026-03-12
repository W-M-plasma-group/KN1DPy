;+
; PrintCR.pro
;
;   Prints a string with carriage return
;
pro printCR,string,rev=rev,overwrite=overwrite
if n_params() lt 1 then string=''  
if keyword_set(rev) then revon
EEOL=''
if keyword_set(overwrite) then begin
  b=bytarr(3)
  b(0)='1B'XB
  b(1)=byte('[')
  b(2)=byte('K')
  EEOL=string(b) ; erase to end of line
endif
print,string+EEOL
if keyword_set(overwrite) then begin
  b=bytarr(3)
  b(0)='1B'XB
  b(1)=byte('[')
  b(2)=byte('A')
  a=string(b) ; up one line
  print,format='("'+a+'",$)'
endif
if keyword_set(rev) then revoff
return
end
