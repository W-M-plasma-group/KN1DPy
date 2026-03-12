;+
; PrintNCR,pro
;
;   Prints a string with no carriage return
;
pro printNCR,string,rev=rev
if keyword_set(rev) then revon
print,format='("'+string+'",$)'
if keyword_set(rev) then revoff
return
end
