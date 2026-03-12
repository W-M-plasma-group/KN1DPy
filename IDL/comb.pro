;+
; Comb.pro
;
Function Comb,data,npts=npts
if not keyword_set(npts) then npts=2000L
if n_elements(data) gt npts then begin
   ii=(findgen(npts)*(n_elements(data)-1))/(npts-1)
   return,data(ii)
endif else begin
   return,data
endelse
end
;-
