; 
; Tek.pro
;
; Sets graphics device to Tektronix 4100
pro tek
if !d.name eq 'PS' then begin
   print,'TEK => closing .ps file'
   device,/close
endif
set_plot,'tek'
device,/tek4100,colors=16,gin_chars=6
!p.font=-1
!x.thick=0.0
!y.thick=0.0
!p.thick=0.0
@setup_colors
end

