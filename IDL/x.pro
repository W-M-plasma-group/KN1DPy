;+
pro x,init=init,xsize=xsize,ysize=ysize,title=title,mag=mag,window=window,colors=colors,ypos=ypos,xpos=xpos,pixmap=pixmap,colortable=colortable,reuse=reuse
@/home/labombard/idl_lib/decw_display.include
mag_default=1.5
key_default,colortable,-1
key_default,init,0
key_default,reuse,0
if reuse then begin
   if type_of(window) eq 0 then message,'WINDOW is not defined!'
   if !d.name ne 'X' then begin
      set_plot,'X' 
   endif
   if xdefined(window) then begin
      xset,window,/refresh
      return
   endif
endif
init_done=n_elements(_active) eq 32
if not init_done then begin
   _window=0
   _colors=replicate(16,32)
   _mag=replicate(mag_default,32)
   _xsize=replicate(640*mag_default,32)
   _ysize=replicate(480*mag_default,32)
   _wtitle=replicate('X',32)
    for ii=0,31 do _wtitle(ii)=_wtitle(ii)+':'+sval(ii)
   _xpos=replicate(60+(_xsize(0)-640*mag_default),32)
   _ypos=replicate(680+(_ysize(0)-480*mag_default),32)
   _pixmap=intarr(32)
   _colortable=intarr(32)-1
   bp=!p & bp.color=1
   _bang_model={p:bp, x:!x, y:!y, z:!z}
   _bang=replicate(_bang_model,32)
   _last_window=_window
   _active=intarr(32)
endif
if type_of(window) eq 0 then begin
   window=_window
endif
_window=window
if type_of(colors) eq 0 then colors=_colors(_window)
if type_of(xsize) eq 0 then xsize=_xsize(_window)
if xsize eq 0 then xsize=640
if type_of(ysize) eq 0 then ysize=_ysize(_window)
if ysize eq 0 then ysize=480
if type_of(title) eq 0 then title=_wtitle(_window)
if type_of(xpos) eq 0 then xpos=_xpos(_window)
if type_of(ypos) eq 0 then ypos=_ypos(_window)
if type_of(pixmap) eq 0 then pixmap=_pixmap(_window)
if type_of(colortable) eq 0 then colortable=_colortable(_window)
if type_of(colortable) ne 0 then _colortable(_window)=colortable
_window=window
_colors(_window)=colors
if type_of(mag) ne 0 then begin
   xsize=640*mag & ysize=480*mag
endif else begin
   mag=mag_default
endelse
_xsize(_window)=xsize
_ysize(_window)=ysize
_mag(_window)=mag
_wtitle(_window)=title
_xpos(_window)=xpos
_ypos(_window)=ypos
_pixmap(_window)=pixmap
_active(_window)=1
if init then return
if _window ne _last_window then begin
  _bang(_last_window).p=!p
  _bang(_last_window).x=!x
  _bang(_last_window).y=!y
  _bang(_last_window).z=!z
  _last_window=_window
endif
if !d.name eq 'PS' then device,/close

if strlen(getenv('DISPLAY')) gt 0 then begin
   set_plot,'X'
   device,retain=2,decompose=0,true=24
   _w=window
   window,_w,colors=_colors(_w),xsize=_xsize(_w),ysize=_ysize(_w),title=_wtitle(_w),retain=2,xpos=_xpos(_w),ypos=_ypos(_w),pixmap=_pixmap(_w)
;   wdelete
   !p.font=-1
   !x.thick=0.0
   !y.thick=0.0
   !p.thick=0.0
endif
@/usr/local/cmod/codes/edge/db/plot/setup_colors
end
