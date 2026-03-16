pro color_set,r,g,b,colors

if !d.name eq 'PS' or !d.name eq 'PRINTER' then begin
r = [255,000,200,000,000,220,240,000,255,100,120,255,200,180,255,255,200,220,80,0,255,100,000,000,110,000,050]
g = [255,000,000,000,150,220,160,240,000,000,120,200,255,200,255,200,255,220,80,0,255,000,000,100,110,120,000]
b = [255,000,000,200,000,000,000,240,255,160,120,200,200,255,180,255,255,220,80,0,255,000,110,000,000,120,080]
endif else begin
r = [000,255,255,000,000,255,255,000,255,120,140,255,200,180,255,255,200,220,80,0,255,100,000,000,110,000,050]
g = [000,255,000,000,255,255,180,255,000,000,140,200,255,200,255,200,255,220,80,0,255,000,000,100,110,120,000]
b = [000,255,000,255,000,000,000,255,255,180,140,200,200,255,180,255,255,220,80,0,255,000,110,000,000,120,080]
endelse

colors = ['background', 'foreground', 'red', 'blue','green', 'yellow', $
          'orange','cyan', 'magenta', 'purple', 'gray','pale_red',$
          'pale_green', 'pale_blue','pale_yellow','pale_purple',$
          'pale_cyan','pale_gray','dark_gray','true_black', 'true_white', $
          'dark_red', 'dark_blue','dark_green', 'dark_yellow',$
          'dark_cyan', 'dark_purple']

tvlct,r,g,b
;black = 0
;white = 1
;red = 2
;blue = 3
;green = 4
;yellow = 5
;orange = 6
;cyan = 7
;magenta = 8
;purple = 9
;gray = 10
;pale_red = 11
;pale_green = 12
;pale_blue = 13
;pale_yellow = 14
;pale_purple = 15
;pale_cyan = 16
;pale_gray = 17
;dark_gray = 18
;true_black = 19
;true_white = 20
;dark_red = 21
;dark_blue = 22
;dark_green = 23
;dark_yellow = 24
;dark_cyan = 25
;dark_purple = 26

return
end

;;;;;;;

function color,color_name

color_set,r,g,b,colors

dum = size(color_name)
if dum(0) eq 0 and dum(1) eq 2 then $
  s_index = color_name $
  else if dum(0) eq 0 and dum(1) eq 7 then begin
    s_index = where(colors eq color_name,count)
    if count ne 1 then begin
    	print,'error specifying color or index!'
    	return,16777215
    endif else $
      s_index=s_index(0)
  endif else print,'error specifying color or index!!'

;code below is not needed when  device,decompose=0
;if !d.n_colors gt 256 then $
;	index = r(s_index) + long(256)*g(s_index) + long(256)*256*b(s_index) else $
	    index = s_index

return,index
end
