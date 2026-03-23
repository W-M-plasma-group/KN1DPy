;
; Psymbol.pro
;
   function psymbol,psym,color=color,nofill=nofill,thick=thick,yscale=yscale,xscale=xscale
   if keyword_set(color) then use_color=1 else use_color=0
   if keyword_set(nofill) then nofill=1 else nofill=0
   if keyword_set(thick) then begin
      thick=thick 
   endif else begin
      thick=!p.thick
   endelse
   key_default,yscale,1.0
   key_default,xscale,1.0
   if abs(psym) lt 8 then return,psym
   closed=0
   case abs(psym) of
      8:goto,none
      9:goto,cross
      10:goto,x
      11:goto,diamond
      12:goto,triangle
      13:goto,square
      14:goto,circle
      15:goto,invert_triangle
      16:goto,Left_triangle
      17:goto,Right_triangle
      18:goto,star
      19:goto,aster
      20:goto,dot
   else:goto,next
   endcase

next:
   closed=1
   case abs(psym) of
      21:goto,diamond
      22:goto,triangle
      23:goto,square
      24:goto,circle
      25:goto,invert_triangle
      26:goto,Left_triangle
      27:goto,Right_triangle
      28:goto,star
   else:goto,abort
   endcase

none:
; return system-defined psym = 0
   return,0

dot:
; return system-defined psym = 3
   return,sign(psym)*3

cross:
; return system-defined psym = 1
   return,sign(psym)*1
x:
; return system-defined psym = 7
   return,sign(psym)*7
aster:
; return system-defined psym = 2
   return,sign(psym)*2
;
; Define symbols with an area of 2 square units

Diamond:
   yy=[ 0.0,1.0,0.0,-1.0, 0.0]
   xx=[-1.0,0.0,1.0, 0.0,-1.0]
   goto,plotit

Triangle:
   height=sqrt(2.0)*3.0^0.25
   halfbase=2.0/height
   below=halfbase/sqrt(3.0)
   above=height-below
   yy=[-below,above,-below,-below]
   xx=[-halfbase,0,halfbase,-halfbase]
   goto,plotit

Square:
   r2=sqrt(2.0)/2.0
   yy=[-r2, r2, r2,-r2,-r2]
   xx=[-r2,-r2, r2, r2,-r2]
   goto,plotit

Circle:
   angle=2*!pi*findgen(49)/48 & radius=sqrt(2.0/!pi)
   xx=radius*cos(angle) & yy=radius*sin(angle)
   goto,plotit

Invert_Triangle:
   height=sqrt(2.0)*3.0^0.25
   halfbase=2.0/height
   below=halfbase/sqrt(3.0)
   above=height-below
   yy=[below,-above,below,below]
   xx=[-halfbase,0,halfbase,-halfbase]
   goto,plotit

Left_Triangle:
   height=sqrt(2.0)*3.0^0.25
   halfbase=2.0/height
   below=halfbase/sqrt(3.0)
   above=height-below
   xx=[below,-above,below,below]
   yy=[-halfbase,0,halfbase,-halfbase]
   goto,plotit

Right_Triangle:
   height=sqrt(2.0)*3.0^0.25
   halfbase=2.0/height
   below=halfbase/sqrt(3.0)
   above=height-below
   xx=[-below,above,-below,-below]
   yy=[-halfbase,0,halfbase,-halfbase]
   goto,plotit

Star:
   r1=sqrt(2/1.04068)
   r2=r1*cos(72*!dtor)/cos(36*!dtor)
   id=indgen(11)
   angle=2*!pi*id/10+0.5*!pi
   even=where(id mod 2 eq 0)
   odd=where(id mod 2 eq 1)
   xx=fltarr(11) & yy=fltarr(11)
   xx(even)=r1*cos(angle(even)) & xx(odd)=r2*cos(angle(odd))
   yy(even)=r1*sin(angle(even)) & yy(odd)=r2*sin(angle(odd))
   goto,plotit


plotit:
;
; Scale symbols to have an area of !pi square units
;
   xx=xx*sqrt(!pi/2.0)*xscale
   yy=yy*sqrt(!pi/2.0)*yscale
;
; set usersym
;
   if closed eq 1 and not nofill then begin
      npts=n_elements(xx)
      if use_color then begin
         usersym,xx,yy,fill=1,color=color,thick=thick
;         usersym,xx(0:npts-2),yy(0:npts-2),fill=1,color=color,thick=thick
      endif else begin
         usersym,xx,yy,fill=1,thick=thick
;         usersym,xx(0:npts-2),yy(0:npts-2),fill=1,thick=thick
      endelse
   endif else begin
      if use_color then begin
         usersym,xx,yy,fill=0,color=color,thick=thick
      endif else begin
         usersym,xx,yy,fill=0,thick=thick
      endelse
   endelse
   if psym lt 0 then return,-8 else return,8

abort:
   return,sign(psym)*1
   end

