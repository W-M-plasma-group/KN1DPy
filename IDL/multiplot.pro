;+
;________________________________________________________________________________
; multiplot.pro
;________________________________________________________________________________
; plots multiple data arrays versus a common 'time' axis
;
Pro Multiplot,X1,Y1,X2,Y2,X3,Y3,X4,Y4,X5,Y5,X6,Y6,X7,Y7,X8,Y8,$
              OX1=OX1,OY1=OY1,OX2=OX2,OY2=OY2,OX3=OX3,OY3=OY3,$
              OX4=OX4,OY4=OY4,OX5=OX5,OY5=OY5,OX6=OX6,OY6=OY6,$
              OX7=OX7,OY7=OY7,OX8=OX8,OY8=OY8,$
              title=title,xtitle=xtitle,xlegend=xlegend,ylegend=ylegend,$
              Rytitle=Rytitle,ytitle=ytitle,xrange=xrange,yrange=yrange,varlabel=varlabel,ovarlabel=ovarlabel,rvarlabel=rvarlabel,$
              psym=psym,opsym=opsym,rpsym=rpsym,Ryrange=Ryrange,$
              ylog=ylog,xlog=xlog,ystyle=ystyle,xstyle=xstyle,logy=logy,logx=logx,npts=npts,$
              color=color,ocolor=ocolor,rcolor=rcolor,ytop=ytop,ybottom=ybottom,xleft=xleft,xright=xright,$
	      noerase=noerase,dbfilter=dbfilter,xticks=xticks,yticks=yticks,ryticks=ryticks,$
	      xtickv=xtickv,ytickv=ytickv,logRy=logRy,xmark=xmark

common multiplotR,RX1,RY1,RX2,RY2,RX3,RY3,RX4,RY4,RX5,RY5,RX6,RY6,RX7,RY7,RX8,RY8
;________________________________________________________________________________
;
; Input Variables:
;       X1,Y1,  X2,Y2,  X3,Y3,... are fltarr(*)s of X,Y data to be plotted in  panel number 1,2,3,...
;		
; Input Keywords:
;	OX1,OY1, OX2,OY2, OX3,OY3,... are fltarr(*)s of X,Y data to be overplotted in panel number 1,2,3...
;	Title, XTitle, xrange, xstyle, ystyle,  ... same as in PLOT
;	ylog (or logy), xlog (or logx) are fltarr(*) arrays corresponding to panel *+1
;	yrange is a fltarr(2,*) giving y-range for panel number *+1
;	ytickv is a fltarr(*,*) giving y-tick values for panel number *+1
;	Rytitle,Ytitle, yticks, psym,opsym,color,ocolor are strarr(*) and fltarr(*) arrays corresponding to panel *+1
;	npts specifies the number of points to use to represent X*,Y* on a plot panel
;	ytop,ybottom,xleft,xright specify the boundaries of the plot (normal units)
;	noerase specifies no erase before plotting.
;       varlabel, ovarlabel - array of variable names to display on plots and oplots
;	dbfilter - if set then remove null values before plotting
;________________________________________________________________________________
;-
@decw_display.include
common multiplot,null
if type_of(logy) ne 0 then ylog=logy
if type_of(logx) ne 0 then xlog=logx
x,/init
Key_Default,title,''
Key_Default,xtitle,''
np=n_params()
if np lt 2 then message,'X and Y must be specified!'
if (np mod 2) eq 1 then message,'Must have an even number of arguments!'
np=np/2
key_default,RX1,[0]
key_default,RY1,[0]
ii=where(RX1 lt 1.0e32 and RY1 lt 1.0e32 and finite(RY1) and finite(RX1)) & RX1=RX1(ii) & RY1=RY1(ii)
key_default,RX2,[0]
key_default,RY2,[0]
ii=where(RX2 lt 1.0e32 and RY2 lt 1.0e32 and finite(RY2) and finite(RX2)) & RX2=RX2(ii) & RY2=RY2(ii)
key_default,RX3,[0]
key_default,RY3,[0]
ii=where(RX3 lt 1.0e32 and RY3 lt 1.0e32 and finite(RY3) and finite(RX3)) & RX3=RX3(ii) & RY3=RY3(ii)
key_default,RX4,[0]
key_default,RY4,[0]
ii=where(RX4 lt 1.0e32 and RY4 lt 1.0e32 and finite(RY4) and finite(RX4)) & RX4=RX4(ii) & RY4=RY4(ii)
key_default,RX5,[0]
key_default,RY5,[0]
ii=where(RX5 lt 1.0e32 and RY5 lt 1.0e32 and finite(RY5) and finite(RX5)) & RX5=RX5(ii) & RY5=RY5(ii)
key_default,RX6,[0]
key_default,RY6,[0]
ii=where(RX6 lt 1.0e32 and RY6 lt 1.0e32 and finite(RY6) and finite(RX6)) & RX6=RX6(ii) & RY6=RY6(ii)
key_default,RX7,[0]
key_default,RY7,[0]
ii=where(RX7 lt 1.0e32 and RY7 lt 1.0e32 and finite(RY7) and finite(RX7)) & RX7=RX7(ii) & RY7=RY7(ii)
key_default,RX8,[0]
key_default,RY8,[0]
ii=where(RX8 lt 1.0e32 and RY8 lt 1.0e32 and finite(RY8) and finite(RX8)) & RX8=RX8(ii) & RY8=RY8(ii)
Key_Default,psym,intarr(np)
Key_Default,rpsym,intarr(np)
Key_Default,opsym,intarr(np)
Key_Default,ytitle,replicate('',np)
Key_Default,Rytitle,replicate('',np)
if type_of(ylog) ne 0 then if n_elements(ylog) ne np then ylog=replicate(ylog(0),np)
Key_Default,ylog,replicate(0,np)
Key_Default,logRy,replicate(0,np)
Key_Default,ocolor,replicate(1,np)
Key_Default,rcolor,replicate(1,np)
Key_Default,ytop,.9
Key_Default,ybottom,.1
Key_Default,xleft,.13
Key_Default,xright,.92
key_default,noerase,0
key_default,dbfilter,0
key_default,xticks,0
key_default,xtickv,[0]
key_default,yticks,intarr(np)
key_default,Ryticks,intarr(np)
key_default,ytickv,fltarr(1,np)
key_default,xmark,[1.7e38]
_xmark=[xmark]
_xlegend=0.02
_ylegend=.8
if keyword_set(xlegend) then _xlegend=xlegend
if keyword_set(ylegend) then _ylegend=ylegend
!p.multi=0 & !p.multi(1)=1 & !p.multi(2)=np
xwindow=xright-xleft & ywindow=ytop-ybottom
dxw=xwindow/1 & dyw=ywindow/np
case np of
1: !p.charsize=1.0
2: !p.charsize=1.2
3: !p.charsize=1.4
4: !p.charsize=1.6
5: !p.charsize=1.8
6: !p.charsize=2.0
7: !p.charsize=2.0
8: !p.charsize=2.0
else: !p.charsize=1.0
endcase
!p.charsize=!p.charsize*_mag(_window)
key_default,color,replicate(1,np)
if type_of(xrange) eq 0 then begin
   xrange=[min(X1),max(X1(where(X1 lt 1.0e32)))]
   inrange=where(X1 lt 1.0e32 and Y1 lt 1.0e32,count)
   if count gt 1 then xrange=[min(X1(inrange)),max(X1(inrange))]
   for i=1,np-1 do begin
      cmd='inrange=where(X'+sval(i+1)+' lt 1.0e32 and Y'+sval(i+1)+' lt 1.0e32,count)'
      status=execute(cmd)
      if count gt 1 then begin
         cmd='xrange=[min([xrange,X'+sval(i+1)+'(inrange)]),max([xrange,X'+sval(i+1)+'(inrange)])]'
         status=execute(cmd)
      endif
   endfor
endif
nx=0
key_default,npts,4*1024
xticklen=0
if !x.ticklen eq 0 then begin
   xticklen=!x.ticklen
   !x.ticklen=(0.02+np*0.01)/(ytop-ybottom)
endif
yticklen=0
if !y.ticklen eq 0 then begin
   yticklen=!y.ticklen
   !y.ticklen=0.015
endif
_xlegend=_xlegend*dxw
_ylegend=_ylegend*dyw
if not keyword_set(yrange) then begin
   yrange=fltarr(2,np)
   for i=0,np-1 do begin
      cmd='inrange=where(X'+sval(i+1)+' ge xrange(0) and X'+sval(i+1)+' le xrange(1) '+$
          'and Y'+sval(i+1)+' lt 1.0e32,count)'
      status=execute(cmd)
      if count gt 0 then begin
         status=execute('yrange(*,i)=[min([Y'+sval(i+1)+'(inrange)]),max([Y'+sval(i+1)+'(inrange)])]')
         ix=0 & iy=0
         status=execute('ix=type_of(OX'+sval(i+1)+')')
         status=execute('iy=type_of(OY'+sval(i+1)+')')
         if ix ne 0 and iy ne 0 then begin
            cmd='inrange=where(OX'+sval(i+1)+' ge xrange(0) and OX'+sval(i+1)+' le xrange(1) '+$
                'and OY'+sval(i+1)+' lt 1.0e32,count)'
            status=execute(cmd)
            if count gt 0 then begin
               status=execute('yrange(*,i)=[min([yrange(*,i),OY'+sval(i+1)+'(inrange)]),max([yrange(*,i),OY'+sval(i+1)+'(inrange)])]')
            endif
         endif
      endif
   endfor
endif
for i=0,np-1 do begin
   ny=np-1-i
   if i eq np-1 then begin
      xcharsize=1.0 
   endif else begin
      xcharsize=1e-3
   endelse
   if i eq 0 then begin
      pTitle=Title
   endif else begin
      pTitle=''
   endelse
   !p.position=[xleft+dxw*nx,ybottom+dyw*ny,xleft+dxw*(nx+1),ybottom+dyw*(ny+1)]
   
   xleg=!p.position(0)+_xlegend
   oxleg=!p.position(2)-_xlegend
   yleg=!p.position(1)+_ylegend
   if dbfilter then begin
      status=execute('dbfilter,X'+sval(i+1)+',Y'+sval(i+1)+',X'+sval(i+1)+',Y'+sval(i+1))
   endif
   ix=0 & iy=0
   status=execute('ix=n_elements(RX'+sval(i+1)+')')
   status=execute('iy=n_elements(RY'+sval(i+1)+')')
   if type_of(ystyle) eq 0 then _ystyle=0 else _ystyle=ystyle
   if ix gt 1 and iy gt 1 then begin
      _ystyle=_ystyle or 8
   endif
   if yrange(1,i) gt yrange(0,i) then begin
      status=execute('plot,comb(X'+sval(i+1)+',npts=npts),comb(Y'+sval(i+1)+',npts=npts),Title=pTitle,XTitle=XTitle,'+$
			'YTitle=YTitle(i),yrange=yrange(*,i),xrange=xrange,psym=psymbol(psym(i)),xcharsize=xcharsize,'+$
                        'ylog=ylog(i),xlog=xlog,ystyle=_ystyle,xstyle=xstyle,/nodata,noerase=noerase,xticks=xticks,'+$
                        'yticks=yticks(i),xtickv=xtickv,ytickv=ytickv(*,i)')
   endif else begin
      status=execute('plot,comb(X'+sval(i+1)+',npts=npts),comb(Y'+sval(i+1)+',npts=npts),Title=pTitle,XTitle=XTitle,'+$
			'YTitle=YTitle(i),xrange=xrange,psym=psymbol(psym(i)),xcharsize=xcharsize,'+$
                        'ylog=ylog(i),xlog=xlog,ystyle=_ystyle,xstyle=xstyle,/nodata,noerase=noerase,xticks=xticks,'+$
                        'yticks=yticks(i),xtickv=xtickv,ytickv=ytickv(*,i)')
   endelse     
   cmd='ii=where(X'+sval(i+1)+' lt 1.0e32 and finite(X'+sval(i+1)+')'+$
                           ' and Y'+sval(i+1)+' lt 1.0e32 and finite(Y'+sval(i+1)+'),count)'
   status=execute(cmd)
   if count gt 0 then begin
      cmd='oplot,comb(X'+sval(i+1)+'(ii),npts=npts),comb(Y'+sval(i+1)+'(ii),npts=npts),psym=psymbol(psym(i)),color=color(i)'
      status=execute(cmd)
   endif
   if keyword_set(varlabel) then xyouts,xleg,yleg,/norm,align=0.0,varlabel(i),charsize=1.0*_mag(_window),color=color(i)
   ix=0 & iy=0
   status=execute('ix=type_of(OX'+sval(i+1)+')')
   status=execute('iy=type_of(OY'+sval(i+1)+')')
   if ix ne 0 and iy ne 0 then begin
      status=execute('oplot,comb(OX'+sval(i+1)+',npts=npts),comb(OY'+sval(i+1)+',npts=npts),psym=psymbol(opsym(i)),color=ocolor(i)')
      if keyword_set(ovarlabel) then xyouts,xleg,yleg-.05,/norm,align=0.0,ovarlabel(i),charsize=1.0*_mag(_window),color=ocolor(i)
   endif
   ix=0 & iy=0
   status=execute('ix=n_elements(RX'+sval(i+1)+')')
   status=execute('iy=n_elements(RY'+sval(i+1)+')')
   if ix gt 1 and iy gt 1 then begin
      RYvar='RY'+sval(i+1)

      cmd='_RYrange=[min('+RYvar+'(where(finite('+RYvar+')))),max('+RYvar+'(where(finite('+RYvar+') and '+RYvar+' lt 1.0e32)))]'
      status=execute(cmd)
      if keyword_set(ryrange) then _ryrange=reform(ryrange(*,i))
      if type_of(ystyle) ne 0 then begin
         status=execute('axis,yaxis=1,ystyle=ystyle,yrange=_ryrange,ytitle=RYtitle(i),yticks=Ryticks(i),ylog=logRy(i),/save')
      endif else begin
         status=execute('axis,yaxis=1,yrange=_ryrange,ytitle=RYtitle(i),yticks=Ryticks(i),ylog=logRy(i),/save')
      endelse
      status=execute('oplot,comb(RX'+sval(i+1)+',npts=npts),comb(RY'+sval(i+1)+',npts=npts),psym=psymbol(rpsym(i)),color=Rcolor(i)')
      if keyword_set(Rvarlabel) then xyouts,oxleg,yleg,/norm,align=1.0,rvarlabel(i),charsize=1.0*_mag(_window),color=rcolor(i)
   endif
   ii=where(_xmark lt 1.0e32,count)
   if count gt 0 then begin
      for k=0,count-1 do begin
         yr=!y.crange
         if ylog(i) then yr=10^yr
         oplot,[_xmark(k),_xmark(k)],yr,linestyle=5
      endfor
   endif
endfor
!p.charsize=0.0
!p.position=0.0
!p.multi=0
!x.ticklen=xticklen
!y.ticklen=yticklen
RX1=[0]
RY1=[0]
RX2=[0]
RY2=[0]
RX3=[0]
RY3=[0]
RX4=[0]
RY4=[0]
RX5=[0]
RY5=[0]
RX6=[0]
RY6=[0]
RX7=[0]
RY7=[0]
RX8=[0]
RY8=[0]
end
