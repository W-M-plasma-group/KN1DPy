;+
;  Create_VrVxMesh.pro
;
; Sets up optimum Vr and Vx velocity space Mesh for Kinetic_Neutrals procedure
;_________________________________________________________________________________
pro Create_VrVxMesh,nv,Ti,vx,vr,Tnorm,E0=E0,ixE0=ixE0,irE0=irE0,Tmax=Tmax
;
; Input:
;	nv - integer, number of elements desired in vr mesh
;	Ti  - fltarr(*), Ti profile
;
;   Keyword:
;	E0  - fltarr, energy where a velocity bin is desired (optional)
;	Tmax- float, ignore Ti above this value
;
; Output:
;	vx  - dblarr(2*nv+1), x velocity grid
;	vr  - dblarr(nv), r velocity grid
;	Tnorm - float, optimum normalization temperature
;
;   Keyword:
;	ixE0 - intarr, returns array elements of vx corresponding to energy E0
;	irE0 - integer, returns array element of vr corresponding to energy E0
;_________________________________________________________________________________
;-
;
; Notes on algorithm
;	
;  Vmin=0.0
;
;  V=a*i + b*i^2
;
; Condition 1
;
; dVdi|0 = dVdi|imax * sqrt(min(Ti)/max(Ti))
; => a = (a+2*b*imax)* sqrt(min(Ti)/max(Ti))
;    a*(1-sqrt(min(Ti)/max(Ti)) = 2*b*imax*sqrt(min(Ti)/max(Ti))
;    a = 2*b*imax*sqrt(min(Ti)/max(Ti))/(1-sqrt(min(Ti)/max(Ti))
;
;    G== 2*imax*sqrt(min(Ti)/max(Ti))/(1-sqrt(min(Ti)/max(Ti))
;    a = G*b
;
; Condition 2
;    V|imax= Vmax = a*imax + b*imax^2
;
;  1 + 2
;
;     b*imax^2 = Vmax - a*imax
;     b*imax^2 + G*b*imax = Vmax
;     b = Vmax/(imax*(imax+G))    
;
   Key_Default,E0,[0.0]
   Key_Default,Tmax,0.0
   _Ti=Ti
   for k=0,n_elements(E0)-1 do if E0(k) gt 0 then _Ti=[_Ti,E0(k)]
   if Tmax gt 0 then begin
      ii=where(_Ti lt Tmax,count)
      if count gt 0 then _Ti=_Ti(ii)
   endif
   maxTi=max(_Ti)
   minTi=min(_Ti)
   Tnorm=1.0D0*total(_Ti)/n_elements(_Ti)
   vmax=3.5D0
   if maxTi-minTi le 0.1*maxTi then begin
      v=Dindgen(nv+1)*Vmax/nv
   endif else begin
      G=2*nv*sqrt(minTi/maxTi)/(1-sqrt(minTi/maxTi))
      b=Vmax/(nv*(nv+G))
      a=G*b
      v=a*dindgen(nv+1)+b*dindgen(nv+1)^2
   endelse
;
; Option: add velocity bins corresponding to E0
;
   v0=0.0
   for k=0,n_elements(E0)-1 do begin
      if E0(k) gt 0.0 then begin
         v0=sqrt(E0(k)/Tnorm)
         ii=where(v gt v0,count)
         if count gt 0 then begin
            v=[v(0:ii(0)-1),v0,v(ii(0):*)]
         endif else begin
            v=[v,v0]
         endelse
      endif
   endfor

   vr=v(1:*)
   vx=[-reverse(vr),vr]

   ixE0=where(abs(vx) eq v0,count)
   if count eq 1 then ixE0=ixE0(0)
   irE0=where(vr eq v0,count)
   if count eq 1 then irE0=irE0(0)

   return
   end
