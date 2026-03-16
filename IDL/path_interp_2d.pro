;+
;________________________________________________________________________________
; Path_Interp_2D.pro
;________________________________________________________________________________
; Interpolates an inputted 2D data array along a specified path or trajectory.
; Sets up call to IDL routine INTERPOLATE
;
Function Path_Interp_2D,P,PX,PY,X,Y,grid=grid
;
;________________________________________________________________________________
; Variables Passed:
;	P 	- Inputted 2D data array
;	PX 	- independent coordinate cooresponding to first dimension of P
;	PY 	- independent coordinate cooresponding to second dimension of P
;	X 	- inputted 1st coordinate of trajectory (vector) 
;	Y 	- inputted 2nd coordinate of trajectory (vector) 
; Returned:
;       if GRID is not set:
;	  array of same size as X containing P interpolated to coordinates X,Y
;       if GRID is set:
;         array of size (n_elements(X),n_elements(Y)) interpolated on a regular X,Y coordinate grid
;________________________________________________________________________________
;-
;
; Check for non-monotonic PX, PY
;
   dPX=PX-shift(PX,1)
   ii=where(dPX(1:*) lt 0,count)
   if count gt 0 then message,'ERROR in PATH_INTERP_2D => PX is non-monotonic!'
   dPY=PY-shift(PY,1)
   ii=where(dPY(1:*) lt 0,count)
   if count gt 0 then message,'ERROR in PATH_INTERP_2D => PY is non-monotonic!'
;
; compute coordinates normalized to array indices
;
   iPX=findgen(n_elements(PX))
   iPY=findgen(n_elements(PY))
   iX=interpol(iPX,float(PX),float(X))
   iY=interpol(iPY,float(PY),float(Y))
;
; Call interpolate
;
   return,interpolate(P,iX,iY,grid=grid)
   end
