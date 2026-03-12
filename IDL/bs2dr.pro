;+
; BS2DR.pro
;
;   Calls public domain procedure B2VAL via FAST_B2VAL to evaluate the derivative of a two-dimensional 
; tensor-product spline given its the tensor-product B-spline representation. 
; Returns a vector of derivative evaluations at the positions specified by 
; vectors X and Y.
;
   Function BS2DR,IXDER,IYDER,X,Y,KXord,KYord,Xknot,Yknot,BSCoef
;
;______________________________________________________________________________
;   IXDER - integer, Order of the derivative in X-direction (input)
;   IYDER - integer, Order of the derivative in Y-direction (input)
;   X - fltarr, containing x coordinates at which the deriv is to be evaluated (input)
;   Y - fltarr, containing y coordinates at which the deriv is to be evaluated (input)
;   KXord- integer, Order of spline in X direction. 2=linear spline, 3=quadratic (input)
;   KYord- integer, Order of spline in Y direction. 2=linear spline, 3=quadratic (input)
;   Xknot - fltarr, X knot sequence (input) -> from B2SIN
;   Yknot - fltarr, Y knot sequence (input) -> from B2SIN
;   BSCoef- fltarr(nx*ny) of tensor product B-spline coefficients (input) -> from B2SIN
;______________________________________________________________________________
;-
   IXD=long(IXDER)
   IYD=long(IYDER)
   Xd=float(X)
   NV=n_elements(Xd)
   Yd=float(Y)
   Kx=long(KXord)
   Ky=long(KYord)
   Xk=float(Xknot)
   NXCOEF=n_elements(Xk)-Kx
   Yk=float(Yknot)
   NYCOEF=n_elements(Yk)-Ky
   BSC=float(BSCoef)
   Result=fltarr(NV)
   work=fltarr(3*max([KX,KY]) + KY)
   status=CALL_EXTERNAL('fast_b2val.so','fast_b2val_',NV,IXD,IYD,Xd,Yd,Kx,Ky,Xk,Yk,NXCOEF,NYCOEF,BSC,RESULT,WORK)
   return,Result
   end
