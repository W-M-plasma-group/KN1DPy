;+
; Integ_BL.pro
;
Function Integ_BL,X,Y,a=a,b=b,x_axis=x_axis,value_only=value_only,rev=rev,from_b_to_a=from_b_to_a
;
; If rev or from_b_to_a is set, integrate in direction from b to a. In this
; case, the last element of the result will be zero (instead of the
; first)
;-
ans=Y & ans(*)=0
n=n_elements(X)
if n ne n_elements(Y) then begin
  message,'number of elements of X and Y do not agree!',/info
  return,ans
endif
if n lt 2 then return,ans
rev_limits=X(n-1) lt X(0)
key_default,a,X(0)
key_default,b,X(n-1)
if rev_limits then begin
   _a=a < X(0)
   _b=b > X(n-1)
   ii=where(X ge _b and X le _a,count)
endif else begin
   _a=a > X(0)
   _b=b < X(n-1)
   ii=where(X ge _a and X le _b,count)
endelse
if count gt 0 then begin
   if X(ii(0)) eq _a then begin
      if X(ii(count-1)) eq _b then begin
         _X=[X(ii)]
         _Y=[Y(ii)]
      endif else begin
         _X=[X(ii),_b]
         _Y=[Y(ii),interpol(Y,X,_b)]
      endelse
   endif else begin
      if X(ii(count-1)) eq _b then begin
         _X=[_a,X(ii)]
         _Y=[interpol(Y,X,_a),Y(ii)]
      endif else begin
         _X=[_a,X(ii),_b]
         _Y=[interpol(Y,X,_a),Y(ii),interpol(Y,X,_b)]
      endelse
   endelse
endif else begin
   _X=[_a,_b]
   _Y=[interpol(Y,X,_a),interpol(Y,X,_b)]
endelse
x_axis=_X
if keyword_set(rev) or keyword_set(from_b_to_a) then begin
   _Y=reverse(_Y)
   _X=-reverse(_X)
endif
ans=_Y & ans(*)=0
for ii=1L,n_elements(_Y)-1 do ans(ii)=ans(ii-1)+(_X(ii)-_X(ii-1))*0.5*(_Y(ii)+_Y(ii-1))
if keyword_set(value_only) then begin
   return,ans(n_elements(ans)-1)
endif else begin
   if keyword_set(rev) or keyword_set(from_b_to_a) then return,reverse(ans)
   return,ans
endelse
end
;-
