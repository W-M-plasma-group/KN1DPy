;+
; Interp_fVrVxX.pro
;
;   Interpolates distribution functions used by Kinetic_Neutrals.pro,
; Kinetic_H.pro, Kinetic_H2.pro, and other related procedures.
;
pro Interp_fVrVxX,fa,Vra,Vxa,Xa,Tnorma,fb,Vrb,Vxb,Xb,Tnormb,warn=warn,debug=debug,correct=correct
;________________________________________________________________________________
;  Input:
;     Input Distribution function 'a'
;	fa	- dblarr(nVra,nVxa,nXa) distribution function
;	Vra	- fltarr(nVra) - radial velocity
;	Vxa	- fltarr(nVxa) - axial velocity
;       Xa	- fltarr(nXa)  - spatial coordinate
;       Tnorma	- float,  Normalization temperature for Vra & Vxa
;
;    Desired phase space coordinates of Output Distribution function 'b'
;	Vrb	- fltarr(nVrb) - radial velocity
;	Vxb	- fltarr(nVxb) - axial velocity
;       Xb	- fltarr(nXb)  - spatial coordinate
;       Tnormb	- float,  Normalization temperature for Vrb & Vxb
;
;  Output:
;     Interpolated Distribution function 'b'
;	fb	- dblarr(nVrb,nVxb,nXb) distribution function
;	          fb is scaled if necessary to make its
;	          digital integral over all velocity space
;		  equal to that of fa.
;
;  Keywords:
;     Input:
;	Warn	- float, acceptable truncation level.
;		  For interpolations outside the phase space set by
;		  (Vra, Vxa, Xa), the values of fb are set to zero.
;		  This may not be acceptable. A test is performed on
;		  fb at the boundaries. If fb at the boundaries is greater
;		  than Warn times the maximum value of fb,
;		  a warning message is generated.
;________________________________________________________________________________
;-
;
   prompt='INTERP_FVRVXX => '
   common INTERP_FVRVXX_internal1,vra1,vxa1,Tnorma1,vrb1,vxb1,Tnormb1,weight1
   common INTERP_FVRVXX_internal2,vra2,vxa2,Tnorma2,vrb2,vxb2,Tnormb2,weight2

   key_default,debug,0
   key_default,correct,1
   nvra=n_elements(vra)
   nvxa=n_elements(vxa)
   nxa=n_elements(xa)
   mH=1.6726231D-27
   q=1.602177D-19	
   mu=1

   fV=sqrt(Tnormb/Tnorma)
;
; Compute Vtha, Vtha2, Vthb and Vthb2
;
   Vtha=sqrt(2*q*Tnorma/(mu*mH))
   Vtha2=vtha*vtha
   Vthb=sqrt(2*q*Tnormb/(mu*mH))
   Vthb2=vthb*vthb

   if n_elements(fa(*,0,0)) ne nvra then begin
      message,'Number of elements in fa(*,0,0) and Vra do not agree!' 
   endif
   if n_elements(fa(0,*,0)) ne nvxa then begin 
      message,'Number of elements in fa(0,*,0) and Vxa do not agree!' 
   endif
   if n_elements(fa(0,0,*)) ne nxa then begin 
      message,'Number of elements in fa(0,0,*) and Xa do not agree!' 
   endif

   oki=where(fV*Vrb le max(Vra) and fV*Vrb ge min(Vra),ni)
   if ni lt 1 then message,'No values of Vrb are within range of Vra'
   i0=oki(0) & i1=oki(ni-1)

   okj=where(fV*Vxb le max(Vxa) and fV*Vxb ge min(Vxa),nj)
   if nj lt 1 then message,'No values of Vxb are within range of Vxa'
   j0=okj(0) & j1=okj(nj-1)

   okk=where(Xb le max(Xa) and Xb ge min(Xa),nk)
   if nk lt 1 then message,'No values of Xb are within range of Xa'
   k0=okk(0) & k1=okk(nk-1)

   nvrb=n_elements(vrb)
   nvxb=n_elements(vxb)
   nxb=n_elements(xb)
   fb=dblarr(nVrb,nVxb,nXb)
   Make_dVr_dVx,vra,vxa,Vr2pidVra,VrVr4pidVra,dVxa,vrL=vraL,vrR=vraR,vxL=vxaL,vxR=vxaR,Vr2Vx2=Vra2Vxa2

   Make_dVr_dVx,vrb,vxb,Vr2pidVrb,VrVr4pidVrb,dVxb,vrL=vrbL,vrR=vrbR,vxL=vxbL,vxR=vxbR,$
                Vol=Vol,Vth_DeltaVx=Vth_DVx,Vx_DeltaVx=Vx_DVx,Vr_DeltaVr=Vr_DVr,Vr2Vx2=Vrb2Vxb2,$
                jpa=jpa,jpb=jpb,jna=jna,jnb=jnb

;
; Determine if Weight was already computed by checking vra_s,vxa_s,Tnorma_s,vrb_s,vxb_s,Tnormb_s
; for cases 1 and 2
;



   w1_active=0
   w1_match=0
   if type_of(vra1) ne 0 then begin
      w1_active=1
      test=0
      ii=where(vra1 ne vra,count) & test=test+count
      ii=where(vxa1 ne vxa,count) & test=test+count
      ii=where(Tnorma1 ne Tnorma,count) & test=test+count
      ii=where(vrb1 ne vrb,count) & test=test+count
      ii=where(vxb1 ne vxb,count) & test=test+count
      ii=where(Tnormb1 ne Tnormb,count) & test=test+count
      if test le 0 then w1_match=1
   endif
   w2_active=0
   w2_match=0
   if type_of(vra2) ne 0 then begin
      w2_active=1
      test=0
      ii=where(vra2 ne vra,count) & test=test+count
      ii=where(vxa2 ne vxa,count) & test=test+count
      ii=where(Tnorma2 ne Tnorma,count) & test=test+count
      ii=where(vrb2 ne vrb,count) & test=test+count
      ii=where(vxb2 ne vxb,count) & test=test+count
      ii=where(Tnormb2 ne Tnormb,count) & test=test+count
      if test le 0 then w2_match=1
   endif
   w_new=0
   if w1_match or w2_match  then begin
      if w1_match then begin
         weight=weight1
         if debug then print,prompt+'using Weight1'
      endif
      if w2_match then begin
         if debug then print,prompt+'using Weight2'
         weight=weight2
      endif
   endif else begin
; 
; If not then compute Weight for this combination of Vr and Vx
;
; Determine Left and Right limits on 'cells' for Vra, Vxa, Vrb, Vxb
;
      if debug then print,prompt+'computing new Weight'
      w_new=1
;
; Set area contributions to Weight array
;
      _weight=dblarr(nvrb,nvxb,nvra,nvxa)
      weight=dblarr(nvrb*nvxb,nvra*nvxa)
      for ib=0,nvrb-1 do begin   
         for jb=0,nvxb-1 do begin   
            for ia=0,nvra-1 do begin   
               vraMin=max([fv*vrbL(ib),vraL(ia)])
               vraMax=min([fv*vrbR(ib),vraR(ia)])
               for ja=0,nvxa-1 do begin   
                  vxaMin=max([fv*vxbL(jb),vxaL(ja)])
                  vxaMax=min([fv*vxbR(jb),vxaR(ja)])
                  if (vraMax gt vraMin) and (vxaMax gt vxaMin) then begin
                     _weight(ib,jb,ia,ja)=2*!pi*(vraMax^2-vraMin^2)*(vxaMax-vxaMin)/(Vr2pidVrb(ib)*dVxb(jb))
                  endif
               endfor
            endfor
         endfor
      endfor
      weight(*)=_weight
   endelse

   fb_xa=dblarr(nvrb*nvxb,nxa)
;
; Determine fb_xa from weight array
;
   _fa=dblarr(nvra*nvxa,nxa)
   _fa(*)=fa
   fb_xa=weight#_fa

;
; Compute _Wxa and _Ea - these are the desired moments of fb, but on the xa grid
;
   na=dblarr(nXa)
   _Wxa=dblarr(nXa)
   _Ea=dblarr(nXa)
   for k=0,nxa-1 do begin
      na(k)=total(Vr2pidVra*(fa(*,*,k)#dVxa))
      if na(k) gt 0 then begin
         _Wxa(k)=sqrt(Tnorma)*total(Vr2pidVra*(fa(*,*,k)#(Vxa*dVxa)))/na(k)
         _Ea(k)=Tnorma*total(Vr2pidVra*((vra2vxa2*fa(*,*,k))#dVxa))/na(k)
      endif
   endfor
;
; Interpolate in x to get fb from fb_xa and to get Wxa, Ea from _Wva, _Ea
;
   Wxa=dblarr(nXb)
   Ea=dblarr(nXb)
   for k=k0,k1 do begin
      kL=locate(xa,xb(k)) > 0
      kR=(kL+1) < (n_elements(xa)-1)
      kL = kL < (kR-1)
      f=(xb(k)-xa(kL))/(xa(kR)-xa(kL))
      fb(*,*,k)=fb_xa(*,kL)+(fb_xa(*,kR)-fb_xa(*,kL))*f
      Wxa(k)=_Wxa(kL)+(_Wxa(kR)-_Wxa(kL))*f
      Ea(k)=_Ea(kL)+(_Ea(kR)-_Ea(kL))*f
   endfor

;
; Correct fb so that it has the same Wx and E moments as fa
;
   if correct then begin
;
; Process each spatial location
;
      AN=dblarr(nvrb,nvxb,2)
      BN=dblarr(nvrb,nvxb,2)
      sgn=[1,-1]
      for k=0,nxb-1 do begin
         allow_neg=0
;
; Compute nb, Wxb, and Eb - these are the current moments of fb
;
         nb=total(Vr2pidVrb*(fb(*,*,k)#dVxb))

         if nb gt 0 then begin
;
; Entry point for iteration
correct:
            nb=total(Vr2pidVrb*(fb(*,*,k)#dVxb))
            Wxb=sqrt(Tnormb)*total(Vr2pidVrb*(fb(*,*,k)#(Vxb*dVxb)))/nb
            Eb=Tnormb*total(Vr2pidVrb*((vrb2vxb2*fb(*,*,k))#dVxb))/nb
;
; Compute Nij from fb, padded with zeros
;
            Nij=dblarr(nvrb+2,nvxb+2)
            Nij(1:nvrb,1:nvxb)=fb(*,*,k)*vol/nb
;
; Set Cutoff and remove Nij very close to zero
;
            cutoff=1.0e-6*max(Nij)
            ii=where(abs(Nij) lt cutoff and abs(Nij) gt 0.0,count)
            if count gt 0 then begin
;if debug then print,'These Nij set to zero:',Nij(ii)
               Nij(ii)=0.0
            endif
            ii=where(Nij(2,*) gt 0.0,count)
            if count lt 1 then begin
;if debug then print,'All Nij with i=1 are equal to zero. Allowing negative values here.'
               allow_neg=1
            endif

	    Nijp1_vx_Dvx=shift(Nij*vx_Dvx,0,-1)
	    Nij_vx_Dvx  =Nij*vx_Dvx
	    Nijm1_vx_Dvx=shift(Nij*vx_Dvx,0,1)
            Nip1j_vr_Dvr=shift(Nij*vr_Dvr,-1,0)
            Nij_vr_Dvr  =Nij*vr_Dvr
            Nim1j_vr_Dvr=shift(Nij*vr_Dvr,1,0)
;
; Compute Ap, Am, Bp, and Bm (0=p 1=m)
;
            _AN=shift(Nij*vth_Dvx,0,1) - Nij*vth_Dvx
            AN(*,*,0)=_AN(1:nvrb,1:nvxb)
            _AN=-shift(Nij*vth_Dvx,0,-1) + Nij*vth_Dvx
            AN(*,*,1)=_AN(1:nvrb,1:nvxb)

            BN(*,jpa+1:jpb,0)=Nijm1_vx_Dvx(1:nvrb,jpa+2:jpb+1)-Nij_vx_Dvx(1:nvrb,jpa+2:jpb+1)
            BN(*,jpa,0)=-Nij_vx_Dvx(1:nvrb,jpa+1)
            BN(*,jnb,0)=Nij_vx_Dvx(1:nvrb,jnb+1)
            BN(*,jna:jnb-1,0)=-Nijp1_vx_Dvx(1:nvrb,jna+1:jnb)+Nij_vx_Dvx(1:nvrb,jna+1:jnb)
            BN(*,*,0)=BN(*,*,0) + Nim1j_vr_Dvr(1:nvrb,1:nvxb)-Nij_vr_Dvr(1:nvrb,1:nvxb)

            BN(*,jpa+1:jpb,1)=-Nijp1_vx_Dvx(1:nvrb,jpa+2:jpb+1)+Nij_vx_Dvx(1:nvrb,jpa+2:jpb+1)
            BN(*,jpa,1)=-Nijp1_vx_Dvx(1:nvrb,jpa+1)
            BN(*,jnb,1)=Nijm1_vx_Dvx(1:nvrb,jnb+1)
            BN(*,jna:jnb-1,1)=Nijm1_vx_Dvx(1:nvrb,jna+1:jnb)-Nij_vx_Dvx(1:nvrb,jna+1:jnb)
            BN(1:nvrb-1,*,1)=BN(1:nvrb-1,*,1) - Nip1j_vr_Dvr(2:nvrb,1:nvxb)+Nij_vr_Dvr(2:nvrb,1:nvxb)
            BN(0,*,1)=BN(0,*,1) - Nip1j_vr_Dvr(1,1:nvxb)
;
; If negative values for Nij must be allowed, then add postive particles to i=0
; and negative particles to i=1 (beta is negative here)
;
            if allow_neg then begin
               BN(0,*,1)=BN(0,*,1) - Nij_vr_Dvr(1,1:nvxb)
               BN(1,*,1)=BN(1,*,1) + Nij_vr_Dvr(1,1:nvxb)
            endif
;
; Remove padded zeros in Nij
;         
            Nij=Nij(1:nvrb,1:nvxb)
;
; Cycle through 4 possibilies of sign(alpha),sign(beta)
;
            TB1=fltarr(2) & TB2=fltarr(2)
            for ia=0,1 do begin
;
; Compute TA1, TA2
;
               TA1=sqrt(Tnormb)*total(AN(*,*,ia)#Vxb)
               TA2=Tnormb*total(vrb2vxb2*AN(*,*,ia))
               for ib=0,1 do begin
;
; Compute TB1, TB2
;
                  if TB1(ib) eq 0 then TB1(ib)=sqrt(Tnormb)*total(BN(*,*,ib)#Vxb)
                  if TB2(ib) eq 0 then TB2(ib)=Tnormb*total(vrb2vxb2*BN(*,*,ib))
                  denom=TA2*TB1(ib)-TA1*TB2(ib)
                  beta=0.0
                  alpha=0.0
                  if denom ne 0.0 and TA1 ne 0.0 then begin
                     beta=(TA2*(Wxa(k)-Wxb)-TA1*(Ea(k)-Eb))/denom
                     alpha=(Wxa(k)-Wxb-TB1(ib)*Beta)/TA1
                  endif
;if debug then print,'Alpha: '+sval(alpha)+'    Beta:'+sval(beta)
                  if alpha*sgn(ia) gt 0.0 and beta*sgn(ib) gt 0.0 then goto,alpha_beta
               endfor
            endfor
Alpha_Beta:
            RHS=AN(*,*,ia)*alpha+BN(*,*,ib)*beta
;
; Are there locations where Nij = 0.0 and RHS is negative?
;
            ii=where(Nij eq 0.0 and RHS lt 0.0,count)
;if count gt 0 and debug then begin
;   print,'There are locations where Nij = 0.0 and RHS is negative'
;   press_return
;endif
;
; Determine limiting scale factor                     
;
            s=1.0
            if not allow_neg then begin
               ii=where(Nij ne 0.0,count)
               if count gt 0 then s=min([1.0/max(-RHS(ii)/Nij(ii)),1.0])
            endif
;if debug then print,k,s
;if debug then plot,fb(0,*,k)
            fb(*,*,k)=nb*(Nij+s*RHS)/vol
;if debug then oplot,fb(0,*,k),color=2
;if debug then press_return
            if s lt 1.0 then goto,correct
         endif
      endfor
   endif
;
   if keyword_set(warn) then begin
;
; Test Boundaries
;
; i0 & i1
;
      big=max(fb)
      i0_error=0
      i1_error=0
      if i0 gt 0 or i1 lt nVrb-1 then begin
         for k=k0,k1 do begin
            for j=j0,j1 do begin
               if (i0_error eq 0) and (i0 gt 0) and (fb(i0,j,k) gt warn*big) then begin
                  message,'Non-zero value of fb detected at min(Vra) boundary',/info
                  i0_error=1
               endif
               if (i1_error eq 0) and (i1 lt nVrb-1) and (fb(i1,j,k) gt warn*big) then begin
                  message,'Non-zero value of fb detected at max(Vra) boundary',/info
                  i1_error=1
               endif
            endfor
         endfor
      endif
;
; j0 & j1
;
      j0_error=0
      j1_error=0
      if j0 gt 0 or j1 lt nVxb-1 then begin
         for k=k0,k1 do begin
            for i=i0,i1 do begin
               if (j0_error eq 0) and (j0 gt 0) and (fb(i,j0,k) gt warn*big) then begin
                  message,'Non-zero value of fb detected at min(Vxa) boundary',/info
                  j0_error=1
               endif
               if (j1_error eq 0) and (j1 lt nVxb-1) and (fb(i,j1,k) gt warn*big) then begin
                  message,'Non-zero value of fb detected at max(Vxa) boundary',/info
                  j1_error=1
               endif
            endfor
         endfor
      endif
;
; k0 & k1
;
      k0_error=0
      k1_error=0
      if k0 gt 0 or k1 lt nXb-1 then begin
         for i=i0,i1 do begin
            for j=j0,j1 do begin
               if (k0_error eq 0) and (k0 gt 0) and (fb(i,j,k0) gt warn*big) then begin
                  message,'Non-zero value of fb detected at min(Xa) boundary',/info
                  k0_error=1
               endif
               if (k1_error eq 0) and (k1 lt nxb-1) and (fb(i,j,k1) gt warn*big) then begin
                  message,'Non-zero value of fb detected at max(Xa) boundary',/info
                  k1_error=1
               endif
            endfor
         endfor
      endif
   endif
;
; Rescale
;
   tot_a=dblarr(nXa)
   for k=0,nXa-1 do tot_a(k)=total(Vr2pidVra*(fa(*,*,k)#dVxa))
   tot_b=dblarr(nXb)
   tot_b(k0:k1)=interpol(tot_a,Xa,Xb(k0:k1))
   ii=where(fb gt 0.0,count)
   if count gt 0 then begin
      min_tot=min(fb(ii))
      for k=k0,k1 do begin
         tot=total(Vr2pidVrb*(fb(*,*,k)#dVxb))
         if tot gt min_tot then begin
            if debug then print,prompt+'Density renormalization factor ='+sval(tot_b(k)/tot)
            fb(*,*,k)=fb(*,*,k)*tot_b(k)/tot
         endif
      endfor
   endif

   if debug then begin
;
; na, Uxa, Ta
;
      na=dblarr(nXa)
      Uxa=dblarr(nXa)
      Ta=dblarr(nXa)
      vr2vx2_ran2=dblarr(nvra,nvxa)
      for k=0,nxa-1 do begin
         na(k)=total(Vr2pidVra*(fa(*,*,k)#dVxa))
         if na(k) gt 0 then begin
            Uxa(k)=vtha*total(Vr2pidVra*(fa(*,*,k)#(Vxa*dVxa)))/na(k)
            for i=0,nvra-1 do vr2vx2_ran2(i,*)=vra(i)^2+(vxa-Uxa(k)/vtha)^2
            Ta(k)=(mu*mH)*vtha2*total(Vr2pidVra*((vr2vx2_ran2*fa(*,*,k))#dVxa))/(3*q*na(k))
         endif
      endfor
;
;
; nb, Uxb, Tb
;
      nb=dblarr(nXb)
      Uxb=dblarr(nXb)
      Tb=dblarr(nXb)
      vr2vx2_ran2=dblarr(nvrb,nvxb)
      for k=0,nxb-1 do begin
         nb(k)=total(Vr2pidVrb*(fb(*,*,k)#dVxb))
         if nb(k) gt 0 then begin
            Uxb(k)=vthb*total(Vr2pidVrb*(fb(*,*,k)#(Vxb*dVxb)))/nb(k)
            for i=0,nvrb-1 do vr2vx2_ran2(i,*)=vrb(i)^2+(vxb-Uxb(k)/vthb)^2
            Tb(k)=(mu*mH)*vthb2*total(Vr2pidVrb*((vr2vx2_ran2*fb(*,*,k))#dVxb))/(3*q*nb(k))
         endif
      endfor
;
      data=[na,nb]
      yrange=[min(data),max(data)]
      plot,xa,na,yrange=yrange,title='Density conserved',/nodata
      oplot,xa,na,color=2
      oplot,xb,nb,color=4
      press_return

      data=[Uxa,Uxb]
      yrange=[min(data),max(data)]
      plot,xa,Uxa,yrange=yrange,title='Ux conserved',/nodata
      oplot,xa,Uxa,color=2
      oplot,xb,Uxb,color=4
      press_return

      data=[Ta,Tb]
      yrange=[min(data),max(data)]
      plot,xa,Ta,yrange=yrange,title='T conserved',/nodata
      oplot,xa,Ta,color=2
      oplot,xb,Tb,color=4
      press_return
   endif

   if w_new then begin
      if w1_active then begin
         if debug then print,prompt+'Storing Weight in Weight2'
         vra2=vra
         vxa2=vxa
         Tnorma2=Tnorma
         vrb2=vrb
         vxb2=vxb
         Tnormb2=Tnormb
         weight2=weight
      endif else begin
         if debug then print,prompt+'Storing Weight in Weight1'
         vra1=vra
         vxa1=vxa
         Tnorma1=Tnorma
         vrb1=vrb
         vxb1=vxb
         Tnormb1=Tnormb
         weight1=weight
      endelse
   endif

   SAVE, fb, FILENAME='fb.sav'


   return
   end
