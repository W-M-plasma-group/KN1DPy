import numpy as np
from warnings import warn
from scipy import interpolate

from Make_dVr_dVx import Make_dVr_dVx
from locate import locate
from sval import sval

from global_vars import mH, q

def interp_fvrvxx(fa,Vra,Vxa,Xa,Tnorma,Vrb,Vxb,Xb,Tnormb,do_warn=None, debug=0, correct=1,g=None): # added global variable argument

    #  Input:
    #     Input Distribution function 'a'
    #	fa	- dblarr(nVra,nVxa,nXa) distribution function
    #	Vra	- fltarr(nVra) - radial velocity
    #	Vxa	- fltarr(nVxa) - axial velocity
    #       Xa	- fltarr(nXa)  - spatial coordinate
    #       Tnorma	- float,  Normalization temperature for Vra & Vxa
    #
    #    Desired phase space coordinates of Output Distribution function 'b'
    #	Vrb	- fltarr(nVrb) - radial velocity
    #	Vxb	- fltarr(nVxb) - axial velocity
    #       Xb	- fltarr(nXb)  - spatial coordinate
    #       Tnormb	- float,  Normalization temperature for Vrb & Vxb
    #
    #  Output:
    #     Interpolated Distribution function 'b'
    #	fb	- dblarr(nVrb,nVxb,nXb) distribution function
    #	          fb is scaled if necessary to make its
    #	          digital integral over all velocity space
    #		  equal to that of fa.
    #
    #  Keywords:
    #     Input:
    #	do_warn	- float, acceptable truncation level.
    #		  For interpolations outside the phase space set by
    #		  (Vra, Vxa, Xa), the values of fb are set to zero.
    #		  This may not be acceptable. A test is performed on
    #		  fb at the boundaries. If fb at the boundaries is greater
    #		  than do_warn times the maximum value of fb,
    #		  a warning message is generated.

    prompt='INTERP_FVRVXX => '

    #   Calls INTERP_FVRVXX_internal1 and INTERP_FVRVXX_internal2 common blocks
    vra1=g.INTERP_FVRVXX_internal1_vra1
    vxa1=g.INTERP_FVRVXX_internal1_vxa1
    Tnorma1=g.INTERP_FVRVXX_internal1_Tnorma1
    vrb1=g.INTERP_FVRVXX_internal1_vrb1
    vxb1=g.INTERP_FVRVXX_internal1_vxb1
    Tnormb1=g.INTERP_FVRVXX_internal1_Tnormb1
    weight1=g.INTERP_FVRVXX_internal1_weight1 # fixed capitalization - GG
    
    vra2=g.INTERP_FVRVXX_internal2_vra2
    vxa2=g.INTERP_FVRVXX_internal2_vxa2
    Tnorma2=g.INTERP_FVRVXX_internal2_Tnorma2
    vrb2=g.INTERP_FVRVXX_internal2_vrb2
    vxb2=g.INTERP_FVRVXX_internal2_vxb2
    Tnormb2=g.INTERP_FVRVXX_internal2_Tnormb2
    weight2=g.INTERP_FVRVXX_internal2_weight2

    nvra=Vra.size
    nvxa=Vxa.size
    nxa=Xa.size

    mu=1

    fV=np.sqrt(Tnormb/Tnorma)

    #   Compute Vtha, Vtha2, Vthb and Vthb2

    Vtha=np.sqrt(2*q*Tnorma/(mu*mH))
    Vtha2=Vtha*Vtha # fixed capitalization
    Vthb=np.sqrt(2*q*Tnormb/(mu*mH))
    Vthb2=Vthb*Vthb # fixed capitalization

    if fa[0,0,:].size!=nvra:
        raise Exception('Number of elements in fa(*,0,0) and Vra do not agree!')
    if fa[0,:,0].size!=nvxa:
        raise Exception('Number of elements in fa(0,*,0) and Vxa do not agree!')
    if fa[:,0,0].size!=nxa:
        raise Exception('Number of elements in fa(0,0,*) and Xa do not agree!')

    oki=np.bitwise_and(fV*Vrb<=max(Vra),fV*Vrb>=min(Vra)).nonzero()[0]
    if oki.size<1:
        raise Exception('No values of Vrb are within range of Vra')
    i0,i1=oki[0],oki[-1]

    okj=np.bitwise_and(fV*Vxb<=max(Vxa),fV*Vxb>=min(Vxa)).nonzero()[0]
    if okj.size<1:
        raise Exception('No values of Vxb are within range of Vxa')
    j0,j1=okj[0],okj[-1]

    okk=np.bitwise_and(Xb<=max(Xa), Xb>=min(Xa)).nonzero()[0] # deleted fV terms because they shouldnt be there - GG
    if okk.size<1:
        raise Exception('No values of Xb are within range of Xa')
    k0,k1=okk[0],okk[-1]

    nvrb=Vrb.size
    nvxb=Vxb.size
    nxb=Xb.size

    fb=np.zeros((nxb,nvxb,nvrb)) # added parenthesis - GG

    make_dvr_dvx_out=Make_dVr_dVx(Vra,Vxa)
    Vr2pidVra,VrVr4pidVra,dVxa,vraL,vraR,vxaL,vxaR=make_dvr_dvx_out[:7]
    Vra2Vxa2=make_dvr_dvx_out[11]
    Vr2pidVrb,VrVr4pidVrb,dVxb,vrbL,vrbR,vxbL,vxbR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,Vrb2Vxb2,jpa,jpb,jna,jnb=Make_dVr_dVx(Vrb,Vxb)

    #   Determine if Weight was already computed by checking vra_s,vxa_s,Tnorma_s,vrb_s,vxb_s,Tnormb_s for cases 1 and 2

    w1_active=0
    w1_match=0

    if not vra1 is None:
        w1_active=1
        test=0

        test+=vra1[vra1!=Vra].size
        test+=vxa1[vxa1!=Vxa].size

        test+=Tnorma1[Tnorma1!=Tnorma].size
        test+=Tnormb1[Tnormb1!=Tnormb].size

        test+=vxb1[vxb1!=Vxb].size
        test+=vrb1[vrb1!=Vrb].size

        if test==0:
            w1_match=1

    w2_active=0
    w2_match=0

    if not vra2 is None:
        w2_active=1
        test=0

        test+=vra2[vra2!=Vra].size
        test+=vxa2[vxa2!=Vxa].size

        test+=Tnorma2[Tnorma1!=Tnorma].size
        test+=Tnormb2[Tnormb1!=Tnormb].size

        test+=vxb2[vxb2!=Vxb].size
        test+=vrb2[vrb2!=Vrb].size

        if test==0:
            w2_match=1

    w_new=0

    if w1_match or w2_match:
        if w1_match:
            weight=weight1
            if debug:
                print(prompt+'using weight1')
        if w2_match:
            weight=weight2
            if debug:
                print(prompt+'using weight2')

    #   If not then compute Weight for this combination of Vr and Vx

    else:

        #   Determine Left and Right limits on 'cells' for Vra, Vxa, Vrb, Vxb

        if debug:
            print(prompt+'computing new weight')
        w_new=1

        #   Set area contributions to Weight array

        _weight=np.zeros((nvxa,nvra,nvxb,nvrb))
        weight=np.zeros((nvra*nvxa,nvrb*nvxb))

        for ib in range(nvrb):
            for jb in range(nvxb):
                for ia in range(nvra):
                    vraMin=max([fV*vrbL[ib],vraL[ia]])
                    vraMax=min([fV*vrbR[ib],vraR[ia]])
                    for ja in range(nvxa):
                        vxaMin=max([fV*vxbL[jb],vxaL[ja]])
                        vxaMax=min([fV*vxbR[jb],vxaR[ja]])
                        if vraMax>vraMin and vxaMax>vxaMin:
                            _weight[ja,ia,jb,ib]=2*np.pi*(vraMax**2-vraMin**2)*(vxaMax-vxaMin)/(Vr2pidVrb[ib]*dVxb[jb])

        weight=np.reshape(_weight,weight.shape) # previous version caused error 

    fb_xa=np.zeros((nxa,nvrb*nvxb)) 

    #   Determine fb_xa from weight array

    _fa=np.zeros((nxa,nvra*nvxa))
    _fa=np.reshape(fa,_fa.shape) # previous version caused error
    fb_xa=np.matmul(_fa,weight)

    #   Compute _Wxa and _Ea - these are the desired moments of fb, but on the xa grid

    na=np.zeros(nxa) # fixed capitalization
    _Wxa=np.zeros(nxa)
    _Ea=np.zeros(nxa)

    for k in range(nxa):
        na[k]=np.sum(Vr2pidVra*np.matmul(dVxa,fa[k,:,:])) # fixed typo - GG
        if na[k]>0:
            _Wxa[k]=np.sqrt(Tnorma)*np.sum(Vr2pidVra*np.matmul((Vxa*dVxa),fa[k,:,:]))/na[k]
            _Ea[k]=Tnorma*np.sum(Vr2pidVra*np.matmul(dVxa,(Vra2Vxa2*fa[k,:,:])))/na[k]

    wxa=np.zeros(nxb) # fixed capitalization
    Ea=np.zeros(nxb)

    for k in range(k0,k1+1):
        kL=np.maximum(locate(Xa,Xb[k]),0) # fixed capitalization
        kR=np.minimum(kL+1,Xa.size-1)
        kL=np.minimum(kL,kR-1)

        f=(Xb[k]-Xa[kL])/(Xa[kR]-Xa[kL])
        fb[k,:,:]=np.reshape(fb_xa[kL,:]+(fb_xa[kR,:]-fb_xa[kL,:])*f,fb[k,:,:].shape) # previous version caused error
        wxa[k]=_Wxa[kL]+(_Wxa[kR]-_Wxa[kL])*f
        Ea[k]=_Ea[kL]+(_Ea[kR]-_Ea[kL])*f

    #   Correct fb so that it has the same Wx and E moments as fa

    if correct:

        #   Process each spatial location

        AN=np.zeros((2,nvxb,nvrb))
        BN=np.zeros((2,nvxb,nvrb))

        sgn=np.array([-1,1])

        for k in range(nxb):
            allow_neg=0

            #   Compute nb, Wxb, and Eb - these are the current moments of fb

            nb=np.sum(Vr2pidVrb*np.matmul(dVxb,fb[k,:,:]))
            if nb>0:
                
                #   Entry point for iteration - 'correct' tag in original code
                #   Since Python doesn't have goto, a while loop was used

                goto_correct=True
                while goto_correct:
                    goto_correct=False
                    nb=np.sum(Vr2pidVrb*np.matmul(dVxb,fb[k,:,:]))
                    Wxb=np.sqrt(Tnormb)*np.sum(Vr2pidVrb*np.matmul(Vxb*dVxb,fb[k,:,:]))/nb
                    Eb=Tnormb*np.sum(Vr2pidVrb*np.matmul(dVxb,Vrb2Vxb2*fb[k,:,:]))/nb

                    #   Compute Nij from fb, padded with zeros

                    Nij=np.zeros((nvxb+2,nvrb+2))
                    Nij[1:-1,1:-1]=fb[k,:,:]*Vol/nb

                    #   Set Cutoff and remove Nij very close to zero

                    cutoff=1.0e-6*np.max(Nij)
                    ii=np.argwhere(np.bitwise_and(abs(Nij)<cutoff,abs(Nij)>0))
                    for i in ii:
                        Nij[tuple(i)]=0
                    if max(Nij[2,:])<=0:
                        allow_neg=1
                    Nijp1_vx_Dvx=np.roll(Nij*Vx_DVx,-1,0)
                    Nij_vx_Dvx=Nij*Vx_DVx
                    Nijm1_vx_Dvx=np.roll(Nij*Vx_DVx,1,0)
                    Nip1j_vr_Dvr=np.roll(Nij*Vr_DVr,-1,1)
                    Nij_vr_Dvr=Nij*Vr_DVr
                    Nim1j_vr_Dvr=np.roll(Nij*Vr_DVr,1,1)

                    #   Compute Ap, Am, Bp, and Bm (0=p 1=m)

                    _AN=np.roll(Nij*Vth_DVx,1,0)-Nij*Vth_DVx
                    AN[0,:,:]=_AN[1:nvxb+1,1:nvrb+1]
                    _AN=-np.roll(Nij*Vth_DVx,-1,0)
                    AN[1,:,:]=_AN[1:nvxb+1,1:nvrb+1]

                    BN[0,jpa+1:jpb+1,:]=Nijm1_vx_Dvx[jpa+2:jpb+2,1:nvrb+1]-Nij_vx_Dvx[jpa+2:jpb+2,1:nvrb+1]
                    BN[0,jpa,:]=-Nij_vx_Dvx[jpa+1,1:nvrb+1]
                    BN[0,jnb,:]=Nij_vx_Dvx[jnb+1,1:nvrb+1]
                    BN[0,jna:jnb,:]=-Nijp1_vx_Dvx[jna+1:jnb+1,1:nvrb+1]+Nij_vx_Dvx[jna+1:jnb+1,1:nvrb+1]
                    BN[0,:,:]=BN[0,:,:]+Nim1j_vr_Dvr[1:nvxb+1,1:nvrb+1]-Nij_vr_Dvr[1:nvxb+1,1:nvrb+1]

                    BN[1,jpa+1:jpb+1,:]=-Nijp1_vx_Dvx[jpa+2:jpb+2,1:nvrb+1]+Nij_vx_Dvx[jpa+2:jpb+2,1:nvrb+1]
                    BN[1,jpa,:]=-Nijp1_vx_Dvx[jpa+1,1:nvrb+1]
                    BN[1,jnb,:]=Nijm1_vx_Dvx[jnb+1,1:nvrb+1]
                    BN[1,jna:jnb,:]=Nijm1_vx_Dvx[jna+1:jnb+1,1:nvrb+1]-Nij_vx_Dvx[jna+1:jnb+1,1:nvrb+1]
                    BN[1,:,1:nvrb]=BN[1,:,1:nvrb]-Nip1j_vr_Dvr[1:nvxb+1,2:nvrb+1]+Nij_vr_Dvr[1:nvxb+1,2:nvrb+1]
                    BN[1,:,0]=BN[1,:,0]-Nip1j_vr_Dvr[1:nvxb+1,1]

                    #   If negative values for Nij must be allowed, then add postive particles to i=0 and negative particles to i=1 (beta is negative here)

                    if allow_neg:
                        BN[1,:,0]=BN[1,:,0]-Nij_vr_Dvr[1:nvxb+1,1]
                        BN[1,:,1]=BN[1,:,1]+Nij_vr_Dvr[1:nvxb+1,1]

                    #   Remove padded zeros in Nij

                    Nij=Nij[1:nvxb+1,1:nvrb+1]

                    #   Cycle through 4 possibilies of sign(alpha),sign(beta)

                    TB1=np.zeros(2); TB2=np.zeros(2)

                    for ia in range(2):

                        #   Compute TA1, TA2

                        TA1=np.sqrt(Tnormb)*np.sum(np.matmul(Vxb,AN[ia,:,:]))
                        TA2=Tnormb*np.sum(Vrb2Vxb2*AN[ia,:,:])
                        for ib in range(2):

                            #   Compute TB1, TB2

                            if TB1[ib]==0:
                                TB1[ib]=np.sqrt(Tnormb)*np.sum(np.matmul(Vxb,BN[ib,:,:]))
                            if TB2[ib]==0:
                                TB2[ib]=Tnormb*np.sum(Vrb2Vxb2*BN[ib,:,:])

                            denom=TA2*TB1[ib]-TA1*TB2[ib]
                            beta=0
                            alpha=0

                            if denom!=0 and TA1!=0:
                                beta=(TA2*(wxa[k]-Wxb)-TA1*(Ea[k]-Eb))/denom # fixed capitalization
                                alpha=(wxa[k]-Wxb-TB1[ib]*beta)/TA1

                            do_break=alpha*sgn[ia]>0 and beta*sgn[ib]>0
                            if do_break:
                                break
                        if do_break:
                            break

                    #   Entry point for 'alpha_beta' tag from original code

                    RHS=AN[ia,:,:]*alpha+BN[ib,:,:]*beta

                    #   Are there locations where Nij = 0.0 and RHS is negative?
                    #   ii=where(Nij eq 0.0 and RHS lt 0.0,count) was in the original code, I don't think it does anything

                    s=1
                    if not allow_neg:
                        ii=np.nonzero(Nij)
                        if np.size(ii)>0:
                            s=min(1/np.max(-RHS[ii]/Nij[ii]),1)
                    fb[k,:,:]=nb*(Nij+s*RHS)/Vol # fixed capitalization

                    goto_correct=s<1

    if do_warn!=None:

        #   Test Boundaries:

        #   i0 & i1 

        big=np.max(fb)

        i0_error=0
        i1_error=0
        if i0>0 or i1<nvrb-1:
            for k in range(k0,k1+1):
                for j in range(j0,j1+1):
                    if (i0_error == 0) and (i0 > 0) and (fb[k,j,i0] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Vra) boundary')
                        i0_error=1
                    if (i1_error == 0) and (i1 < nvrb-1) and (fb[k,j,i1] > do_warn*big): # fixed capitalization
                        warn('Non-zero value of fb detected at max(Vra) boundary')
                        i1_error=1

        #   j0 & j1

        j0_error=0
        j1_error=0
        if j0>0 or j1<nvxb-1:
            for k in range(k0,k1+1):
                for i in range(i0,i1+1):
                    if (j0_error == 0) and (j0 > 0) and (fb[k,j0,i] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Vxa) boundary')
                        j0_error=1
                    if (j1_error == 0) and (j1 < nvxb-1) and (fb[k,j1,i] > do_warn*big): # fixed capitalization
                        warn('Non-zero value of fb detected at max(Vxa) boundary')
                        j1_error=1

        #   k0 & k1

        k0_error=0
        k1_error=1
        if k0 > 0 or k1 < nxb-1:
            for i in range(i0,i1+1):
                for j in range(j0,j1+1):
                    if (k0_error == 0) and (k0 > 0) and (fb[k0,j,i] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Xa) boundary')
                        k0_error=1
                    if (k1_error == 0) and (k1 < nxb-1) and (fb[k1,j,i] > do_warn*big):
                        warn('Non-zero value of fb detected at max(Xa) boundary')
                        k1_error=1

    #   Rescale

    tot_a=np.zeros(nxa)
    for k in range(nxa):
        tot_a[k]=np.sum(Vr2pidVra*np.matmul(dVxa,fa[k,:,:]))
    tot_b=np.zeros(nxb)
    tot_b[k0:k1+1]=interpolate.interp1d(Xa,tot_a,fill_value="extrapolate")(Xb[k0:k1+1])
    ii=np.argwhere(fb>0)
    if ii.size>0: # replaced fb with ii
        min_tot=np.min(np.array([fb[tuple(i)] for i in ii]))
        for k in range(k0,k1+1):
            tot=np.sum(Vr2pidVrb*np.matmul(dVxb,fb[k,:,:]))
            if tot>min_tot:
                if debug:
                    print(prompt+'Density renormalization factor ='+sval(tot_b[k]/tot))

    if debug:

        #   na, Uxa, Ta

        na=np.zeros(nxa)
        Uxa=np.zeros(nxa)
        Ta=np.zeros(nxa)
        vr2vx2_ran2=np.zeros((nvxa,nvra)) # fixed np.zeros() call

        for k in range(nxa):
            na[k]=np.sum(Vr2pidVra*np.matmul(dVxa,fa[k,:,:])) # fixed capitalization
            if na[k]>0:
                Uxa[k]=Vtha*np.sum(Vr2pidVra*np.matmul(Vxa*dVxa,fa[k,:,:]))/na[k]
                for i in range(nvra):
                    vr2vx2_ran2[:,i]=Vra[i]**2+(Vxa-Uxa[k]/Vtha)**2 # fixed capitalization
                Ta[k]=mu*mH*Vtha2*np.sum(Vr2pidVra*np.matmul(dVxa,vr2vx2_ran2*fa[k,:,:]))/(3*q*na[k])

        #   nb, Uxb, Tb

        nb=np.zeros(nxb)
        Uxb=np.zeros(nxb)
        Tb=np.zeros(nxb)
        vr2vx2_ran2=np.zeros((nvxb,nvrb)) # fixed np.zeros() call

        for k in range(nxb):
            nb[k]=np.sum(Vr2pidVrb*np.matmul(dVxb,fb[k,:,:]))
            if nb[k]>0:
                Uxb[k]=Vthb*np.sum(Vr2pidVrb*np.matmul(Vxb*dVxb,fb[k,:,:]))/nb[k] # fixed typo
                for i in range(nvrb):
                    vr2vx2_ran2[:,i]=Vrb[i]**2+(Vxb-Uxb[k]/Vthb)**2 # fixed capitalization
                Tb[k]=mu*mH*Vthb*np.sum(Vr2pidVrb*np.matmul(dVxb,vr2vx2_ran2*fb[k,:,:]))/(3*q*nb[k])

        #   Plotting stuff was here in the original code
        #   May be added later, but has been left out for now
    #   update INTERP_FVRVXX_internal1 and INTERP_FVRVXX_internal2 common blocks
    if w_new:
        if w1_active:
            if debug:
                print(prompt+'Storing Weight in Weight2')
            g.INTERP_FVRVXX_internal2_vra2=Vra
            g.INTERP_FVRVXX_internal2_vxa2=Vxa
            g.INTERP_FVRVXX_internal2_Tnorma2=Tnorma
            g.INTERP_FVRVXX_internal2_vrb2=Vrb
            g.INTERP_FVRVXX_internal2_vxb2=Vxb
            g.INTERP_FVRVXX_internal2_Tnormb2=Tnormb
            g.INTERP_FVRVXX_internal2_weight2=weight
        else:
            if debug:
                print(prompt+'Storing Weight in Weight1')
            g.INTERP_FVRVXX_internal1_vra1=Vra
            g.INTERP_FVRVXX_internal1_vxa1=Vxa
            g.INTERP_FVRVXX_internal1_Tnorma1=Tnorma
            g.INTERP_FVRVXX_internal1_vrb1=Vrb
            g.INTERP_FVRVXX_internal1_vxb1=Vxb
            g.INTERP_FVRVXX_internal1_Tnormb1=Tnormb
            g.INTERP_FVRVXX_internal1_weight1=weight

    return fb
