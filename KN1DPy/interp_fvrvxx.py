import numpy as np
import copy
from warnings import warn
from scipy import interpolate

from .make_dvr_dvx import make_dvr_dvx
from .locate import locate
from .sval import sval
from .kinetic_mesh import kinetic_mesh

from .common import constants as CONST
from .common.INTERP_FVRVXX import INTERP_FVRVXX_internal

def interp_fvrvxx(fa, mesh_a : kinetic_mesh, mesh_b : kinetic_mesh, internal : INTERP_FVRVXX_internal, 
                  do_warn=None, debug=0, correct=1, debug_flag = 0): #NOTE Debug flag added to mark specific function calls, remove later

    # NOTE Passing full mesh may be unneccessary, works for this code, but would need to be refactored for other projects
    # NOTE If removing mesh passing, will need to change code, mesh variables have replaced several variables here
    
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
    vra1 = internal.segment1.vra
    vxa1 = internal.segment1.vxa
    Tnorma1 = internal.segment1.Tnorma
    vrb1 = internal.segment1.vrb
    vxb1 = internal.segment1.vxb
    Tnormb1 = internal.segment1.Tnormb
    weight1 = internal.segment1.weight # fixed capitalization - GG
    
    vra2 = internal.segment2.vra
    vxa2 = internal.segment2.vxa
    Tnorma2 = internal.segment2.Tnorma
    vrb2 = internal.segment2.vrb
    vxb2 = internal.segment2.vxb
    Tnormb2 = internal.segment2.Tnormb
    weight2 = internal.segment2.weight

    nvra = mesh_a.vr.size
    nvxa = mesh_a.vx.size
    nxa = mesh_a.x.size

    mu=1

    fV=np.sqrt(mesh_b.Tnorm/mesh_a.Tnorm)

    #   Compute Vtha, Vtha2, Vthb and Vthb2

    Vtha = np.sqrt((2*CONST.Q*mesh_a.Tnorm) / (mu*CONST.H_MASS))
    Vtha2 = Vtha*Vtha # fixed capitalization
    Vthb = np.sqrt((2*CONST.Q*mesh_b.Tnorm) / (mu*CONST.H_MASS))
    Vthb2 = Vthb*Vthb # fixed capitalization

    if fa[:,0,0].size != nvra:
        raise Exception('Number of elements in fa(*,0,0) and Vra do not agree!')
    if fa[0,:,0].size != nvxa:
        raise Exception('Number of elements in fa(0,*,0) and Vxa do not agree!')
    if fa[0,0,:].size != nxa:
        raise Exception('Number of elements in fa(0,0,*) and Xa do not agree!')

    oki = np.where((fV*mesh_b.vr <= max(mesh_a.vr)) & (fV*mesh_b.vr >= min(mesh_a.vr)))[0]
    if len(oki) < 1:
        raise Exception('No values of Vrb are within range of Vra')
    i0, i1 = oki[0], oki[-1]

    okj = np.where((fV*mesh_b.vx <= max(mesh_a.vx)) & (fV*mesh_b.vx >= min(mesh_a.vx)))[0]
    if okj.size < 1:
        raise Exception('No values of Vxb are within range of Vxa')
    j0, j1 = okj[0], okj[-1]

    okk = np.where((mesh_b.x <= max(mesh_a.x)) & (mesh_b.x >= min(mesh_a.x)))[0]
    if okk.size < 1:
        raise Exception('No values of Xb are within range of Xa')
    k0, k1 = okk[0], okk[-1]

    # print(i0, i1, j0, j1, k0, k1)
    # input()

    nvrb = mesh_b.vr.size
    nvxb = mesh_b.vx.size
    nxb = mesh_b.x.size

    fb = np.zeros((nvrb, nvxb, nxb))

    make_dvr_dvx_out = make_dvr_dvx(mesh_a.vr,mesh_a.vx)
    Vr2pidVra,VrVr4pidVra,dVxa,vraL,vraR,vxaL,vxaR = make_dvr_dvx_out[:7]
    Vra2Vxa2 = make_dvr_dvx_out[11]

    make_dvr_dvx_out = make_dvr_dvx(mesh_b.vr,mesh_b.vx)
    Vr2pidVrb,VrVr4pidVrb,dVxb,vrbL,vrbR,vxbL,vxbR,Vol,Vth_DVx,Vx_DVx,Vr_DVr,Vrb2Vxb2,jpa,jpb,jna,jnb = make_dvr_dvx_out

    #   Determine if Weight was already computed by checking vra_s,vxa_s,Tnorma_s,vrb_s,vxb_s,Tnormb_s for cases 1 and 2

    w1_active = 0
    w1_match = 0

    if (not vra1 is None) and (vra1.size == mesh_a.vr.size): #NOTE Might need a better fix for this, but works for now
        w1_active = 1
        test = 0

        test += vra1[vra1 != mesh_a.vr].size
        test += vxa1[vxa1 != mesh_a.vx].size

        test += Tnorma1[Tnorma1 != mesh_a.Tnorm].size
        test += Tnormb1[Tnormb1 != mesh_b.Tnorm].size

        test += vxb1[vxb1 != mesh_b.vx].size
        test += vrb1[vrb1 != mesh_b.vr].size

        if test == 0:
            w1_match = 1

    w2_active = 0
    w2_match = 0

    if (not vra2 is None) and (vra2.size == mesh_a.vr.size):
        w2_active = 1
        test = 0

        test += vra2[vra2 != mesh_a.vr].size
        test += vxa2[vxa2 != mesh_a.vx].size

        test += Tnorma2[Tnorma1 != mesh_a.Tnorm].size
        test += Tnormb2[Tnormb1 != mesh_b.Tnorm].size

        test += vxb2[vxb2 != mesh_b.vx].size
        test += vrb2[vrb2 != mesh_b.vr].size

        if test == 0:
            w2_match = 1

    w_new = 0

    if w1_match or w2_match:
        if w1_match:
            weight = weight1
            if debug:
                print(prompt+'using weight1')
        if w2_match:
            weight = weight2
            if debug:
                print(prompt+'using weight2')

    #   If not then compute Weight for this combination of Vr and Vx

    else:

        #   Determine Left and Right limits on 'cells' for Vra, Vxa, Vrb, Vxb

        if debug:
            print(prompt+'computing new weight')
        w_new = 1

        #   Set area contributions to Weight array

        _weight = np.zeros((nvrb,nvxb,nvra,nvxa))
        weight = np.zeros((nvrb*nvxb,nvra*nvxa))

        for ib in range(nvrb):
            for jb in range(nvxb):
                for ia in range(nvra):
                    vraMin = max([fV*vrbL[ib], vraL[ia]])
                    vraMax = min([fV*vrbR[ib], vraR[ia]])
                    for ja in range(nvxa):
                        vxaMin = max([fV*vxbL[jb], vxaL[ja]])
                        vxaMax = min([fV*vxbR[jb], vxaR[ja]])
                        if vraMax > vraMin and vxaMax > vxaMin:
                            _weight[ib,jb,ia,ja] = 2*np.pi*(vraMax**2 - vraMin**2)*(vxaMax - vxaMin) / (Vr2pidVrb[ib]*dVxb[jb])

        weight = np.reshape(_weight, weight.shape, order = 'F') # previous version caused error

    # print("weight", weight.T[np.nonzero(weight.T)].size)
    # input()

    fb_xa = np.zeros((nvrb*nvxb,nxa)) 

    #   Determine fb_xa from weight array

    _fa = np.reshape(fa, (nvra*nvxa, nxa), order = 'F')
    fb_xa = weight @ _fa

    # print("fb_xa", fb_xa)
    # print(fb_xa.shape)
    # input()

    #   Compute _Wxa and _Ea - these are the desired moments of fb, but on the xa grid

    na = np.zeros(nxa)
    _Wxa = np.zeros(nxa)
    _Ea = np.zeros(nxa)

    for k in range(nxa):
        na[k] = np.sum(Vr2pidVra*(fa[:,:,k] @ dVxa))
        if na[k] > 0:
            _Wxa[k] = np.sqrt(mesh_a.Tnorm)*np.sum(Vr2pidVra*(fa[:,:,k] @ (mesh_a.vx*dVxa)))/na[k]
            if debug_flag:
                print("_wxa", _Wxa[k])
                print("a", fa[:,:,k].T)
                print("b", (mesh_a.vx*dVxa))
                print("c", (fa[:,:,k] @ (mesh_a.vx*dVxa)))
                input()
            _Ea[k] = mesh_a.Tnorm*np.sum(Vr2pidVra*((Vra2Vxa2*fa[:,:,k]) @ dVxa))/na[k]
    # print("na", na)
    # input()

    wxa = np.zeros(nxb)
    Ea = np.zeros(nxb)

    for k in range(k0, k1+1):
        kL = np.maximum(locate(mesh_a.x, mesh_b.x[k]), 0)
        kR = np.minimum(kL+1, mesh_a.x.size-1)
        kL = np.minimum(kL, kR-1)

        f = (mesh_b.x[k] - mesh_a.x[kL]) / (mesh_a.x[kR] - mesh_a.x[kL])
        fb[:,:,k] = np.reshape((fb_xa[:,kL] + (fb_xa[:,kR] - fb_xa[:,kL])*f), fb[:,:,k].shape, order='F')
        wxa[k] = _Wxa[kL] + (_Wxa[kR] - _Wxa[kL])*f
        if debug_flag:
            print("klkr", kL, kR)
            print("wxa1", wxa[k])
            print("wxa2", _Wxa[kL])
            print("wxa3", _Wxa[kR])
            input()
        Ea[k]=_Ea[kL] + (_Ea[kR] - _Ea[kL])*f
    # print("fb", fb)
    print("wxa", wxa)
    # print("Ea", Ea)
    # input()

    #   Correct fb so that it has the same Wx and E moments as fa

    if correct:

        #   Process each spatial location

        AN = np.zeros((nvrb, nvxb, 2))
        BN = np.zeros((nvrb, nvxb, 2))

        sgn = np.array([1, -1])

        for k in range(nxb):
            allow_neg = 0

            #   Compute nb, Wxb, and Eb - these are the current moments of fb

            nb = np.sum(Vr2pidVrb*(fb[:,:,k] @ dVxb))
            # print("nb", nb)
            # input()
            if nb > 0:
                
                #   Entry point for iteration - 'correct' tag in original code
                #   Since Python doesn't have goto, a while loop was used

                goto_correct = True
                while goto_correct:
                    goto_correct = False
                    nb = np.sum(Vr2pidVrb*(fb[:,:,k] @ dVxb))
                    Wxb = np.sqrt(mesh_b.Tnorm)*np.sum(Vr2pidVrb*(fb[:,:,k] @ (mesh_b.vx*dVxb))) / nb
                    Eb = mesh_b.Tnorm*np.sum(Vr2pidVrb*((Vrb2Vxb2*fb[:,:,k]) @ dVxb)) / nb
                    # print("nb", nb)
                    # print("Wxb", Wxb)
                    # print("Eb", Eb)
                    # input()

                    #   Compute Nij from fb, padded with zeros

                    Nij = np.zeros((nvrb+2,nvxb+2))
                    Nij[1:nvrb+1,1:nvxb+1] = fb[:,:,k]*Vol / nb

                    #   Set Cutoff and remove Nij very close to zero

                    cutoff = 1.0e-6*np.max(Nij)
                    ii = np.where((abs(Nij) < cutoff) & (abs(Nij) > 0))
                    if ii[0].size > 0:
                        Nij[ii] = 0.0

                    if max(Nij[2,:]) <= 0:
                        allow_neg = 1

                    Nijp1_vx_Dvx = np.roll(Nij*Vx_DVx, shift=-1, axis=1)
                    Nij_vx_Dvx = Nij*Vx_DVx
                    Nijm1_vx_Dvx = np.roll(Nij*Vx_DVx, shift=1, axis=1)
                    Nip1j_vr_Dvr = np.roll(Nij*Vr_DVr, shift=-1, axis=0)
                    Nij_vr_Dvr = Nij*Vr_DVr
                    Nim1j_vr_Dvr = np.roll(Nij*Vr_DVr, shift=1, axis=0)

                    #   Compute Ap, Am, Bp, and Bm (0=p 1=m)

                    _AN                 = np.roll(Nij*Vth_DVx, shift=1, axis=1) - Nij*Vth_DVx
                    AN[:,:,0]           = copy.copy(_AN[1:nvrb+1,1:nvxb+1])
                    
                    _AN                 = -np.roll(Nij*Vth_DVx, shift=-1, axis=1) + Nij*Vth_DVx
                    AN[:,:,1]           = copy.copy(_AN[1:nvrb+1,1:nvxb+1])

                    BN[:,jpa+1:jpb+1,0] =  Nijm1_vx_Dvx[1:nvrb+1,jpa+2:jpb+2] - Nij_vx_Dvx[1:nvrb+1,jpa+2:jpb+2]
                    BN[:,jpa,0]         = -Nij_vx_Dvx[1:nvrb+1,jpa+1]
                    BN[:,jnb,0]         =  Nij_vx_Dvx[1:nvrb+1,jnb+1]
                    BN[:,jna:jnb,0]     = -Nijp1_vx_Dvx[1:nvrb+1,jna+1:jnb+1] + Nij_vx_Dvx[1:nvrb+1,jna+1:jnb+1]
                    BN[:,:,0]           =  BN[:,:,0] + Nim1j_vr_Dvr[1:nvrb+1,1:nvxb+1] - Nij_vr_Dvr[1:nvrb+1,1:nvxb+1]

                    BN[:,jpa+1:jpb+1,1] = -Nijp1_vx_Dvx[1:nvrb+1,jpa+2:jpb+2] + Nij_vx_Dvx[1:nvrb+1,jpa+2:jpb+2]
                    BN[:,jpa,1]         = -Nijp1_vx_Dvx[1:nvrb+1,jpa+1]
                    BN[:,jnb,1]         =  Nijm1_vx_Dvx[1:nvrb+1,jnb+1]
                    BN[:,jna:jnb,1]     =  Nijm1_vx_Dvx[1:nvrb+1,jna+1:jnb+1] - Nij_vx_Dvx[1:nvrb+1,jna+1:jnb+1]
                    BN[1:nvrb,:,1]      =  BN[1:nvrb,:,1] - Nip1j_vr_Dvr[2:nvrb+1,1:nvxb+1] + Nij_vr_Dvr[2:nvrb+1,1:nvxb+1]
                    BN[0,:,1]           =  BN[0,:,1] - Nip1j_vr_Dvr[1,1:nvxb+1]

                    # print("AN", AN.T)
                    # print("BN", BN.T)
                    # input()

                    #   If negative values for Nij must be allowed, then add postive particles to i=0 and negative particles to i=1 (beta is negative here)

                    if allow_neg:
                        BN[0,:,1] = BN[0,:,1] - Nij_vr_Dvr[1,1:nvxb+1]
                        BN[1,:,1] = BN[1,:,1] + Nij_vr_Dvr[1,1:nvxb+1]

                    #   Remove padded zeros in Nij

                    Nij = Nij[1:nvrb+1,1:nvxb+1]
                    

                    #   Cycle through 4 possibilies of sign(alpha),sign(beta)

                    TB1 = np.zeros(2, float)
                    TB2 = np.zeros(2, float)

                    for ia in range(2):
                        # print("a")
                        #   Compute TA1, TA2

                        TA1 = np.sqrt(mesh_b.Tnorm)*np.sum((AN[:,:,ia] @ mesh_b.vx))
                        TA2 = mesh_b.Tnorm*np.sum(Vrb2Vxb2*AN[:,:,ia])
                        for ib in range(2):
                            # print("b")

                            #   Compute TB1, TB2

                            if TB1[ib] == 0:
                                TB1[ib] = np.sqrt(mesh_b.Tnorm)*np.sum((BN[:,:,ib] @ mesh_b.vx))
                            if TB2[ib] == 0:
                                TB2[ib] = mesh_b.Tnorm*np.sum(Vrb2Vxb2*BN[:,:,ib])

                            denom = TA2*TB1[ib] - TA1*TB2[ib]
                            if debug_flag:
                                print("denom", denom)
                                # print(TA2*TB1[ib])
                                # print(TA1*TB2[ib])
                                # print(TB2[ib])
                                # print(BN[:,:,ib].T)
                                input()
                            beta = 0
                            alpha = 0

                            if (denom != 0) and (TA1 != 0):
                                beta = (TA2*(wxa[k] - Wxb) - TA1*(Ea[k] - Eb))/denom # fixed capitalization
                                alpha = (wxa[k] - Wxb - TB1[ib]*beta)/TA1

                            do_break = ((alpha*sgn[ia]) > 0) and ((beta*sgn[ib]) > 0)
                            # print("break", do_break)
                            # print(Wxb)
                            if do_break:
                                break
                        if do_break:
                            break

                    #   Entry point for 'alpha_beta' tag from original code

                    RHS = AN[:,:,ia]*alpha + BN[:,:,ib]*beta

                    #   Are there locations where Nij = 0.0 and RHS is negative?
                    #   ii=where(Nij eq 0.0 and RHS lt 0.0,count) was in the original code, I don't think it does anything

                    s = 1
                    if not allow_neg:
                        ii = np.nonzero(Nij)
                        if np.size(ii) > 0:
                            s = min(1/np.max(-RHS[ii]/Nij[ii]),1)
                    
                    fb[:,:,k] = nb*(Nij + s*RHS)/Vol # fixed capitalization

                    # print("goto", goto_correct)
                    goto_correct = (s < 1)
    if(debug_flag != 0):
        ii = np.nonzero(fb.reshape(fb.size, order='F'))
        print("fSHAnz", ii)
        input()

    if do_warn != None:

        #   Test Boundaries:

        #   i0 & i1 

        big = np.max(fb)

        i0_error = 0
        i1_error = 0
        if (i0 > 0) or (i1 < nvrb-1):
            for k in range(k0, k1+1):
                for j in range(j0, j1+1):
                    if (i0_error == 0) and (i0 > 0) and (fb[i0,j,k] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Vra) boundary')
                        i0_error = 1
                    if (i1_error == 0) and (i1 < nvrb-1) and (fb[i1,j,k] > do_warn*big): # fixed capitalization
                        warn('Non-zero value of fb detected at max(Vra) boundary')
                        i1_error = 1

        #   j0 & j1

        j0_error = 0
        j1_error = 0
        if (j0 > 0) or (j1 < nvxb-1):
            for k in range(k0, k1+1):
                for i in range(i0, i1+1):
                    if (j0_error == 0) and (j0 > 0) and (fb[i,j0,k] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Vxa) boundary')
                        j0_error=1
                    if (j1_error == 0) and (j1 < nvxb-1) and (fb[i,j1,k] > do_warn*big): # fixed capitalization
                        warn('Non-zero value of fb detected at max(Vxa) boundary')
                        j1_error=1

        #   k0 & k1

        k0_error = 0
        k1_error = 0
        if (k0 > 0) or (k1 < nxb-1):
            for i in range(i0, i1+1):
                for j in range(j0, j1+1):
                    if (k0_error == 0) and (k0 > 0) and (fb[i,j,k0] > do_warn*big):
                        warn('Non-zero value of fb detected at min(Xa) boundary')
                        k0_error = 1
                    if (k1_error == 0) and (k1 < nxb-1) and (fb[i,j,k1] > do_warn*big):
                        warn('Non-zero value of fb detected at max(Xa) boundary')
                        k1_error = 1

    #   Rescale

    tot_a = np.zeros(nxa)
    for k in range(nxa):
        tot_a[k] = np.sum(Vr2pidVra*(fa[:,:,k] @ dVxa))
    tot_b = np.zeros(nxb)
    tot_b[k0:k1+1] = interpolate.interp1d(mesh_a.x,tot_a,fill_value="extrapolate")(mesh_b.x[k0:k1+1])
    ii = np.where(fb>0)
    if ii[0].size > 0: # replaced fb with ii
        min_tot = np.min(np.array(fb[ii])) #(np.array([fb[tuple(i)] for i in ii]))
        for k in range(k0, k1+1):
            tot = np.sum(Vr2pidVrb*(fb[:,:,k] @ dVxb))
            if tot > min_tot:
                if debug:
                    print(prompt+'Density renormalization factor ='+sval(tot_b[k]/tot))
                fb[:,:,k] = fb[:,:,k]*tot_b[k]/tot

    if debug:

        #   na, Uxa, Ta

        na = np.zeros(nxa)
        Uxa = np.zeros(nxa)
        Ta = np.zeros(nxa)
        vr2vx2_ran2 = np.zeros((nvra,nvxa)) # fixed np.zeros() call

        for k in range(nxa):
            na[k] = np.sum(Vr2pidVra*(fa[:,:,k] @ dVxa)) # fixed capitalization
            if na[k] > 0:
                Uxa[k] = Vtha*np.sum(Vr2pidVra*(fa[k,:,:] @ (mesh_a.vx*dVxa)))/na[k]
                for i in range(nvra):
                    vr2vx2_ran2[i,:] = mesh_a.vr[i]**2 + (mesh_a.vx - Uxa[k]/Vtha)**2 # fixed capitalization
                Ta[k] = mu*CONST.H_MASS*Vtha2*np.sum(Vr2pidVra*((vr2vx2_ran2*fa[:,:,k]) @ dVxa))/(3*CONST.Q*na[k])

        #   nb, Uxb, Tb

        nb = np.zeros(nxb)
        Uxb = np.zeros(nxb)
        Tb = np.zeros(nxb)
        vr2vx2_ran2 = np.zeros((nvrb,nvxb)) # fixed np.zeros() call

        for k in range(nxb):
            nb[k] = np.sum(Vr2pidVrb*(fb[:,:,k] @ dVxb))
            if nb[k] > 0:
                Uxb[k] = Vthb*np.sum(Vr2pidVrb*(fb[:,:,k] @ (mesh_b.vx*dVxb)))/nb[k] # fixed typo
                for i in range(nvrb):
                    vr2vx2_ran2[i,:] = mesh_b.vr[i]**2 + (mesh_b.vx - Uxb[k]/Vthb)**2 # fixed capitalization
                Tb[k] = mu*CONST.H_MASS*Vthb2*np.sum(Vr2pidVrb*((vr2vx2_ran2*fb[:,:,k]) @ dVxb))/(3*CONST.Q*nb[k])

        #   Plotting stuff was here in the original code
        #   May be added later, but has been left out for now
    #   update INTERP_FVRVXX_internal1 and INTERP_FVRVXX_internal2 common blocks

    if w_new:
        if w1_active:
            if debug:
                print(prompt+'Storing Weight in Weight2')
            internal.segment2.vra = mesh_a.vr
            internal.segment2.vxa = mesh_a.vx
            internal.segment2.Tnorma = mesh_a.Tnorm
            internal.segment2.vrb = mesh_b.vr
            internal.segment2.vxb = mesh_b.vx
            internal.segment2.Tnormb = mesh_b.Tnorm
            internal.segment2.weight = weight
        else:
            if debug:
                print(prompt+'Storing Weight in Weight1')
            internal.segment1.vra = mesh_a.vr
            internal.segment1.vxa = mesh_a.vx
            internal.segment1.Tnorma = mesh_a.Tnorm
            internal.segment1.vrb = mesh_b.vr
            internal.segment1.vxb = mesh_b.vx
            internal.segment1.Tnormb = mesh_b.Tnorm
            internal.segment1.weight = weight

    # print("fb", fb)
    # input()

    return fb
