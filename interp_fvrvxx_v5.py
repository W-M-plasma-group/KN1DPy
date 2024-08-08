import numpy as np
from warnings import warn
from scipy import interpolate

from Make_dVr_dVx import Make_dVr_dVx
from locate import locate
from sval import sval

from global_vars import mH, q

def interp_fvrvxx(fa, Vra, Vxa, Xa, Tnorma, Vrb, Vxb, Xb, Tnormb, do_warn=None, debug=0, correct=1, g=None):
    prompt = 'INTERP_FVRVXX => '

    # Calls INTERP_FVRVXX_internal1 and INTERP_FVRVXX_internal2 common blocks
    vra1 = g.INTERP_FVRVXX_internal1_vra1
    vxa1 = g.INTERP_FVRVXX_internal1_vxa1
    Tnorma1 = g.INTERP_FVRVXX_internal1_Tnorma1
    vrb1 = g.INTERP_FVRVXX_internal1_vrb1
    vxb1 = g.INTERP_FVRVXX_internal1_vxb1
    Tnormb1 = g.INTERP_FVRVXX_internal1_Tnormb1
    weight1 = g.INTERP_FVRVXX_internal1_weight1

    vra2 = g.INTERP_FVRVXX_internal2_vra2
    vxa2 = g.INTERP_FVRVXX_internal2_vxa2
    Tnorma2 = g.INTERP_FVRVXX_internal2_Tnorma2
    vrb2 = g.INTERP_FVRVXX_internal2_vrb2
    vxb2 = g.INTERP_FVRVXX_internal2_vxb2
    Tnormb2 = g.INTERP_FVRVXX_internal2_Tnormb2
    weight2 = g.INTERP_FVRVXX_internal2_weight2

    nvra = Vra.size
    nvxa = Vxa.size
    nxa = Xa.size

    mu = 1

    fV = np.sqrt(Tnormb / Tnorma)

    Vtha = np.sqrt(2 * q * Tnorma / (mu * mH))
    Vtha2 = Vtha * Vtha
    Vthb = np.sqrt(2 * q * Tnormb / (mu * mH))
    Vthb2 = Vthb * Vthb

    if fa.shape != (nxa, nvxa, nvra):
        raise ValueError('Shape of fa does not match dimensions of Vra, Vxa, Xa')

    oki = np.bitwise_and(fV * Vrb <= max(Vra), fV * Vrb >= min(Vra)).nonzero()[0]
    if oki.size < 1:
        raise ValueError('No values of Vrb are within range of Vra')
    i0, i1 = oki[0], oki[-1]

    okj = np.bitwise_and(fV * Vxb <= max(Vxa), fV * Vxb >= min(Vxa)).nonzero()[0]
    if okj.size < 1:
        raise ValueError('No values of Vxb are within range of Vxa')
    j0, j1 = okj[0], okj[-1]

    okk = np.bitwise_and(Xb <= max(Xa), Xb >= min(Xa)).nonzero()[0]
    if okk.size < 1:
        raise ValueError('No values of Xb are within range of Xa')
    k0, k1 = okk[0], okk[-1]

    nvrb = Vrb.size
    nvxb = Vxb.size
    nxb = Xb.size

    fb = np.zeros((nxb, nvxb, nvrb), dtype=np.float64)

    make_dvr_dvx_out = Make_dVr_dVx(Vra, Vxa)
    Vr2pidVra, VrVr4pidVra, dVxa, vraL, vraR, vxaL, vxaR = make_dvr_dvx_out[:7]
    Vra2Vxa2 = make_dvr_dvx_out[11]
    Vr2pidVrb, VrVr4pidVrb, dVxb, vrbL, vrbR, vxbL, vxbR, Vol, Vth_DVx, Vx_DVx, Vr_DVr, Vrb2Vxb2, jpa, jpb, jna, jnb = Make_dVr_dVx(Vrb, Vxb)

    w1_active = 0
    w1_match = 0

    if vra1 is not None:
        w1_active = 1
        test = 0

        test += vra1[~np.isin(vra1, Vra)].size
        test += vxa1[~np.isin(vxa1, Vxa)].size
        test += Tnorma1[Tnorma1 != Tnorma].size
        test += Tnormb1[Tnormb1 != Tnormb].size
        test += vxb1[~np.isin(vxb1, Vxb)].size
        test += vrb1[~np.isin(vrb1, Vrb)].size

        if test == 0:
            w1_match = 1

    w2_active = 0
    w2_match = 0

    if vra2 is not None:
        w2_active = 1
        test = 0

        test += vra2[~np.isin(vra2, Vra)].size
        test += vxa2[~np.isin(vxa2, Vxa)].size
        test += Tnorma2[Tnorma2 != Tnorma].size
        test += Tnormb2[Tnormb2 != Tnormb].size
        test += vxb2[~np.isin(vxb2, Vxb)].size
        test += vrb2[~np.isin(vrb2, Vrb)].size

        if test == 0:
            w2_match = 1

    w_new = 0

    if w1_match or w2_match:
        if w1_match:
            weight = weight1
            if debug:
                print(prompt + 'using weight1')
        if w2_match:
            weight = weight2
            if debug:
                print(prompt + 'using weight2')
    else:
        if debug:
            print(prompt + 'computing new weight')
        w_new = 1

        _weight = np.zeros((nvxa, nvra, nvxb, nvrb), dtype=np.float64)
        weight = np.zeros((nvra * nvxa, nvrb * nvxb), dtype=np.float64)

        for ib in range(nvrb):
            for jb in range(nvxb):
                for ia in range(nvra):
                    vraMin = max([fV * vrbL[ib], vraL[ia]])
                    vraMax = min([fV * vrbR[ib], vraR[ia]])
                    for ja in range(nvxa):
                        vxaMin = max([fV * vxbL[jb], vxaL[ja]])
                        vxaMax = min([fV * vxbR[jb], vxaR[ja]])
                        if vraMax > vraMin and vxaMax > vxaMin:
                            _weight[ja, ia, jb, ib] = 2 * np.pi * (vraMax ** 2 - vraMin ** 2) * (vxaMax - vxaMin) / (Vr2pidVrb[ib] * dVxb[jb])

        weight = np.reshape(_weight, weight.shape)

    fb_xa = np.zeros((nxa, nvrb * nvxb), dtype=np.float64)

    _fa = np.zeros((nxa, nvra * nvxa), dtype=np.float64)
    _fa = np.reshape(fa, _fa.shape)
    fb_xa = np.matmul(_fa, weight)

    na = np.zeros(nxa, dtype=np.float64)
    _Wxa = np.zeros(nxa, dtype=np.float64)
    _Ea = np.zeros(nxa, dtype=np.float64)

    for k in range(nxa):
        na[k] = np.sum(Vr2pidVra * np.matmul(dVxa, fa[k, :, :]))
        if na[k] > 0:
            _Wxa[k] = np.sqrt(Tnorma) * np.sum(Vr2pidVra * np.matmul((Vxa * dVxa), fa[k, :, :])) / na[k]
            _Ea[k] = Tnorma * np.sum(Vr2pidVra * np.matmul(dVxa, (Vra2Vxa2 * fa[k, :, :]))) / na[k]

    wxa = np.zeros(nxb, dtype=np.float64)
    Ea = np.zeros(nxb, dtype=np.float64)

    for k in range(k0, k1 + 1):
        kL = int(np.maximum(locate(Xa, Xb[k]), 0))  # Convert to int
        kR = int(np.minimum(kL + 1, Xa.size - 1))   # Convert to int
        kL = np.minimum(kL, kR - 1)

        f = (Xb[k] - Xa[kL]) / (Xa[kR] - Xa[kL])
        fb[k, :, :] = np.reshape(fb_xa[kL, :] + (fb_xa[kR, :] - fb_xa[kL, :]) * f, fb[k, :, :].shape)
        wxa[k] = _Wxa[kL] + (_Wxa[kR] - _Wxa[kL]) * f
        Ea[k] = _Ea[kL] + (_Ea[kR] - _Ea[kL]) * f

    if correct:
        AN = np.zeros((2, nvxb, nvrb), dtype=np.float64)
        BN = np.zeros((2, nvxb, nvrb), dtype=np.float64)

        sgn = np.array([-1, 1])

        for k in range(nxb):
            allow_neg = 0

            nb = np.sum(Vr2pidVrb * np.matmul(dVxb, fb[k, :, :]))
            if nb > 0:
                goto_correct = True
                while goto_correct:
                    goto_correct = False
                    nb = np.sum(Vr2pidVrb * np.matmul(dVxb, fb[k, :, :]))
                    Wxb = np.sqrt(Tnormb) * np.sum(Vr2pidVrb * np.matmul(Vxb * dVxb, fb[k, :, :])) / nb
                    Eb = Tnormb * np.sum(Vr2pidVrb * np.matmul(dVxb, Vrb2Vxb2 * fb[k, :, :])) / nb

                    Nij = np.zeros((nvxb + 2, nvrb + 2), dtype=np.float64)
                    Nij[1:-1, 1:-1] = fb[k, :, :] * Vol / nb

                    cutoff = 1.0e-6 * np.max(Nij)
                    ii = np.argwhere(np.bitwise_and(abs(Nij) < cutoff, abs(Nij) > 0))
                    for i in ii:
                        Nij[tuple(i)] = 0
                    if max(Nij[2, :]) <= 0:
                        allow_neg = 1
                    Nijp1_vx_Dvx = np.roll(Nij * Vx_DVx, -1, 0)
                    Nij_vx_Dvx = Nij * Vx_DVx
                    Nijm1_vx_Dvx = np.roll(Nij * Vx_DVx, 1, 0)
                    Nip1j_vr_Dvr = np.roll(Nij * Vr_DVr, -1, 1)
                    Nij_vr_Dvr = Nij * Vr_DVr
                    Nim1j_vr_Dvr = np.roll(Nij * Vr_DVr, 1, 1)

                    _AN = np.roll(Nij * Vth_DVx, 1, 0) - Nij * Vth_DVx
                    AN[0, :, :] = _AN[1:nvxb + 1, 1:nvrb + 1]
                    _AN = -np.roll(Nij * Vth_DVx, -1, 0)
                    AN[1, :, :] = _AN[1:nvxb + 1, 1:nvrb + 1]

                    BN[0, jpa + 1:jpb + 1, :] = Nijm1_vx_Dvx[jpa + 2:jpb + 2, 1:nvrb + 1] - Nij_vx_Dvx[jpa + 2:jpb + 2, 1:nvrb + 1]
                    BN[0, jpa, :] = -Nij_vx_Dvx[jpa + 1, 1:nvrb + 1]
                    BN[0, jnb, :] = Nij_vx_Dvx[jnb + 1, 1:nvrb + 1]
                    BN[0, jna:jnb, :] = -Nijp1_vx_Dvx[jna + 1:jnb + 1, 1:nvrb + 1] + Nij_vx_Dvx[jna + 1:jnb + 1, 1:nvrb + 1]
                    BN[0, :, :] = BN[0, :, :] + Nim1j_vr_Dvr[1:nvxb + 1, 1:nvrb + 1] - Nij_vr_Dvr[1:nvxb + 1, 1:nvrb + 1]

                    BN[1, jpa + 1:jpb + 1, :] = -Nijp1_vx_Dvx[jpa + 2:jpb + 2, 1:nvrb + 1] + Nij_vx_Dvx[jpa + 2:jpb + 2, 1:nvrb + 1]
                    BN[1, jpa, :] = -Nijp1_vx_Dvx[jpa + 1, 1:nvrb + 1]
                    BN[1, jnb, :] = Nijm1_vx_Dvx[jnb + 1, 1:nvrb + 1]
                    BN[1, jna:jnb, :] = Nijm1_vx_Dvx[jna + 1:jnb + 1, 1:nvrb + 1] - Nij_vx_Dvx[jna + 1:jnb + 1, 1:nvrb + 1]
                    BN[1, :, 1:nvrb] = BN[1, :, 1:nvrb] - Nip1j_vr_Dvr[1:nvxb + 1, 2:nvrb + 1] + Nij_vr_Dvr[1:nvxb + 1, 2:nvrb + 1]
                    BN[1, :, 0] = BN[1, :, 0] - Nip1j_vr_Dvr[1:nvxb + 1, 1]

                    if allow_neg:
                        BN[1, :, 0] = BN[1, :, 0] - Nij_vr_Dvr[1:nvxb + 1, 1]
                        BN[1, :, 1] = BN[1, :, 1] + Nij_vr_Dvr[1:nvxb + 1, 1]

                    Nij = Nij[1:nvxb + 1, 1:nvrb + 1]

                    TB1 = np.zeros(2, dtype=np.float64)
                    TB2 = np.zeros(2, dtype=np.float64)

                    for ia in range(2):
                        TA1 = np.sqrt(Tnormb) * np.sum(np.matmul(Vxb, AN[ia, :, :]))
                        TA2 = Tnormb * np.sum(Vrb2Vxb2 * AN[ia, :, :])
                        for ib in range(2):
                            if TB1[ib] == 0:
                                TB1[ib] = np.sqrt(Tnormb) * np.sum(np.matmul(Vxb, BN[ib, :, :]))
                            if TB2[ib] == 0:
                                TB2[ib] = Tnormb * np.sum(Vrb2Vxb2 * BN[ib, :, :])

                            denom = TA2 * TB1[ib] - TA1 * TB2[ib]
                            beta = 0
                            alpha = 0

                            if denom != 0 and TA1 != 0:
                                beta = (TA2 * (wxa[k] - Wxb) - TA1 * (Ea[k] - Eb)) / denom
                                alpha = (wxa[k] - Wxb - TB1[ib] * beta) / TA1

                            do_break = alpha * sgn[ia] > 0 and beta * sgn[ib] > 0
                            if do_break:
                                break
                        if do_break:
                            break

                    RHS = AN[ia, :, :] * alpha + BN[ib, :, :] * beta

                    s = 1
                    if not allow_neg:
                        ii = np.nonzero(Nij)
                        if np.size(ii) > 0:
                            s = min(1 / np.max(-RHS[ii] / Nij[ii]), 1)
                    fb[k, :, :] = nb * (Nij + s * RHS) / Vol

                    goto_correct = s < 1

    if do_warn is not None:
        big = np.max(fb)

        if i0 > 0 or i1 < nvrb - 1:
            for k in range(k0, k1 + 1):
                for j in range(j0, j1 + 1):
                    if i0 > 0 and fb[k, j, i0] > do_warn * big:
                        warn('Non-zero value of fb detected at min(Vra) boundary')
                    if i1 < nvrb - 1 and fb[k, j, i1] > do_warn * big:
                        warn('Non-zero value of fb detected at max(Vra) boundary')

        if j0 > 0 or j1 < nvxb - 1:
            for k in range(k0, k1 + 1):
                for i in range(i0, i1 + 1):
                    if j0 > 0 and fb[k, j0, i] > do_warn * big:
                        warn('Non-zero value of fb detected at min(Vxa) boundary')
                    if j1 < nvxb - 1 and fb[k, j1, i] > do_warn * big:
                        warn('Non-zero value of fb detected at max(Vxa) boundary')

        if k0 > 0 or k1 < nxb - 1:
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    if k0 > 0 and fb[k0, j, i] > do_warn * big:
                        warn('Non-zero value of fb detected at min(Xa) boundary')
                    if k1 < nxb - 1 and fb[k1, j, i] > do_warn * big:
                        warn('Non-zero value of fb detected at max(Xa) boundary')

    tot_a = np.zeros(nxa, dtype=np.float64)
    for k in range(nxa):
        tot_a[k] = np.sum(Vr2pidVra * np.matmul(dVxa, fa[k, :, :]))
    tot_b = np.zeros(nxb, dtype=np.float64)
    tot_b[k0:k1 + 1] = interpolate.interp1d(Xa, tot_a, fill_value="extrapolate")(Xb[k0:k1 + 1])
    ii = np.argwhere(fb > 0)
    if ii.size > 0:
        min_tot = np.min(np.array([fb[tuple(i)] for i in ii]))
        for k in range(k0, k1 + 1):
            tot = np.sum(Vr2pidVrb * np.matmul(dVxb, fb[k, :, :]))
            if tot > min_tot:
                if debug:
                    print(prompt + 'Density renormalization factor =' + sval(tot_b[k] / tot))

    if debug:
        na = np.zeros(nxa, dtype=np.float64)
        Uxa = np.zeros(nxa, dtype=np.float64)
        Ta = np.zeros(nxa, dtype=np.float64)
        vr2vx2_ran2 = np.zeros((nvxa, nvra), dtype=np.float64)

        for k in range(nxa):
            na[k] = np.sum(Vr2pidVra * np.matmul(dVxa, fa[k, :, :]))
            if na[k] > 0:
                Uxa[k] = Vtha * np.sum(Vr2pidVra * np.matmul(Vxa * dVxa, fa[k, :, :])) / na[k]
                for i in range(nvra):
                    vr2vx2_ran2[:, i] = Vra[i] ** 2 + (Vxa - Uxa[k] / Vtha) ** 2
                Ta[k] = mu * mH * Vtha2 * np.sum(Vr2pidVra * np.matmul(dVxa, vr2vx2_ran2 * fa[k, :, :])) / (3 * q * na[k])

        nb = np.zeros(nxb, dtype=np.float64)
        Uxb = np.zeros(nxb, dtype=np.float64)
        Tb = np.zeros(nxb, dtype=np.float64)
        vr2vx2_ran2 = np.zeros((nvxb, nvrb), dtype=np.float64)

        for k in range(nxb):
            nb[k] = np.sum(Vr2pidVrb * np.matmul(dVxb, fb[k, :, :]))
            if nb[k] > 0:
                Uxb[k] = Vthb * np.sum(Vr2pidVrb * np.matmul(Vxb * dVxb, fb[k, :, :])) / nb[k]
                for i in range(nvrb):
                    vr2vx2_ran2[:, i] = Vrb[i] ** 2 + (Vxb - Uxb[k] / Vthb) ** 2
                Tb[k] = mu * mH * Vthb * np.sum(Vr2pidVrb * np.matmul(dVxb, vr2vx2_ran2 * fb[k, :, :])) / (3 * q * nb[k])

    if w_new:
        if w1_active:
            if debug:
                print(prompt + 'Storing Weight in Weight2')
            g.INTERP_FVRVXX_internal2_vra2 = Vra
            g.INTERP_FVRVXX_internal2_vxa2 = Vxa
            g.INTERP_FVRVXX_internal2_Tnorma2 = Tnorma
            g.INTERP_FVRVXX_internal2_vrb2 = Vrb
            g.INTERP_FVRVXX_internal2_vxb2 = Vxb
            g.INTERP_FVRVXX_internal2_Tnormb2 = Tnormb
            g.INTERP_FVRVXX_internal2_weight2 = weight
        else:
            if debug:
                print(prompt + 'Storing Weight in Weight1')
            g.INTERP_FVRVXX_internal1_vra1 = Vra
            g.INTERP_FVRVXX_internal1_vxa1 = Vxa
            g.INTERP_FVRVXX_internal1_Tnorma1 = Tnorma
            g.INTERP_FVRVXX_internal1_vrb1 = Vrb
            g.INTERP_FVRVXX_internal1_vxb1 = Vxb
            g.INTERP_FVRVXX_internal1_Tnormb1 = Tnormb
            g.INTERP_FVRVXX_internal1_weight1 = weight

    return fb

