import copy
import numpy               as np
import matplotlib.pyplot   as plt
from mpl_toolkits.mplot3d  import Axes3D
from scipy.interpolate     import interp1d

from p000_variables        import VARIABLES
from p001_Make_dVr_dVx     import Make_dVr_dVx
from p028_locate           import locate

from p000_common_variables import obj_interp_fvrvxx_1
from p000_common_variables import obj_interp_fvrvxx_2

from numba import njit
from numba import prange
from tqdm import tqdm





@njit(parallel=True)
def compute_weight(_weight, f_v, vrL_b, vrR_b, vxL_b, vxR_b,
                vrL_a, vrR_a, vxL_a, vxR_a, Vr2pidVr_b, dVx_b,
                nvr_b, nvx_b, nvr_a, nvx_a):

    for ib in prange(nvr_b):
        for jb in prange(nvx_b):
            for ia in prange(nvr_a):
                vraMin = max(f_v * vrL_b[ib], vrL_a[ia])
                vraMax = min(f_v * vrR_b[ib], vrR_a[ia])
                for ja in prange(nvx_a):
                    vxaMin = max(f_v * vxL_b[jb], vxL_a[ja])
                    vxaMax = min(f_v * vxR_b[jb], vxR_a[ja])

                    if (vraMax > vraMin) and (vxaMax > vxaMin):
                        _weight[ib, jb, ia, ja] = (
                            2 * np.pi * (vraMax**2 - vraMin**2) * (vxaMax - vxaMin)
                        ) / (Vr2pidVr_b[ib] * dVx_b[jb])

    return _weight


def interp_f_vr_vx_X(
        f_a :    np.ndarray,        # Distrib. Function: f_a = np.array((len(vr_a),len(vx_a),len(x_a)))
        vr_a:    np.ndarray,        # Radial velocity
        vx_a:    np.ndarray,        # Axial  velocity
        x_a :    np.ndarray,        # Spatial coordinate
        Tnorm_a: float,

        #f_b:    np.ndarray,        # Distrib. Function: f_b = np.array((len(vr_b),len(vx_b),len(x_b)))
        vr_b:    np.ndarray,        # Radial velocity
        vx_b:    np.ndarray,        # Axial  velocity
        x_b :    np.ndarray,        # Spatial coordinate
        Tnorm_b: float,

        debug   = 0,    
        correct = 1,
        warn: float = None # float, acceptable truncation level.
    ):
    
    nvr_a = len(vr_a)
    nvx_a = len(vx_a)
    nx_a  = len(x_a)

    #print(f'f_a.shape: {f_a.shape}')
    #print(f'nvr_a: {nvr_a}')    
    #print(f'nvx_a: {nvx_a}')
    #print(f'nx_a:  {nx_a}')

    mH = VARIABLES.mH
    q  = VARIABLES.q

    mu = 1
    #--------------------------------------------------------------------------------
    f_v = np.sqrt(Tnorm_b/Tnorm_a)
    # Compute vth_a, vth_a2, vth_b, vth_b2
    vth_a  = np.sqrt((2*q*Tnorm_a)/(mu*mH))
    vth_a2 = vth_a**2
    
    vth_b  = np.sqrt((2*q*Tnorm_b)/(mu*mH))
    vth_b2 = vth_b**2

    #--------------------------------------------------------------------------------
    # Comprobar si las dimensiones de fa coinciden con las de Vra, Vxa y Xa
    if f_a.shape[0] != nvr_a:
        raise ValueError("Number of elements in fa(*,0,0) and vr_a do not agree!")

    if f_a.shape[1] != nvx_a:
        raise ValueError("Number of elements in fa(0,*,0) and vx_a do not agree!")

    if f_a.shape[2] != nx_a:
        raise ValueError("Number of elements in fa(0,0,*) and x_a  do not agree!")
    
    #--------------------------------------------------------------------------------
    # Verificar valores de vr_b dentro del rango de vr_a
    oki = np.where((f_v * vr_b <= np.nanmax(vr_a)) & (f_v * vr_b >= np.nanmin(vr_a)))[0]
    if len(oki) < 1:
        raise ValueError("No values of vr_b are within range of vr_a")
    i0, i1 = oki[0], oki[-1]

    # Verificar valores de vx_b dentro del rango de vx_a
    okj = np.where((f_v * vx_b <= np.nanmax(vx_a)) & (f_v * vx_b >= np.nanmin(vx_a)))[0]
    if len(okj) < 1:
        raise ValueError("No values of vx_b are within range of vx_a")
    j0, j1 = okj[0], okj[-1]

    # Verificar valores de Xb dentro del rango de Xa
    okk = np.where((x_b <= np.nanmax(x_a)) & (x_b >= np.nanmin(x_a)))[0]
    #print(okk)
    if len(okk) < 1:
        raise ValueError("No values of x_b are within range of x_a")
    k0, k1 = okk[0], okk[-1]
    #print('okk',okk)
    #print('k0',k0,'k1',k1)
    #print('x_b',x_b)
    #print('here_3')

    nvr_b = len(vr_b)
    nvx_b = len(vx_b)
    nx_b  = len(x_b)

    f_b = np.zeros((nvr_b, nvx_b, nx_b), dtype=np.float64)
    #print(f'f_b.shape: {f_b.shape}')
    #print(f'nvr_b: {nvr_b}')    
    #print(f'nvx_b: {nvx_b}')
    #print(f'nx_b:  {nx_b}')

    #print('makedvrdvx_a')
    Vr2pidVr_a, VrVr4pidVr_a, dVx_a, vrL_a, vrR_a, vxL_a, vxR_a,\
    vol_a, vth_Deltavx_a, vx_Deltavx_a, vr_Deltavr_a, vr2vx2_a,\
    jpa_a, jpb_a, jna_a, jnb_a = Make_dVr_dVx(vr_a,vx_a)
    #print('makedvrdvx_b')
    Vr2pidVr_b, VrVr4pidVr_b, dVx_b, vrL_b, vrR_b, vxL_b, vxR_b,\
    vol_b, vth_Deltavx_b, vx_Deltavx_b, vr_Deltavr_b, vr2vx2_b,\
    jpa_b, jpb_b, jna_b, jnb_b = Make_dVr_dVx(vr_b,vx_b)

    w1_active = 0
    w1_match  = 0

    if obj_interp_fvrvxx_1.vra1 is not None:
        w1_active = 1
        test = 0
        # ii = len(np.where(obj_interp_fvrvxx_1.vra1 != vr_a)[0])
        ii = np.nansum(obj_interp_fvrvxx_1.vra1    != vr_a     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_1.vxa1    != vx_a     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_1.Tnorma1 != Tnorm_a  )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_1.vrb1    != vr_b     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_1.vxb1    != vx_b     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_1.Tnormb1 != Tnorm_b  )
        test += ii
        if test <= 0:
            w1_match = 1

    w2_active = 0
    w2_match  = 0

    if obj_interp_fvrvxx_2.vra2 is not None:
        w2_active = 1
        test = 0
        ii = np.nansum(obj_interp_fvrvxx_2.vra2    != vr_a     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_2.vxa2    != vx_a     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_2.Tnorma2 != Tnorm_a  )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_2.vrb2    != vr_b     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_2.vxb2    != vx_b     )
        test += ii
        ii = np.nansum(obj_interp_fvrvxx_2.Tnormb2 != Tnorm_b  )
        test += ii
        if test <= 0:
            w2_match = 1

    w_new = 0
    if w1_match or w2_match:
        if w1_match:
            weight = obj_interp_fvrvxx_1.weight1
            if debug:
                print('using Weight1')
        if w2_match:
            weight = obj_interp_fvrvxx_2.weight2
            if debug:
                print('using Weight2')
    
    else:
        w_new = 1        
        if debug:
            print('Computing new Weight')

        _weight = np.zeros((nvr_b,  nvx_b, nvr_a,  nvx_a), dtype=np.float64)
        weight  = np.zeros((nvr_b * nvx_b, nvr_a * nvx_a), dtype=np.float64)

        # _weight = np.zeros((nvr_b, nvx_b, nvr_a, nvx_a), dtype=np.float64, order='F')
        # weight  = _weight.reshape((nvr_b * nvx_b, nvr_a * nvx_a), order='F')
        
        _weight = compute_weight(_weight, f_v, vrL_b, vrR_b, vxL_b, vxR_b,
                                vrL_a, vrR_a, vxL_a, vxR_a, Vr2pidVr_b, dVx_b,
                                nvr_b, nvx_b, nvr_a, nvx_a)
        # for ib in tqdm(range(nvr_b)):
        #     for jb in tqdm(range(nvx_b)):
        #         for ia in range(nvr_a):
        #             vraMin = max(f_v * vrL_b[ib], vrL_a[ia])
        #             vraMax = min(f_v * vrR_b[ib], vrR_a[ia])
        #             for ja in range(nvx_a):
        #                 vxaMin = max(f_v * vxL_b[jb], vxL_a[ja])
        #                 vxaMax = min(f_v * vxR_b[jb], vxR_a[ja])

        #                 if (vraMax > vraMin) and (vxaMax> vxaMin):
        #                     _weight[ib, jb, ia, ja] = (2 * np.pi * (vraMax**2 - vraMin**2) * (vxaMax - vxaMin)) / (Vr2pidVr_b[ib] * dVx_b[jb])

        # weight = np.transpose(_weight).reshape((nvr_b * nvx_b, nvr_a * nvx_a))
        weight = _weight.reshape((nvr_b * nvx_b, nvr_a * nvx_a), order='F')
        #print('weight:',weight.shape)
        
    # Determine fb_xa from weight array
    fb_xa  = np.zeros((                 nvr_b*nvx_b, nx_a))
    #print(f'fb_xa.shape01:{fb_xa.shape}')
    _fa    = np.zeros((                 nvr_a*nvx_a, nx_a))
    #print(f'f_a-------------------------------------------: {f_a.shape}')
    # _fa    = np.transpose(f_a).reshape((nvr_a*nvx_a, nx_a))
    _fa = f_a.reshape((nvr_a * nvx_a, nx_a), order='F')

    # _fa    = f_a.reshape((nvr_a*nvx_a, nx_a))
    #print('weight:',weight.shape)
    #print('_fa:',_fa.shape)
    fb_xa  = np.matmul(weight, _fa)
    #print(f'fb_xa.shape02:{fb_xa.shape}')
    
    # Compute _Wxa and _Ea - these are the desired moments of fb, but on the xa grid
    na     = np.zeros(nx_a)
    _Wxa   = np.zeros(nx_a)
    _Ea    = np.zeros(nx_a)

    for k in range(nx_a):
        # #print(f'k:{k}')
        na[k] = np.nansum(Vr2pidVr_a*(np.matmul(f_a[:,:,k],dVx_a)))
        if na[k] > 0:
            _Wxa[k] = np.sqrt(Tnorm_a) * np.nansum(Vr2pidVr_a * np.matmul(f_a[:, :, k], (vx_a * dVx_a))) / na[k]
            _Ea[k] = Tnorm_a *   np.nansum(Vr2pidVr_a * np.matmul((vr2vx2_a[:, :] * f_a[:, :, k]), dVx_a)) / na[k]




    # Interpolate in x to get fb from fb_xa and to get Wxa, Ea from _Wxa, _Ea
    Wxa    = np.zeros(nx_b)
    Ea     = np.zeros(nx_b)
    #print('#'*100)
    #print(k0, k1)
    for k in range(k0,k1+1):
        kL = max(locate(x_a,x_b[k])[0],0)
        kR = min(kL + 1, nx_a - 1)
        kL = min(kL,     kR   - 1)
        f  = (x_b[k] - x_a[kL]) / (x_a[kR] - x_a[kL])
        # #print('f_b[:,:,k]:',f_b[:,:,k].shape)
        aux_f_b = fb_xa[:,kL] + ((fb_xa[:,kR] - fb_xa[:,kL]) * f)
        #print('aux_f_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb:',aux_f_b.shape)
        # f_b[:,:,k] = np.transpose(aux_f_b).reshape((nvr_b, nvx_b))
        # f_b[:,:,k] = aux_f_b.reshape((nvr_b, nvx_b), order='F')
        f_b[:, :, k] = aux_f_b.reshape((nvr_b, nvx_b), order='F')
        #print('f_b:',f_b.shape)

        Wxa[k] = _Wxa[kL] + ((_Wxa[kR] - _Wxa[kL]) * f)
        Ea[k]  =  _Ea[kL] + (( _Ea[kR] -  _Ea[kL]) * f)



    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(f_b.shape[0]):
        x = x_b
        y = np.full_like(x_b, i)  # eje Y: índice i constante
        z = f_b[i, 0, :]
        ax.plot(x, y, z, color='b', alpha=0.5)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    for i in range(f_b.shape[0]):
        x = vx_b
        y = np.full_like(vx_b, i)  # eje Y: índice i constante
        z = f_b[i, :, 0]
        ax.plot(x, y, z, color='b', alpha=0.5)
    plt.tight_layout()
    plt.show()

    #print("f_b[:,0,0] =", f_b[:,0,0])
    #print("f_b[0,:,0] =", f_b[0,:,0])
    #print("f_b[0,0,:] =", f_b[0,0,:])
    #print("Hay NaNs?", np.isnan(f_b).any())
    #print("Hay Infs?", np.isinf(f_b).any())
    #print("Min/Max f_b[:,0,0]:", np.nanmin(f_b[:,0,0]), np.nanmax(f_b[:,0,0]))
    #print('#'*70)
    #print('vr_b',vr_b)
    #print('vx_b',vx_b)
    #print('#'*70)
    #print('comparación entre fa y fb')
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))  # 1 fila, 3 columnas
    # Primer scatter
    # axs[0].scatter(vr_b, f_b[:, 0, 0])
    axs[0].plot(vr_b, f_b[:, 0, 0])
    axs[0].set_title('vr_b vs f_b[:,0,0]')
    axs[0].grid(True)

    # Segundo scatter
    # axs[1].scatter(vx_b, f_b[0, :, 0])
    axs[1].plot(vx_b, f_b[0, :, 0])
    axs[1].set_title('vx_b vs f_b[0,:,0]')
    axs[1].grid(True)

    # Tercer scatter
    axs[2].scatter(x_b, f_b[0, 0, :])
    axs[2].plot(x_b, f_b[0, 0, :])
    axs[2].set_title('x_b vs f_b[0,0,:]')
    axs[2].grid(True)

    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))  # figsize ajusta el tamaño total del gráfico

    # Gráfico 1: x_a vs na
    axes[0].plot(x_a, na, color='blue', label='na')
    axes[0].set_title('x_a vs na')
    axes[0].set_xlabel('x_a')
    axes[0].set_ylabel('na')
    axes[0].grid(True)
    axes[0].legend()

    # Gráfico 2: x_a vs _Wxa
    axes[1].plot(x_a, _Wxa, color='green', label='xa_Wxa')
    axes[1].plot(x_b,  Wxa, color='blue', label='xb_Wxa')
    axes[1].set_title('x_a vs _Wxa')
    axes[1].set_xlabel('x_a')
    axes[1].set_ylabel('_Wxa')
    axes[1].grid(True)
    axes[1].legend()

    # Gráfico 3: x_a vs _Ea
    axes[2].plot(x_a, _Ea, color='red', label='xa_Ea')
    axes[2].plot(x_b,  Ea, color='blue', label='xb_Ea')
    axes[2].set_title('x_a vs _Ea')
    axes[2].set_xlabel('x_a')
    axes[2].set_ylabel('_Ea')
    axes[2].grid(True)
    axes[2].legend()
    # Ajustar el diseño para evitar superposiciones
    plt.tight_layout()
    plt.show()

    #print(f'f_b[:, 0, 0]:{f_b[:, 0, 0]}')
    #print(f'f_b[0, :, 0]:{f_b[0, :, 0]}')
    #print(f'f_b[0, 0, :]:{f_b[0, 0, :]}')
    #print('pre_if')
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))  # 1 fila, 3 columnas
    # Primer scatter
    axs[0].scatter(vr_b, f_b[:, 0, 0])
    axs[0].plot(vr_b, f_b[:, 0, 0])
    axs[0].set_title('021vr_b vs f_b[:,0,0]')
    axs[0].grid(True)

    # Segundo scatter
    axs[1].scatter(vx_b, f_b[0, :, 0])
    axs[1].plot(vx_b, f_b[0, :, 0])
    axs[1].set_title('021vx_b vs f_b[0,:,0]')
    axs[1].grid(True)

    # Tercer scatter
    axs[2].scatter(x_b, f_b[0, 0, :])
    axs[2].plot(x_b, f_b[0, 0, :])
    axs[2].set_title('021x_b vs f_b[0,0,:]')
    axs[2].grid(True)

    plt.tight_layout()
    plt.show() 
    if correct:
        #print('Dentro de correct')
        AN = np.zeros((nvr_b,nvx_b,2))
        BN = np.zeros((nvr_b,nvx_b,2))
        sgn = [1,-1]
        
        for k in range(0,nx_b):
            # #print('dentro del range k de 0 a nx_b')
            s=0
            allow_neg = False
            #  Compute nb, Wxb, and Eb - these are the current moments of fb}
            nb = np.nansum(Vr2pidVr_b*(np.matmul(f_b[:,:,k],dVx_b)))
            #print('nb>0:',nb>0,'nb:',nb)
            if nb > 0:
                exit_for_ib = False
                contador = 0
                while s<1:
                    contador = contador +1
                    #print(f'Dentro del while, k:{k},\t contador del for:{contador}')
                    # Compute Wxb, and Eb - these are the current moments of fb
                    nb  = np.nansum(Vr2pidVr_b*(np.matmul(f_b[:,:,k],dVx_b)))                       
                    Wxb = np.sqrt(Tnorm_b) * np.nansum(Vr2pidVr_b * (np.matmul(f_b[:, :, k], (vx_b * dVx_b)))) / nb
                    Eb  = Tnorm_b * np.nansum(Vr2pidVr_b * ((np.matmul((vr2vx2_b * f_b[:, :, k]), dVx_b)))) / nb
                    
                    # Compute Nij from fb, padded with zeros, this is used as a normalization
                    Nij = np.zeros((nvr_b+2,nvx_b+2))
                    Nij[1:nvr_b+1, 1:nvx_b+1] = f_b[:,:,k] * vol_b/nb
                    #print(f'Nij: {Nij}')
                    # Set Cutoff and remove Nij very close to zero
                    cutoff = 1.0e-6 * np.nanmax(Nij)
                    mask = (np.abs(Nij) < cutoff) & (np.abs(Nij) > 0.0)
                    Nij[mask] = 0.0
                    #print(f'Nij: {Nij}')
                    # ii = np.where((np.abs(Nij) < cutoff) & (np.abs(Nij) > 0.0))[0]
                    # ii = np.where((np.abs(Nij) < cutoff) & (np.abs(Nij) > 0.0))
                    # #print('ii_1',ii)
                    # if len(ii) > 0:
                    #     Nij[ii] = 0.0
                    ii = np.where(Nij[2,:] > 0.0)[0]
                    # #print('ii_2',len(ii))
                    if len(ii) < 1:
                        allow_neg = True

                    auxNij        = Nij*vx_Deltavx_b
                    Nijp1_vx_Dvx  = np.roll(auxNij,  shift=-1, axis=1)
                    Nij_vx_Dvx    = Nij*vx_Deltavx_b
                    Nijm1_vx_Dvx  = np.roll(auxNij,  shift= 1, axis=1)
                    auxNij2       = Nij*vr_Deltavr_b
                    Nip1j_vr_Dvr  = np.roll(auxNij2, shift=-1, axis=0)
                    Nij_vr_Dvr    = Nij*vr_Deltavr_b
                    Nim1j_vr_Dvr  = np.roll(auxNij2, shift= 1, axis=0)

                    # Compute Ap, Am, Bp, and Bm (0=p 1=m)
                    aux_AN        = np.zeros((nvx_b+2,nvr_b+2), dtype=np.float64)
                    aux_AN        = Nij*vth_Deltavx_b
                    
                    _AN           = np.zeros((nvx_b+2,nvr_b+2), dtype=np.float64)
                    _AN           = np.roll(aux_AN, shift=1, axis=1) - aux_AN
                    AN[:,:,0]     = copy.copy(_AN[1:nvr_b+1,1:nvx_b+1])                    
                    
                    _AN           = np.zeros((nvx_b+2,nvr_b+2), dtype=np.float64)
                    _AN           =-np.roll(aux_AN, shift=-1, axis=1) + aux_AN
                    AN[:,:,1]     = copy.copy(_AN[1:nvr_b+1,1:nvx_b+1])

                    BN[:,jpa_b+1:jpb_b+1,0] =  Nijm1_vx_Dvx[1:nvr_b+1,jpa_b+2:jpb_b+2] - Nij_vx_Dvx[1:nvr_b+1,jpa_b+2:jpb_b+2]
                    BN[:,jpa_b,0]           = -Nij_vx_Dvx[1:nvr_b+1,jpa_b+1]
                    BN[:,jnb_b,0]           =  Nij_vx_Dvx[1:nvr_b+1,jnb_b+1]
                    BN[:,jna_b:jnb_b,0]     = -Nijp1_vx_Dvx[1:nvr_b+1,jna_b+1:jnb_b+1] + Nij_vx_Dvx[1:nvr_b+1,jna_b+1:jnb_b+1]
                    BN[:,:,0]               =  BN[:,:,0] + Nim1j_vr_Dvr[1:nvr_b+1,1:nvx_b+1] - Nij_vr_Dvr[1:nvr_b+1,1:nvx_b+1]

                    BN[:,jpa_b+1:jpb_b+1,1] = -Nijp1_vx_Dvx[1:nvr_b+1,jpa_b+2:jpb_b+2] + Nij_vx_Dvx[1:nvr_b+1,jpa_b+2:jpb_b+2]
                    BN[:,jpa_b,1]           = -Nijp1_vx_Dvx[1:nvr_b+1,jpa_b+1]
                    BN[:,jnb_b,1]           =  Nijm1_vx_Dvx[1:nvr_b+1,jnb_b+1]
                    BN[:,jna_b:jnb_b,1]     =  Nijm1_vx_Dvx[1:nvr_b+1,jna_b+1:jnb_b+1] - Nij_vx_Dvx[1:nvr_b+1,jna_b+1:jnb_b+1]
                    BN[1:nvr_b,:,1]         =  BN[1:nvr_b,:,1] - Nip1j_vr_Dvr[2:nvr_b+1,1:nvx_b+1] + Nij_vr_Dvr[2:nvr_b+1,1:nvx_b+1]
                    BN[0,:,1]               =  BN[0,:,1] - Nip1j_vr_Dvr[1,1:nvx_b+1]

                    #  If negative values for Nij must be allowed, then add postive particles to i=0
                    #  and negative particles to i=1 (beta is negative here)
                    if allow_neg:
                        BN[0,:,1] = BN[0,:,1] - Nij_vr_Dvr[1,1:nvx_b+1]
                        BN[1,:,1] = BN[1,:,1] - Nij_vr_Dvr[1,1:nvx_b+1]
                    
                    # Remove padded zeros in Nij
                    Nij = Nij[1:nvr_b+1,1:nvx_b+1]
                    # Compute TA1, TA2
                    TB1 = np.zeros(2)
                    TB2 = np.zeros(2)
                    
                    try:
                        #print('dentro del try')
                        for ia in range(0,2):
                            TA1 = np.sqrt(Tnorm_b)*np.nansum(np.matmul(AN[:,:,ia],vx_b))
                            # #print('TA1',TA1)
                            TA2 = Tnorm_b*np.nansum(vr2vx2_b*AN[:,:,ia])
                            for ib in range(0,2):
                                # Compute TB1, TB2
                                if TB1[ib] == 0:
                                    TB1[ib] = np.sqrt(Tnorm_b)*np.nansum(np.matmul(BN[:,:,ib],vx_b))
                                if TB2[ib] == 0:
                                    TB2[ib] = Tnorm_b * np.nansum(vr2vx2_b*BN[:,:,ib])
                                denom = (TA2*TB1[ib]) - (TA1*TB2[ib])
                                beta  = 0.0
                                alpha = 0.0
                                # #print('denom', denom)
                                if (denom != 0.0) and (TA1 != 0.0):
                                    beta  = ( (TA2*(Wxa[k]-Wxb)) - (TA1*(Ea[k]-Eb) )) / denom
                                    alpha = ( Wxa[k]-Wxb-TB1[ib]*beta ) / TA1
                                #print(alpha*sgn[ia] > 0.0,'\t',alpha*sgn[ia])
                                #print(beta*sgn[ia]  > 0.0,'\t',beta*sgn[ia] )
                                #print('cond.: ',(alpha*sgn[ia] > 0.0) and (beta*sgn[ia] > 0.0))
                                if (alpha*sgn[ia] > 0.0) and (beta*sgn[ia] > 0.0):
                                    raise StopIteration
                                # #print('ib:',ib)
                                # #print('ia:',ia)
                                if ib==1 and ia==1:
                                    #print(f'salimos del bucle ib: {ib}')
                                    exit_for_ib = True
                                    break
                            if exit_for_ib:
                                #print(f'salimos del bucle ia: {ia}')
                                break
                        if exit_for_ib:
                            #print(f'salimos del try')
                            break

                    except StopIteration:
                        #print('dentro del excepTTTTTTTTTTTTTTTTTTTT')
                        RHS = (AN[:,:,ia]*alpha) + (BN[:,:,ib]*beta)
                        ii = np.where((Nij == 0.0) & (RHS < 0.0))[0]
                        if len(ii) > 0:
                            print('There are locations where Nij = 0.0 and RHS is negative')
                        s = 1.0
                        #print('ALLOW_NEGATIVEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE',allow_neg)
    
                        if not allow_neg:
                            ii = np.where(Nij != 0.0)[0]
                            if len(ii) > 0:
                                s = min(1.0/np.nanmax(-RHS[ii]/Nij[ii]),1.0)
                        f_b[:,:,k] = nb * (Nij + s*RHS)/vol_b
                        print(f'f_b[:, 0, 0]:{f_b[:, 0, 0]}')
                        print(f'f_b[0, :, 0]:{f_b[0, :, 0]}')
                        print(f'f_b[0, 0, :]:{f_b[0, 0, :]}')
                if exit_for_ib:
                    print('Los bucles terminaron y salimos del try-except y del while')


    fig, axs = plt.subplots(1, 3, figsize=(15, 4))  # 1 fila, 3 columnas
    # Primer scatter
    axs[0].scatter(vr_b, f_b[:, 0, 0])
    axs[0].plot(vr_b, f_b[:, 0, 0])
    axs[0].set_title('021vr_b vs f_b[:,0,0]')
    axs[0].grid(True)

    # Segundo scatter
    axs[1].scatter(vx_b, f_b[0, :, 0])
    axs[1].plot(vx_b, f_b[0, :, 0])
    axs[1].set_title('021vx_b vs f_b[0,:,0]')
    axs[1].grid(True)

    # Tercer scatter
    axs[2].scatter(x_b, f_b[0, 0, :])
    axs[2].plot(x_b, f_b[0, 0, :])
    axs[2].set_title('021x_b vs f_b[0,0,:]')
    axs[2].grid(True)

    plt.tight_layout()
    plt.show() 

    if warn is not None:   # If warn is something, any value but None
        # i0  & i1
        big = np.nanmax(f_b)
        i0_error = 0
        i1_error = 0
        if (i0 >0) or (i1 < nvr_b-1):
            #print(f'i0:{i0}, {i0>0} and i1:{i1},{i1 < nvr_b-1}')
            # #print(f'k0,k1+1:{k0,k1+1}')
            for k in range(k0,k1+1):
                # #print(f'j0,j1+1:{j0,j1+1}')
                for j in range(j0,j1+1):
                    # #print(f'(i0_error == 0) and (i0 > 0)       and f_b[i0,j,k] > warn*big:{(i0_error == 0) and (i0 > 0)       and f_b[i0,j,k] > warn*big}')
                    if (i0_error == 0) and (i0 > 0)       and f_b[i0,j,k] > warn*big:
                        #print('Non-zero value of fb detected at min(vr_a) boundary')
                        i0_error = 1
                    # #print(f'(i1_error == 0) and (i1 < nvr_b-1) and f_b[i1,j,k] > warn*big:{(i1_error == 0) and (i1 < nvr_b-1) and f_b[i1,j,k] > warn*big}')
                    if (i1_error == 0) and (i1 < nvr_b-1) and f_b[i1,j,k] > warn*big:
                        #print('Non-zero value of fb detected at max(vr_a) boundary')
                        i1_error = 1
                    
        # j0  & j1
        j0_error = 0
        j1_error = 0        
        if (j0 > 0) or (j1 < nvx_b - 1):
            for k in range(k0, k1 + 1):  # Itera sobre k
                for i in range(i0, i1 + 1):  # Itera sobre i
                    if (j0_error == 0) and (j0 > 0)         and (f_b[i, j0, k] > warn * big):
                        #print('Non-zero value of fb detected at min(vx_a) boundary')
                        j0_error = 1
                    if (j1_error == 0) and (j1 < nvx_b - 1) and (f_b[i, j1, k] > warn * big):
                        #print('Non-zero value of fb detected at max(vx_a) boundary')
                        j1_error = 1

        # k0  & k1
        k0_error = 0
        k1_error = 0
        if (k0 > 0) or (k1 < nx_b - 1):
            for i in range(i0, i1 + 1):  # Itera sobre i
                for j in range(j0, j1 + 1):  # Itera sobre j
                    if (k0_error == 0) and (k0 > 0) and (f_b[i, j, k0] > warn * big):
                        #print('Non-zero value of fb detected at min(x_a) boundary')
                        k0_error = 1
                    if (k1_error == 0) and (k1 < nx_b - 1) and (f_b[i, j, k1] > warn * big):
                        #print('Non-zero value of fb detected at max(x_a) boundary')
                        k1_error = 1

    # Re-scale
    tot_a = np.zeros(nx_a)
    tot_b = np.zeros(nx_b)
    for k in range(nx_a):
        tot_a[k] = np.nansum(Vr2pidVr_a*(np.matmul(f_a[:,:,k],dVx_a)))
    
    
    interpolator = interp1d(x_a, tot_a, kind='linear', fill_value="extrapolate")
    tot_b[k0:k1 + 1] = interpolator(x_b[k0:k1 + 1])


    plt.plot(x_a,tot_a,linewidth=1)
    plt.plot(x_b,tot_b,linewidth=2)
    plt.title('xa vs total')
    plt.grid(True)
    plt.show()

    ii = np.where(f_b > 0.0)
    valores_positivos = f_b[ii]

    #print(f'f_b[:, 0, 0]:{f_b[:, 0, 0]}')
    #print(f'f_b[0, :, 0]:{f_b[0, :, 0]}')
    #print(f'f_b[0, 0, :]:{f_b[0, 0, :]}')
    #print('f_b.shape:',f_b.shape)
    #print('$'*88)
    #print('$'*88)
    #print('$'*88)
    #print('$'*88)
    # #print(f'ii: {ii}')
    # #print(f'ii.shape: {ii.shape}')
    #print('f_b[ii]', f_b[ii])
    # if len(ii) > 0:
    if valores_positivos.size > 0:
        # #print('f_b[ii]:',f_b[ii].shape)
        # min_tot = np.nanmin(f_b[ii])
        # #print('min_tot:',min_tot)
        # min_tot = np.nanmin(valores_positivos)
        min_tot = np.nanmin(valores_positivos)
        #print("Mínimo valor positivo:", min_tot)
        for k in range(k0,k1+1):
            tot = np.nansum(Vr2pidVr_b*(np.matmul(f_b[:,:,k],dVx_b)))
            if tot > min_tot:
                f_b[:,:,k] = f_b[:,:,k] * tot_b[k]/tot

    #print('&'*88)
    #print('&'*88)
    #print('&'*88)
    #print('&'*88)
    #print(f'f_b[:, 0, 0]:{f_b[:, 0, 0]}')
    #print(f'f_b[0, :, 0]:{f_b[0, :, 0]}')
    #print(f'f_b[0, 0, :]:{f_b[0, 0, :]}')
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))  # 1 fila, 3 columnas
    # Primer scatter
    axs[0].scatter(vr_b, f_b[:, 0, 0])
    axs[0].plot(vr_b, f_b[:, 0, 0])
    axs[0].set_title('001vr_b vs f_b[:,0,0]')
    axs[0].grid(True)

    # Segundo scatter
    axs[1].scatter(vx_b, f_b[0, :, 0])
    axs[1].plot(vx_b, f_b[0, :, 0])
    axs[1].set_title('001vx_b vs f_b[0,:,0]')
    axs[1].grid(True)

    # Tercer scatter
    axs[2].scatter(x_b, f_b[0, 0, :])
    axs[2].plot(x_b, f_b[0, 0, :])
    axs[2].set_title('001x_b vs f_b[0,0,:]')
    axs[2].grid(True)

    plt.tight_layout()
    plt.show()
    if debug:
        n_a  = np.zeros(nx_a)
        ux_a = np.zeros(nx_a)
        T_a  = np.zeros(nx_a)
        vr2_vx2_ran2 = np.zeros((nvr_a,nvx_a))
        for k in range(0,nx_a):
            n_a[k] = np.nansum(Vr2pidVr_a * (np.matmul(f_a[:, :, k],dVx_a)))
            if n_a[k] > 0:
                ux_a[k] = vth_a*np.nansum(Vr2pidVr_a*(np.matmul(f_a[:,:,k],(vx_a*dVx_a)))) / n_a[k]
                # for i in range(0,nvr_a):
                #     vr2_vx2_ran2[i,:] = (vr_a[i]**2) + ((vx_a-ux_a[k])/(vth_a))**2
                vr2_vx2_ran2 = (vr_a[:, None]**2) + ((vx_a - ux_a[k]) / vth_a)**2
                T_a[k]  = (mu * mH) * vth_a2 * (np.nansum(Vr2pidVr_a * (np.matmul((vr2_vx2_ran2 * f_a[:,:,k]), dVx_a)))) / (3 * q * n_a[k])

        n_b  = np.zeros(nx_b)
        ux_b = np.zeros(nx_b)
        T_b  = np.zeros(nx_b)
        vr2_vx2_ran2 = np.zeros((nvr_b, nvx_b))
        for k in range(nx_b):
            n_b[k] = np.nansum(Vr2pidVr_b * (np.matmul(f_b[:, :, k], dVx_b)))
            if n_b[k] > 0:
                ux_b[k] = vth_b * np.nansum(Vr2pidVr_b * np.matmul(f_b[:, :, k], (vx_b * dVx_b))) / n_b[k]
                vr2_vx2_ran2 = (vr_b[:, None]**2) + ((vx_b - ux_b[k]) / vth_b)**2
                T_b[k]  = (mu * mH) * vth_b2 * (np.nansum(Vr2pidVr_b * (np.matmul((vr2_vx2_ran2 * f_b[:,:,k]), dVx_b)))) / (3 * q * n_b[k])
    
    # # Crear una figura con 2 filas y 3 columnas
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))  # 2 filas, 3 columnas

    # Primer gráfico: vr_a vs f_a[:, 0, 0]
    axs[0, 0].scatter(vr_a, f_a[:, 0, 0], color='blue', label='Scatter')
    axs[0, 0].plot(vr_a, f_a[:, 0, 0], color='red', label='Line')
    axs[0, 0].set_title('vr_a vs f_a[:,0,0]')
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    # Segundo gráfico: vx_a vs f_a[0, :, 0]
    axs[0, 1].scatter(vx_a, f_a[0, :, 0], color='blue', label='Scatter')
    axs[0, 1].plot(vx_a, f_a[0, :, 0], color='red', label='Line')
    axs[0, 1].set_title('vx_a vs f_a[0,:,0]')
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # Tercer gráfico: x_a vs f_a[0, 0, :]
    axs[0, 2].scatter(x_a, f_a[0, 0, :], color='blue', label='Scatter')
    axs[0, 2].plot(x_a, f_a[0, 0, :], color='red', label='Line')
    axs[0, 2].set_title('x_a vs f_a[0,0,:]')
    axs[0, 2].grid(True)
    axs[0, 2].legend()

    # Cuarto gráfico: vr_b vs f_b[:, 0, 0]
    axs[1, 0].scatter(vr_b, f_b[:, 0, 0], color='blue', label='Scatter')
    axs[1, 0].plot(vr_b, f_b[:, 0, 0], color='red', label='Line')
    axs[1, 0].set_title('vr_b vs f_b[:,0,0]')
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    # Quinto gráfico: vx_b vs f_b[0, :, 0]
    axs[1, 1].scatter(vx_b, f_b[0, :, 0], color='blue', label='Scatter')
    axs[1, 1].plot(vx_b, f_b[0, :, 0], color='red', label='Line')
    axs[1, 1].set_title('vx_b vs f_b[0,:,0]')
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    # Sexto gráfico: x_b vs f_b[0, 0, :]
    axs[1, 2].scatter(x_b, f_b[0, 0, :], color='blue', label='Scatter')
    axs[1, 2].plot(x_b, f_b[0, 0, :], color='red', label='Line')
    axs[1, 2].set_title('x_b vs f_b[0,0,:]')
    axs[1, 2].grid(True)
    axs[1, 2].legend()

    # Ajustar el diseño para evitar superposiciones
    plt.tight_layout()
    plt.show()


    if w_new:
        if w1_active:
            obj_interp_fvrvxx_2.vra2    = vr_a
            obj_interp_fvrvxx_2.vxa2    = vx_a
            obj_interp_fvrvxx_2.Tnorma2 = Tnorm_a  
            obj_interp_fvrvxx_2.vrb2    = vr_b     
            obj_interp_fvrvxx_2.vxb2    = vx_b     
            obj_interp_fvrvxx_2.Tnormb2 = Tnorm_b 
            obj_interp_fvrvxx_2.weight2 = weight
            print('Storing Weight in Weight2')
        else:
            obj_interp_fvrvxx_1.vra1    = vr_a     
            obj_interp_fvrvxx_1.vxa1    = vx_a     
            obj_interp_fvrvxx_1.Tnorma1 = Tnorm_a  
            obj_interp_fvrvxx_1.vrb1    = vr_b     
            obj_interp_fvrvxx_1.vxb1    = vx_b     
            obj_interp_fvrvxx_1.Tnormb1 = Tnorm_b  
            obj_interp_fvrvxx_1.weight1 = weight
            print('Storing Weight in Weight1')

    #print('$'*150)
    return f_b

if __name__ == '__main__':
    # # Datos de entrada
    # datos = np.load('datos_interp_fvrvxx.npz')
    datos = np.load('datos_interp_fvrvxx_20250423_204806.npz')
    

    f_b_tested  = datos['fHM']
    f_a  = datos['fH'].T
    vr_a = datos['vrA']
    vx_a = datos['vxA']
    x_a  = datos['xH']
    Tnorm_a = datos['TnormA']

    vr_b = datos['vrM']
    vx_b = datos['vxM']
    x_b  = datos['xH2']
    Tnorm_b = datos['TnormM']
    warn = datos['do_warn']

    

    #print(obj_interp_fvrvxx_1.vra1   ) 
    #print(obj_interp_fvrvxx_1.vxa1   ) 
    #print(obj_interp_fvrvxx_1.Tnorma1) 
    #print(obj_interp_fvrvxx_1.vrb1   ) 
    #print(obj_interp_fvrvxx_1.vxb1   ) 
    #print(obj_interp_fvrvxx_1.Tnormb1)

    f_b = interp_f_vr_vx_X(
            f_a     =f_a,
            vr_a    =vr_a,
            vx_a    =vx_a,
            x_a     =x_a,
            Tnorm_a =Tnorm_a,
            vr_b    =vr_b,
            vx_b    =vx_b,
            x_b     =x_b,
            Tnorm_b =Tnorm_b,
            debug   =1,
            correct =1,
            warn    =1e-3
        )

    #print('$'*5)
    #print(f_b.shape)
    #print(f_a.shape)

    # plt.plot(f_b[:,0,0])
    # plt.plot(f_b[:,0,0],'o')
    # plt.show()
    #print(obj_interp_fvrvxx_1.vra1   ) 
    #print(obj_interp_fvrvxx_1.vxa1   ) 
    #print(obj_interp_fvrvxx_1.Tnorma1) 
    #print(obj_interp_fvrvxx_1.vrb1   ) 
    #print(obj_interp_fvrvxx_1.vxb1   ) 
    #print(obj_interp_fvrvxx_1.Tnormb1)
    # f_b = interp_f_vr_vx_X(
    #         f_a     =f_a,
    #         vr_a    =vr_a,
    #         vx_a    =vx_a,
    #         x_a     =x_a,
    #         Tnorm_a =Tnorm_a,
    #         vr_b    =vr_b,
    #         vx_b    =vx_b,
    #         x_b     =x_b,
    #         Tnorm_b =Tnorm_b,
    #         debug   =1,
    #         correct =0,
    #         warn    =1e-3
    #     )
