from p000_variables     import VARIABLES
from p001_Make_dVr_dVx  import Make_dVr_dVx
from tqdm import tqdm
import numpy as np
# import sympy as sp
import copy
import time


def create_maxwellian_include(vr:       np.ndarray, 
                              vx:       np.ndarray, 
                              vx_shift: np.ndarray,
                              Tmaxwell: np.ndarray,
                              Tnorm:    np.float32,
                              mu:       np.float32,
                              mol:      np.int32
                              ):
    q   = VARIABLES.q
    mH  = VARIABLES.mH

    Vr2pidVr, VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,jpa,jpb,jna,jnb = Make_dVr_dVx(vr,vx)

    nvr = len(vr)
    nvx = len(vx)
    nx  = len(vx_shift)
    vth = np.sqrt(2*q*Tnorm/(mu*mH))
    vth2= vth *vth
    vth3= vth2*vth
    vr2vx2_ran2=np.zeros((nvr,nvx), dtype = np.float64)
    
    
    AN  = np.zeros((nvr,nvx,2), dtype=np.float64)
    BN  = np.zeros((nvr,nvx,2), dtype=np.float64)
    sgn = [1,-1]
    maxwell = np.zeros((nvr,nvx,nx),dtype=np.float32)
    #################################################################
    for k in tqdm(range(nx),desc=f'k'):
        if Tmaxwell[k]>0.0:
            # # for i in tqdm(range(nvr),desc=f'i'):
            # for i in range(nvr):
            #     arg = -((vr[i]**2 + ((vx) - (vx_shift[k] / vth))**2) * mol * Tnorm / Tmaxwell[k])  # Error de \pm 3e-3
            #     arg = np.where(np.logical_and((-80 < arg),(arg < 0.0)), arg, -80)
            #     maxwell[i,:,k] = np.exp(arg)
            arg              = -((vr[:, np.newaxis]**2 + (vx - (vx_shift[k] / vth))**2) * mol * Tnorm / Tmaxwell[k])
            # print(f'vr[:,np.newaxis]:{vr[:,np.newaxis].shape}')
            # print(f'vr:{vr.shape}')
            # print(f'(vx - (vx_shift[k] / vth)):{(vx - (vx_shift[k] / vth)).shape}')
            # print(f'arg:{arg.shape}')
            arg              = np.where(np.logical_and((-80 < arg), (arg < 0.0)), arg, -80)
            maxwell[:, :, k] = np.exp(arg)
            # print(maxwell[:,:,k].shape)

            # variable = np.array(sp.Matrix(maxwell[:, :, k]) * sp.Matrix(dVx), dtype=float).flatten()
            # maxwell[:,:,k] = copy.copy(maxwell[:,:,k]/(np.nansum(Vr2pidVr*variable)))

            # This change reduced 1.5s/iteration
            # variable = np.sum(maxwell[:, :, k] * dVx, axis=1)
            variable = np.matmul(maxwell[:, :, k], dVx)
            maxwell[:, :, k] = maxwell[:, :, k] / np.nansum(Vr2pidVr * variable)
            
            # Compute desired moments
            WxD = copy.copy(vx_shift[k])
            ED = (WxD**2)+3*q*(Tmaxwell[k])/(mol*mu*mH)
            # Compute present moments of Maxwell, WxMax, and EMax 
            # these two next variables don't interact with the code
            # # vx_dvx = vx*dVx
            # # maxwell_2d  = copy.copy(maxwell[:, :, k])
            
            # This change also reduced the time in other 1.5s/ it
            # max_2d_xd   = np.array(sp.Matrix(maxwell[:, :, k]) * sp.Matrix(vx * dVx), dtype=float).flatten()
            # max_2d_xd = np.sum(maxwell[:, :, k] * vx * dVx, axis=1)
            max_2d_xd = np.matmul(maxwell[:, :, k],(vx * dVx))            
            WxMax       = vth  * (np.nansum(Vr2pidVr * (max_2d_xd) ))
            EMax        = vth2 * (np.nansum(Vr2pidVr*(np.matmul((vr2vx2*maxwell[:,:,k]),dVx))))
            # Compute Nij from Maxwell, padded with zeros
            Nij = np.zeros((nvr+2,nvx+2), dtype=np.float64)  # The order of the error is maintained with float64, but not for float 32
            Nij[1:nvr+1, 1:nvx+1]   = maxwell[:,:,k]*vol
            auxNij                  = Nij*vx_Deltavx
            Nijp1_vx_Dvx            = np.roll(auxNij,  shift=-1, axis=1)
            Nij_vx_Dvx              = Nij*vx_Deltavx
            Nijm1_vx_Dvx            = np.roll(auxNij,  shift= 1, axis=1)
            auxNij2                 = Nij*vr_Deltavr
            Nip1j_vr_Dvr            = np.roll(auxNij2, shift=-1, axis=0)
            Nij_vr_Dvr              = Nij*vr_Deltavr
            Nim1j_vr_Dvr            = np.roll(auxNij2, shift= 1, axis=0)
            # Compute Ap, Am, Bp, and Bm (0=p 1=m)
            aux_AN = np.zeros((nvx+2,nvr+2),dtype=np.float64)
            aux_AN = Nij*vth_Deltavx

            _AN         = np.zeros((nvx+2,nvr+2), dtype=np.float64)
            _AN         = np.roll(aux_AN, shift=1, axis=1) - aux_AN
            AN[:,:,0]   = copy.copy(_AN[1:nvr+1,1:nvx+1])
            
            _AN         = np.zeros((nvx+2,nvr+2), dtype=np.float64)
            _AN         = -np.roll(aux_AN, shift=-1, axis=1) + aux_AN
            AN[:,:,1]   = copy.copy(_AN[1:nvr+1,1:nvx+1])

            BN[:,jpa+1:jpb+1,0] =  Nijm1_vx_Dvx[1:nvr+1,jpa+2:jpb+2] - Nij_vx_Dvx[1:nvr+1,jpa+2:jpb+2]
            BN[:,jpa,0]         = -Nij_vx_Dvx[1:nvr+1,jpa+1]
            BN[:,jnb,0]         =  Nij_vx_Dvx[1:nvr+1,jnb+1]
            BN[:,jna:jnb,0]     = -Nijp1_vx_Dvx[1:nvr+1,jna+1:jnb+1] + Nij_vx_Dvx[1:nvr+1,jna+1:jnb+1]
            BN[:,:,0]           =  BN[:,:,0] + Nim1j_vr_Dvr[1:nvr+1,1:nvx+1] - Nij_vr_Dvr[1:nvr+1,1:nvx+1]

            BN[:,jpa+1:jpb+1,1] = -Nijp1_vx_Dvx[1:nvr+1,jpa+2:jpb+2] + Nij_vx_Dvx[1:nvr+1,jpa+2:jpb+2]
            BN[:,jpa,1]         = -Nijp1_vx_Dvx[1:nvr+1,jpa+1]
            BN[:,jnb,1]         =  Nijm1_vx_Dvx[1:nvr+1,jnb+1]
            BN[:,jna:jnb,1]     =  Nijm1_vx_Dvx[1:nvr+1,jna+1:jnb+1] - Nij_vx_Dvx[1:nvr+1,jna+1:jnb+1]
            BN[1:nvr,:,1]       =  BN[1:nvr,:,1] - Nip1j_vr_Dvr[2:nvr+1,1:nvx+1] + Nij_vr_Dvr[2:nvr+1,1:nvx+1]
            BN[0,:,1]           =  BN[0,:,1] - Nip1j_vr_Dvr[1,1:nvx+1]

            # Remove padded zeros in Nij
            Nij = Nij[1:nvr+1,1:nvx+1]

            # Cycle through 4 possibilies of sign(a_Max),sign(b_Max)
            TB1 = np.zeros(2)#, dtype=np.float32)
            TB2 = np.zeros(2)#, dtype=np.float32)
            ia=0
            while ia<2:
                # Compute TA1, TA2    
                # aux_TA1 = np.array(sp.Matrix(AN[:,:,ia]*sp.Matrix(vx)),dtype=float).flatten()
                # TA1 = vth*np.sum(aux_TA1)
                aux_TA1 = np.sum(np.matmul(AN[:, :, ia], vx))
                TA1 = vth*aux_TA1
                TA2 = vth2*np.sum(vr2vx2*AN[:,:,ia])
                ib  = 0
                while ib<2:
                    # Compute TB1, TB2
                    if TB1[ib]==0:
                        # aux_TB1 = np.array(sp.Matrix(BN[:,:,ib])*sp.Matrix(vx),dtype=float).flatten()
                        # TB1[ib] = vth*np.sum(aux_TB1)
                        aux_TB1 = np.sum(np.matmul(BN[:, :, ib], vx))
                        TB1[ib] = vth*aux_TB1

                        # print('TB1:',TB1, 'ib:',ib)
                    if TB2[ib]==0:
                        TB2[ib] = vth2*np.sum(vr2vx2*BN[:,:,ib])
                    denom = TA2*TB1[ib]-TA1*TB2[ib]
                    # print('denom:',denom, 'ib:',ib)
                    b_Max = 0.0
                    a_Max = 0.0
                    if denom!=0 and TA1 !=0:
                        b_Max = (TA2*(WxD-WxMax) - TA1*(ED-EMax))/denom
                        a_Max = (WxD-WxMax-TB1[ib]*b_Max)/TA1
                    if a_Max*sgn[ia] > 0 and b_Max*sgn[ib]>0:
                        maxwell[:,:,k] = (Nij + AN[:,:,ia]*a_Max + BN[:,:,ib]*b_Max)/vol
                        # np.savetxt('salida_23.dat', maxwell[:,:,k], fmt='%.4e',delimiter='\t')
                        ia = 2
                        ib = 2
                    ib = ib + 1
                ia= ia + 1

            # aux_maxwell = np.array(sp.Matrix(maxwell[:,:,k])*sp.Matrix(dVx),dtype=float).flatten()
            # maxwell[:,:,k] = maxwell[:,:,k]/np.sum(Vr2pidVr*aux_maxwell)

            aux_maxwell = np.sum(Vr2pidVr * (np.matmul(maxwell[:, :, k], dVx)))
            maxwell[:,:,k] = maxwell[:,:,k]/aux_maxwell
            
    return maxwell
                    
if __name__ == "__main__":
    import time

    from d003_Maxwellian_include import dataaset_d003_1
    from d003_Maxwellian_include import dataaset_d003_2
    from d003_Maxwellian_include import dataaset_d003_3
    from d003_Maxwellian_include import dataaset_d003_4
    from d003_Maxwellian_include import dataaset_d003_5
    from d003_Maxwellian_include import dataaset_d003_6
    from d003_Maxwellian_include import dataaset_d003_7


    def select_class(dataset_selected:object):
        vr              = dataset_selected.vr
        vx              = dataset_selected.vx
        vx_shift        = dataset_selected.vx_shift
        Tmaxwell        = dataset_selected.Tmaxwell
        Tnorm           = dataset_selected.Tnorm
        mu              = dataset_selected.mu
        mol             = dataset_selected.mol
        name_save_fig   = dataset_selected.name_save_fig
        return vr, vx, vx_shift,Tmaxwell,Tnorm,mu,mol, name_save_fig
    
    # Here you can select the dataset you want, from dataset_1 to dataset_7
    vr, vx, vx_shift,Tmaxwell,Tnorm,mu,mol,name_save_fig = select_class(dataaset_d003_2)
    start_time = time.time()
    maxwell = create_maxwellian_include(vr, vx, vx_shift,Tmaxwell,Tnorm,mu,mol)
    end_time = time.time()
    total_time = end_time - start_time
    
    print(total_time)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    from matplotlib.colors import LogNorm
    from matplotlib.animation import FuncAnimation, PillowWriter

    # output_gif_path = f"output_2D_{name_save_fig}.gif"

    # # Crear figura para la animación
    # fig, ax = plt.subplots(figsize=(6, 6))

    # # Verificar si hay valores negativos o ceros
    # if np.any(maxwell <= 0):
    #     maxwell[maxwell <= 0] = 1e-10  # Reemplazar ceros y valores negativos

    # # Configurar la imagen inicial con LogNorm
    # valor_min = np.min(maxwell[maxwell > 0])  # El valor mínimo positivo
    # valor_max = np.max(maxwell)
    # im = ax.imshow(maxwell[:, :, 0], cmap="coolwarm", norm=LogNorm(vmin=valor_min, vmax=valor_max))
    # ax.set_title("Layer 0")
    # plt.colorbar(im, ax=ax, label="Value")

    # # Función de actualización para cada frame
    # def update(frame):
    #     im.set_data(maxwell[:, :, frame])
    #     ax.set_title(f"Layer {frame}")

    # # Crear la animación
    # anim = FuncAnimation(fig, update, frames=maxwell.shape[2], interval=100)

    # # Guardar el GIF
    # anim.save(output_gif_path, writer=PillowWriter(fps=10))
    # plt.show()

    #####################################################################################
    #####################################################################################
    #####################################################################################
    #####################################################################################
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors    import LogNorm
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.animation import FuncAnimation, PillowWriter

    # Parámetros para el GIF
    output_gif_path = f"output_3D_{name_save_fig}.gif"

    # Crear figura para la animación 3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Verificar si hay valores negativos o ceros
    if np.any(maxwell <= 0):
        maxwell[maxwell <= 0] = 1e-10  # Reemplazar ceros y valores negativos

    # Crear malla de velocidades (vr, vx)
    vr = np.linspace(0, 1, maxwell.shape[0])  # Ajusta según tus datos reales
    vx = np.linspace(-1, 1, maxwell.shape[1])  # Ajusta según tus datos reales
    Vr, Vx = np.meshgrid(vr, vx, indexing='ij')

    # Configurar límites de color
    valor_min = np.min(maxwell[maxwell > 0])  # Valor mínimo positivo
    valor_max = np.max(maxwell)

    # Función de actualización para cada frame
    def update(frame):
        ax.clear()  # Limpiar el gráfico anterior
        surf = ax.plot_surface(
            Vr, Vx, maxwell[:, :, frame],
            cmap="coolwarm",
            norm=LogNorm(vmin=valor_min, vmax=valor_max),
            rstride=1, cstride=1, linewidth=0, antialiased=False
        )
        ax.set_title(f"Layer {frame}")
        ax.set_xlabel("Radial Velocity ($v_r$)")
        ax.set_ylabel("Axial Velocity ($v_x$)")
        ax.set_zlabel("Maxwellian Distribution")
        ax.set_zlim(valor_min, valor_max)  # Fijar límites en el eje Z
        return surf,

    # Crear la animación
    anim = FuncAnimation(fig, update, frames=maxwell.shape[2], interval=100)

    # Guardar el GIF
    anim.save(output_gif_path, writer=PillowWriter(fps=10))

    # Mostrar la animación
    plt.show()

    ##################################################
    ##################################################
    ##################################################
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import LogNorm
    # import numpy as np
    # from matplotlib.animation import FuncAnimation, PillowWriter

    # # Parámetros para el GIF
    # output_gif_path = f"output_2D_{name_save_fig}_2.gif"

    # # Crear figura para la animación
    # fig, ax = plt.subplots(figsize=(6, 6))

    # # Verificar si hay valores negativos o ceros
    # if np.any(maxwell <= 0):
    #     maxwell[maxwell <= 0] = 1e-10  # Reemplazar ceros y valores negativos

    # # Configurar la imagen inicial con LogNorm
    # valor_min = np.min(maxwell[maxwell > 0])  # El valor mínimo positivo
    # valor_max = np.max(maxwell)

    # # Crear malla de velocidades (vr, vx)
    # vr = np.linspace(0, 1, maxwell.shape[0])  # Ajusta según tus datos reales
    # vx = np.linspace(-1, 1, maxwell.shape[1])  # Ajusta según tus datos reales

    # # Número de puntos en el dominio espacial
    # nx = len(vx_shift)  # Esto define el número de frames en la animación

    # # Configurar la imagen inicial con LogNorm
    # im = ax.imshow(
    #     maxwell[:, :, 0],
    #     cmap="coolwarm",
    #     norm=LogNorm(vmin=valor_min, vmax=valor_max),
    #     extent=[vx[0], vx[-1], vr[0], vr[-1]],  # Extensión para etiquetar ejes
    #     origin="lower",  # Origen en la esquina inferior izquierda
    #     aspect="auto"    # Ajustar la relación de aspecto
    # )
    # ax.set_title("Layer 0")
    # ax.set_xlabel("Axial Velocity ($v_x$)")  # Etiqueta del eje x
    # ax.set_ylabel("Radial Velocity ($v_r$)")  # Etiqueta del eje y

    # # Añadir colorbar con etiqueta
    # cbar = plt.colorbar(im, ax=ax, label="Maxwellian Distribution Value (Log Scale)")

    # # Función de actualización para cada frame
    # def update(frame):
    #     im.set_data(maxwell[:, :, frame])
    #     ax.set_title(f"Layer {frame}")
    #     return im,

    # # Crear la animación
    # anim = FuncAnimation(fig, update, frames=nx, interval=100)  # Usar nx como número de frames

    # # Guardar el GIF
    # anim.save(output_gif_path, writer=PillowWriter(fps=10))

    # # Mostrar la animación
    # plt.show()