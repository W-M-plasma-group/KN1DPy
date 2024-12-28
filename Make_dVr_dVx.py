# # Make_dVr_dVx
# #   Constructs velocity space differentials for distrobution functions 
# # used by Kinetic_Neutrals, Kinetic_H2, Kinetic_H, and other related
# # procedures 
# #
# # Gwendolyn Galleher 

# import numpy as np

# def Make_dVr_dVx(Vr, Vx): # For this to work inputs must be np arrays so I might need to ammend this later to make sure all inputs are arrays
#     # Vr = np.asarray(Vr)
#     # Vx = np.asarray(Vx)
    
#     # # Convertir a listas
#     # Vr_list = Vr.tolist()
#     # Vx_list = Vx.tolist()
    
#     # # Crear las cadenas en el formato deseado
#     # Vr_str = f"[{', '.join(map(str, Vr_list))}]"
#     # Vx_str = f"[{', '.join(map(str, Vx_list))}]"
    
#     # # Abrir el archivo en modo de escritura
#     # with open('variables_Make_dvr_dvx.txt', 'w') as file:
#     #     file.write(f'vr= {Vr_str}, tipo: {type(Vr_list)}\n')
#     #     file.write(f'vx= {Vx_str}, tipo: {type(Vx_list)}\n')
    
#     # print('--------------->Vr= ', Vr)
#     # print('--------------->Vx= ', Vx)
#     # Vr=np.array(Vr)
#     # Vx=np.array(Vx) # sets inputs to np.array if they aren't already; resolves earlier comment - nh
    
#     # Determine velocity space differentials 
#     nVr = np.size(Vr)
#     nVx = np.size(Vx)

#     # for Vr first 
#     _Vr = np.concatenate([Vr, [2 * Vr[nVr-1] - Vr[nVr-2]]])  # this is an array and Vr(nVr-1) is calling the last cell of Vr
#     Vr_mid = np.concatenate([[0.0], 0.5 * (_Vr + np.roll(_Vr, -1))]) # changed to np.concatenate - nh
    
#     VrR = np.roll(Vr_mid, -1)
#     VrL = Vr_mid

#     Vr2pidVr = np.pi * ((VrR ** 2) - (VrL ** 2))
#     Vr2pidVr = Vr2pidVr[0 : nVr]  # makes it the same length as Vr 

#     VrVr4pidVr = (4/3) * np.pi * ((VrR ** 3) - (VrL ** 3))
#     VrVr4pidVr = VrVr4pidVr[0 : nVr]
#     VrR = VrR[0 : nVr]
#     VrL = VrL[0 : nVr] 

#     # now for Vx
#     _Vx = np.concatenate([[2 * Vx[0] - Vx[1]], Vx, [2 * Vx[nVx - 1] - Vx[nVx - 2]]]) # changed to np.concatenate - nh
#     VxR = 0.5 * (np.roll(_Vx, -1) + _Vx)
#     VxL = 0.5 * (np.roll(_Vx, 1) + _Vx)
#     dVx = VxR[1: nVx+1] - VxL[1:nVx+1]
#     VxR = VxR[1: nVx+1]
#     VxL = VxL[1 : nVx+1]

#     # compute volume elements 
#     vol = np.zeros((nVx, nVr), float)
#     for i in range(0, nVr):
#         vol[:,i] = Vr2pidVr[i] * dVx # fixed minor indexing bugs - nh
#     print('vol: ',vol.shape, type(vol))
#     #compute DeltaVx, DeltaVr
#     DeltaVx = VxR - VxL
#     DeltaVr = VrR - VrL

#     # compute vth_Deltavx, vx_Deltavx, vr_Deltavr, padded with zeros
#     Vth_DeltaVx = np.zeros((nVx + 2, nVr + 2), float)
#     Vx_DeltaVx = np.zeros((nVx + 2, nVr + 2), float)
#     Vr_DeltaVr = np.zeros((nVx + 2, nVr + 2), float) # replaced np.array with np.zeros - nh
#     for i in range(1, nVr+1):
#         Vth_DeltaVx[1 : nVx+1,i] = 1.0/DeltaVx
#         Vx_DeltaVx[1 : nVx+1,i] = Vx/DeltaVx
#     for j in range(1, nVx+1):
#         Vr_DeltaVr[j,1 : nVr+1] = Vr/DeltaVr 
    
#     #compute v^2
#     Vr2Vx2 = np.zeros((nVx, nVr), float)
#     for i in range(0, nVr):
#         Vr2Vx2[:,i] = (Vr[i] ** 2) + (Vx ** 2)

#     # Determine indice range of positive and negative Vx 
#     jpa=jpb=jna=jnb=-1
#     jp = np.argwhere(Vx > 0)
#     if jp.size>0:
#         jpa = jp[0][0]; jpb = jp[np.size(jp) - 1][0]
#     jn = np.argwhere(Vx < 0)
#     if jn.size>0:
#         jna = jn[0][0]; jnb = jn[np.size(jn) - 1][0] # modified section to return -1 if jp or jn is empty (previously this raised an error) - nh
    
#     # changed return line to provide an output as a list - nh
    
#     return Vr2pidVr,VrVr4pidVr,dVx,VrL,VrR,VxL,VxR,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2,jpa,jpb,jna,jnb
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

import numpy as np
import copy
from tqdm import tqdm


''' 
Version 01
Authors: Julio Balbin, Carlo Becerra
Date: August 17th, 2024
'''
def Make_dVr_dVx(vr: np.ndarray
                ,vx: np.ndarray):
    ''' 
    vr & vx are from:
    -> function: create_vr_vx_mesh

    Constructs velocity space differentials for distribution functions
    used by Kinetic_Neutrals.pro, Kinetic_H2.pro, Kinetic_H2.pro

    The abs(difference) between vol variable in python & IDL is in the range: 
    <(np.float64(1.433306555442826e-06), np.float64(6.70904781382986e-12)>
    '''

    # nvr & nvx are taking the shape of vr & vx respectively
    nvr = vr.size
    nvx = vx.size

    # Calculations for r-dimension
    _vr     = np.append(vr, 2 * vr[-1] - vr[-2])
    vr_mid  = np.concatenate(([0.0], 0.5 * (_vr + np.roll(_vr, -1))))

    vrR = np.roll(vr_mid, -1)[0:nvr]
    vrL = copy.copy(vr_mid)[0:nvr]

    Vr2pidVr    =         np.pi * (vrR**2 - vrL**2)
    VrVr4pidVr  = (4/3) * np.pi * (vrR**3 - vrL**3)

    # Calculations for x-dimension
    _vx = np.concatenate(([2 * vx[0] - vx[1]], vx, [2 * vx[-1] - vx[-2]]))

    vxR = 0.5 * (np.roll(_vx, -1) + _vx)[1:nvx+1]
    vxL = 0.5 * (np.roll(_vx,  1) + _vx)[1:nvx+1]

    dVx = vxR - vxL

    # Calc. volumen
    vol = np.zeros((nvr, nvx), dtype=np.float64)
    for i in tqdm(range(nvr), desc=f"Make_dVr_dVx: calc. vol"): 
        vol[i, :] = Vr2pidVr[i] * dVx
        
    Deltavx = vxR - vxL
    Deltavr = vrR - vrL

    vth_Deltavx =np.zeros((nvr+2,nvx+2))
    vx_Deltavx  =np.zeros((nvr+2,nvx+2))
    vr_Deltavr  =np.zeros((nvr+2,nvx+2))

    for j in tqdm(range(1,nvr+1),desc=f"Make_dVr_dVx: calc. vth & vx"):
        vth_Deltavx[j,1:nvx+1] = 1.0/Deltavx    # vth_Deltavx(i,1:nvx)=1.0/Deltavx
        vx_Deltavx[ j,1:nvx+1] =  vx/Deltavx     #  vx_Deltavx(i,1:nvx)= vx/Deltavx

    for k in tqdm(range(1,nvx+1),desc=f"Make_dVr_dVx: calc. vr"):
        vr_Deltavr[ 1:nvr+1,k] =  vr/Deltavr     #  vr_Deltavr(1:nvr,j)=vr/Deltavr
    
    # Compute v^2
    vr2vx2=np.zeros((nvr,nvx), dtype=np.float64)
    for l in range(0,nvr):
        vr2vx2[l,:] = vr[l]**2 + vx**2

    # vx's positive index
    jp = copy.copy(np.where(vx>0)[0])   # This saves the positives index of vx
    jpa = int(jp[0])                    # This saves the first index of jp
    jpb = int(jp[len(jp)-1])            # This saves the last index of jp

    # vx's negative index   
    jn = copy.copy(np.where(vx<0)[0])   # This saves the negatives index of vx
    jna = int(jn[0])                    # This saves the first index of jp
    jnb = int(jn[len(jn)-1])            # This saves the last index of jp

    return Vr2pidVr, VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,\
           vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,\
           jpa,jpb,jna,jnb



if __name__ == "__main__":

    vr = np.array([0.07025444373403009, 0.14139327496887916, 0.21341649370454718, 
                0.2863240999410342, 0.3601160936783402, 0.43479247491646505, 
                0.510353243655409, 0.5867983998951719, 0.6641279436357537, 
                0.7423418748771545, 0.8214401936193743, 0.901422899862413, 
                0.9822899936062708, 1.0640414748509475, 1.146677343596443, 
                1.230197599842758, 1.3146022435898914, 1.3998912748378438, 
                1.4860646935866155, 1.5731224998362061, 1.6610646935866153, 
                1.749891274837844, 1.8396022435898913, 1.9301975998427576, 
                2.0216773435964432, 2.1140414748509473, 2.207289993606271, 
                2.301422899862413, 2.396440193619374, 2.4923418748771544, 
                2.5891279436357535, 2.6867983998951717, 2.785353243655409, 
                2.8847924749164653, 2.9851160936783403, 3.086324099941034, 
                3.1884164937045467, 3.291393274968879, 3.39525444373403, 3.5])

    vx = np.array([-3.5, -3.39525444373403, -3.291393274968879, -3.1884164937045467, 
                -3.086324099941034, -2.9851160936783403, -2.8847924749164653, 
                -2.785353243655409, -2.6867983998951717, -2.5891279436357535, 
                -2.4923418748771544, -2.396440193619374, -2.301422899862413, 
                -2.207289993606271, -2.1140414748509473, -2.0216773435964432, 
                -1.9301975998427576, -1.8396022435898913, -1.749891274837844, 
                -1.6610646935866153, -1.5731224998362061, -1.4860646935866155, 
                -1.3998912748378438, -1.3146022435898914, -1.230197599842758, 
                -1.146677343596443, -1.0640414748509475, -0.9822899936062708, 
                -0.901422899862413, -0.8214401936193743, -0.7423418748771545, 
                -0.6641279436357537, -0.5867983998951719, -0.510353243655409, 
                -0.43479247491646505, -0.3601160936783402, -0.2863240999410342, 
                -0.21341649370454718, -0.14139327496887916, -0.07025444373403009, 
                0.07025444373403009, 0.14139327496887916, 0.21341649370454718, 
                0.2863240999410342, 0.3601160936783402, 0.43479247491646505, 
                0.510353243655409, 0.5867983998951719, 0.6641279436357537, 
                0.7423418748771545, 0.8214401936193743, 0.901422899862413, 
                0.9822899936062708, 1.0640414748509475, 1.146677343596443, 
                1.230197599842758, 1.3146022435898914, 1.3998912748378438, 
                1.4860646935866155, 1.5731224998362061, 1.6610646935866153, 
                1.749891274837844, 1.8396022435898913, 1.9301975998427576, 
                2.0216773435964432, 2.1140414748509473, 2.207289993606271, 
                2.301422899862413, 2.396440193619374, 2.4923418748771544, 
                2.5891279436357535, 2.6867983998951717, 2.785353243655409, 
                2.8847924749164653, 2.9851160936783403, 3.086324099941034, 
                3.1884164937045467, 3.291393274968879, 3.39525444373403, 3.5])

    Vr2pidVr, VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,\
    vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,\
    jpa,jpb,jna,jnb = Make_dVr_dVx(vr,vx)

    vol_db = np.loadtxt('draft_Make_dVr_dVx_vol.dat')
    vol_db = vol_db.reshape((40,80))
    diferencia = np.abs(vol_db -vol)
    
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    plt.imshow(diferencia, cmap='coolwarm', interpolation='nearest', norm=LogNorm())
    plt.colorbar(label='Difference')
    # plt.title("Mapa de diferencias entre vol1 y vol2")
    plt.savefig("draft_Make_dVr_dVx_vol_difference.pdf")
    plt.savefig("draft_Make_dVr_dVx_vol_difference.svg")
    plt.savefig("draft_Make_dVr_dVx_vol_difference.png")
    # plt.show()

    print('max: \t',diferencia.max())
    print('min: \t',diferencia.min())
    print("###################")
    print('Test successfull!!!')
    print("###################")