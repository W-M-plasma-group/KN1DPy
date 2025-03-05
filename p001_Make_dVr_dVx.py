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
    """
    Constructs velocity space differentials for distribution functions.

    Parameters:
    -----------
    vr : np.ndarray
        Array of radial velocities.
    vx : np.ndarray
        Array of axial velocities.

    Returns:
    --------
    Vr2pidVr : np.ndarray
        Differential volume element for radial velocities.
    VrVr4pidVr : np.ndarray
        Differential volume element for radial velocities (higher order).
    dVx : np.ndarray
        Differential for axial velocities.
    vrL, vrR : np.ndarray
        Left and right boundaries for radial velocities.
    vxL, vxR : np.ndarray
        Left and right boundaries for axial velocities.
    vol : np.ndarray
        Volume elements in velocity space.
    vth_Deltavx, vx_Deltavx, vr_Deltavr : np.ndarray
        Auxiliary quantities for kinetic equations.
    vr2vx2 : np.ndarray
        Squared magnitude of the velocity.
    jpa, jpb, jna, jnb : int
        Indices for positive and negative axial velocities.
    """

    # nvr & nvx are taking the shape of vr & vx respectively
    nvr = vr.size
    nvx = vx.size

    # Calculations for r-dimension
    _vr    = np.append(vr, 2 * vr[-1] - vr[-2])
    vr_mid = np.concatenate(([0.0], 0.5 * (_vr + np.roll(_vr, -1))))

    vrR = np.roll(vr_mid, -1)[0:nvr]
    vrL = copy.copy(vr_mid)[0:nvr]

    Vr2pidVr    =         np.pi * (vrR**2 - vrL**2)
    VrVr4pidVr  = (4/3) * np.pi * (vrR**3 - vrL**3)

    # Calculations for x-dimension
    _vx = np.concatenate(([2 * vx[0] - vx[1]], vx, [2 * vx[-1] - vx[-2]]))

    vxR = 0.5 * (np.roll(_vx, -1) + _vx)[1:nvx+1]
    vxL = 0.5 * (np.roll(_vx,  1) + _vx)[1:nvx+1]

    dVx = vxR - vxL

    # Calc. volume
    # --------------------------------------------------------------------
    vol = np.zeros((nvr, nvx), dtype=np.float64)
    # # for i in tqdm(range(nvr), desc=f"Make_dVr_dVx: calc. vol"): 
    # for i in range(nvr):     
    #     vol[i, :] = Vr2pidVr[i] * dVx
    vol = Vr2pidVr[:, np.newaxis] * dVx
    # --------------------------------------------------------------------
    Deltavx = vxR - vxL
    Deltavr = vrR - vrL
    # --------------------------------------------------------------------
    vth_Deltavx = np.zeros((nvr+2,nvx+2))
    vx_Deltavx  = np.zeros((nvr+2,nvx+2))
    vr_Deltavr  = np.zeros((nvr+2,nvx+2))
    # for j in tqdm(range(1,nvr+1),desc=f"Make_dVr_dVx: calc. vth & vx"):
    # for j in range(1,nvr+1):
        # vth_Deltavx[j,1:nvx+1] = 1.0/Deltavx     # vth_Deltavx(i,1:nvx)=1.0/Deltavx
        # vx_Deltavx[ j,1:nvx+1] =  vx/Deltavx     #  vx_Deltavx(i,1:nvx)= vx/Deltavx
    # for k in tqdm(range(1,nvx+1),desc=f"Make_dVr_dVx: calc. vr"):
    # for k in range(1,nvx+1):
    #     vr_Deltavr[ 1:nvr+1,k] =  vr/Deltavr     #  vr_Deltavr(1:nvr,j)=vr/Deltavr
    vth_Deltavx[1:nvr+1, 1:nvx+1] = 1.0 / Deltavx
    vx_Deltavx[ 1:nvr+1, 1:nvx+1] = vx / Deltavx 
    vr_Deltavr[ 1:nvr+1, 1:nvx+1] = vr[:, np.newaxis] / Deltavr[:, np.newaxis]
    # --------------------------------------------------------------------

    # Compute v^2
    vr2vx2=np.zeros((nvr,nvx), dtype=np.float64)
    # for l in range(0,nvr):
    #     vr2vx2[l,:] = vr[l]**2 + vx**2
    vr2vx2 = vr[:, np.newaxis]**2 + vx**2 

    # vx's positive index
    jp = copy.copy(np.where(vx>0)[0])   # This saves the positives index of vx
    jpa = int(jp[0])                    # This saves the first index of jp
    jpb = int(jp[len(jp)-1])            # This saves the last index of jp

    # vx's negative index   
    jn = copy.copy(np.where(vx<0)[0])   # This saves the negatives index of vx
    jna = int(jn[0])                    # This saves the first index of jp
    jnb = int(jn[len(jn)-1])            # This saves the last index of jp

    return Vr2pidVr,VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,\
           vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,\
           jpa,jpb,jna,jnb



if __name__ == "__main__":

    # vr = np.array([ 0.07025444373403009,    0.14139327496887916,    0.21341649370454718,  
    #                 0.2863240999410342,     0.3601160936783402,     0.43479247491646505,  
    #                 0.510353243655409,      0.5867983998951719,     0.6641279436357537,  
    #                 0.7423418748771545,     0.8214401936193743,     0.901422899862413,  
    #                 0.9822899936062708,     1.0640414748509475,     1.146677343596443,  
    #                 1.230197599842758,      1.3146022435898914,     1.3998912748378438,  
    #                 1.4860646935866155,     1.5731224998362061,     1.6610646935866153,  
    #                 1.749891274837844,      1.8396022435898913,     1.9301975998427576,  
    #                 2.0216773435964432,     2.1140414748509473,     2.207289993606271,  
    #                 2.301422899862413,      2.396440193619374,      2.4923418748771544,  
    #                 2.5891279436357535,     2.6867983998951717,     2.785353243655409,  
    #                 2.8847924749164653,     2.9851160936783403,     3.086324099941034,  
    #                 3.1884164937045467,     3.291393274968879,      3.39525444373403, 3.5])
    # vx = np.array([-3.5,                   -3.39525444373403,      -3.291393274968879,  -3.1884164937045467,  
    #                -3.086324099941034,     -2.9851160936783403,    -2.8847924749164653,  
    #                -2.785353243655409,     -2.6867983998951717,    -2.5891279436357535,  
    #                -2.4923418748771544,    -2.396440193619374,     -2.301422899862413,  
    #                -2.207289993606271,     -2.1140414748509473,    -2.0216773435964432,  
    #                 -1.9301975998427576,   -1.8396022435898913,    -1.749891274837844,  
    #                 -1.6610646935866153, -  1.5731224998362061,    -1.4860646935866155,  
    #                 -1.3998912748378438,   -1.3146022435898914,    -1.230197599842758,  
    #                 -1.146677343596443,    -1.0640414748509475,    -0.9822899936062708,  
    #                 -0.901422899862413,    -0.8214401936193743,    -0.7423418748771545,  
    #                 -0.6641279436357537,   -0.5867983998951719,    -0.510353243655409,  
    #                 -0.43479247491646505,  -0.3601160936783402,    -0.2863240999410342,  
    #                 -0.21341649370454718,  -0.14139327496887916,   -0.07025444373403009,  
    #                  0.07025444373403009,   0.14139327496887916,    0.21341649370454718,  
    #                  0.2863240999410342,    0.3601160936783402,     0.43479247491646505,  
    #                  0.510353243655409,     0.5867983998951719,     0.6641279436357537,  
    #                  0.7423418748771545,    0.8214401936193743,     0.901422899862413,  
    #                  0.9822899936062708,    1.0640414748509475,     1.146677343596443,  
    #                  1.230197599842758,     1.3146022435898914,     1.3998912748378438,  
    #                  1.4860646935866155,    1.5731224998362061,     1.6610646935866153,  
    #                  1.749891274837844,     1.8396022435898913,     1.9301975998427576,  
    #                  2.0216773435964432,    2.1140414748509473,     2.207289993606271,  
    #                  2.301422899862413,     2.396440193619374,      2.4923418748771544,  
    #                  2.5891279436357535,    2.6867983998951717,     2.785353243655409,  
    #                  2.8847924749164653,    2.9851160936783403,     3.086324099941034,  
    #                  3.1884164937045467,    3.291393274968879,      3.39525444373403,   3.5])
    
    # Vr2pidVr, VrVr4pidVr,dVx,vrL,vrR,vxL,vxR,\
    # vol,vth_Deltavx,vx_Deltavx,vr_Deltavr,vr2vx2,\
    # jpa,jpb,jna,jnb = Make_dVr_dVx(vr,vx)

    # # print(Vr2pidVr)
    # print(jpa,jna)
    # print(jpb,jnb)
    # print("###################")
    # print('Test successfull!!!')
    # print("###################")
    # print(vrL)
    # print(vrR)
    # print('Test successfull!!!')

    vr = np.array([ 0.07025444373403009,    0.14139327496887916,    0.21341649370454718,  
                    0.2863240999410342,     0.3601160936783402,     0.43479247491646505,  
                    0.510353243655409,      0.5867983998951719,     0.6641279436357537,  
                    0.7423418748771545,     0.8214401936193743,     0.901422899862413,  
                    0.9822899936062708,     1.0640414748509475,     1.146677343596443,  
                    1.230197599842758,      1.3146022435898914,     1.3998912748378438,  
                    1.4860646935866155,     1.5731224998362061,     1.6610646935866153,  
                    1.749891274837844,      1.8396022435898913,     1.9301975998427576,  
                    2.0216773435964432,     2.1140414748509473,     2.207289993606271,  
                    2.301422899862413,      2.396440193619374,      2.4923418748771544,  
                    2.5891279436357535,     2.6867983998951717,     2.785353243655409,  
                    2.8847924749164653,     2.9851160936783403,     3.086324099941034,  
                    3.1884164937045467,     3.291393274968879,      3.39525444373403, 3.5])
    vx = np.array([-3.5,                   -3.39525444373403,      -3.291393274968879,  -3.1884164937045467,  
                   -3.086324099941034,     -2.9851160936783403,    -2.8847924749164653,  
                   -2.785353243655409,     -2.6867983998951717,    -2.5891279436357535,  
                   -2.4923418748771544,    -2.396440193619374,     -2.301422899862413,  
                   -2.207289993606271,     -2.1140414748509473,    -2.0216773435964432,  
                    -1.9301975998427576,   -1.8396022435898913,    -1.749891274837844,  
                    -1.6610646935866153, -  1.5731224998362061,    -1.4860646935866155,  
                    -1.3998912748378438,   -1.3146022435898914,    -1.230197599842758,  
                    -1.146677343596443,    -1.0640414748509475,    -0.9822899936062708,  
                    -0.901422899862413,    -0.8214401936193743,    -0.7423418748771545,  
                    -0.6641279436357537,   -0.5867983998951719,    -0.510353243655409,  
                    -0.43479247491646505,  -0.3601160936783402,    -0.2863240999410342,  
                    -0.21341649370454718,  -0.14139327496887916,   -0.07025444373403009,  
                     0.07025444373403009,   0.14139327496887916,    0.21341649370454718,  
                     0.2863240999410342,    0.3601160936783402,     0.43479247491646505,  
                     0.510353243655409,     0.5867983998951719,     0.6641279436357537,  
                     0.7423418748771545,    0.8214401936193743,     0.901422899862413,  
                     0.9822899936062708,    1.0640414748509475,     1.146677343596443,  
                     1.230197599842758,     1.3146022435898914,     1.3998912748378438,  
                     1.4860646935866155,    1.5731224998362061,     1.6610646935866153,  
                     1.749891274837844,     1.8396022435898913,     1.9301975998427576,  
                     2.0216773435964432,    2.1140414748509473,     2.207289993606271,  
                     2.301422899862413,     2.396440193619374,      2.4923418748771544,  
                     2.5891279436357535,    2.6867983998951717,     2.785353243655409,  
                     2.8847924749164653,    2.9851160936783403,     3.086324099941034,  
                     3.1884164937045467,    3.291393274968879,      3.39525444373403,   3.5])
    
    # Llamada a la función
    Vr2pidVr, VrVr4pidVr, dVx, vrL, vrR, vxL, vxR,\
    vol, vth_Deltavx, vx_Deltavx, vr_Deltavr, vr2vx2,\
    jpa, jpb, jna, jnb = Make_dVr_dVx(vr, vx)

    # Imprimir dimensiones y tipos de datos de cada salida
    outputs = {
        "Vr2pidVr": Vr2pidVr,
        "VrVr4pidVr": VrVr4pidVr,
        "dVx": dVx,
        "vrL": vrL,
        "vrR": vrR,
        "vxL": vxL,
        "vxR": vxR,
        "vol": vol,
        "vth_Deltavx": vth_Deltavx,
        "vx_Deltavx": vx_Deltavx,
        "vr_Deltavr": vr_Deltavr,
        "vr2vx2": vr2vx2,
        "jpa": jpa,
        "jpb": jpb,
        "jna": jna,
        "jnb": jnb
    }

    for name, value in outputs.items():
        if isinstance(value, np.ndarray):
            print(f"{name}: shape={value.shape}, dtype={value.dtype}")
        else:
            print(f"{name}: value={value}, type={type(value)}")

    print("###################")
    print('Test successfull!!!')
    print("###################")