#from cmath import sqrt     replaced with np.sqrt - nh
from reverse import * # fixed import
import numpy as np

# sets up optimum Vr and Vx velocity space mesh for Kinetic_Neutrals procedure 
# Input: 
#   nv - Integer, number of elements desired in vr mesh
#   Ti - arrray, Ti profile
#   E0 - array, energy where a velocity is desired ( optional )
#   Tmax - float, ignore Ti above this value
#
# Gwendolyn Galleher 

# def create_VrVxMesh(nv, Ti, E0 = np.array([0.0]), Tmax = 0.0): # removed unused input parameters
#     print('######################################')
#     print('nv:',nv)
#     print('Ti:',Ti)
#     _Ti = np.array(Ti) 
#     _Ti = np.concatenate([_Ti, E0[E0>0]]) ## Check point <-- tracing (08/18/2024)
#     ####
#     ##
#     ##
#     if Tmax > 0:
#         _Ti = _Ti[_Ti<Tmax] # simplified previous lines by replacing loops - nh
#     maxTi = np.max(_Ti)
#     minTi = min(_Ti)
#     Tnorm = np.mean(_Ti)
#     vmax = 3.5
#     if maxTi - minTi <= 0.1 * maxTi: # changed < to <= like from IDL code - nh
#         v = np.arange(nv+1)*vmax/nv # added *vmax/nv from IDL code - nh // deleted float type because it caused bug - GG
#     else:
#         G = 2 * nv * np.sqrt(minTi / maxTi) / (1- np.sqrt(minTi / maxTi))
#         b = vmax / (nv * (nv + G))
#         a = G * b 
#         v = a * np.arange(nv +1) + b * (np.arange(nv+1) ** 2)

#         # Option: add velocity bins corresponding to E0 
        
#     v0=0 # this line was missing - nh
#     for k in range(np.size(E0)): # fixed range - nh
#         if E0[k] > 0.0:
#             v0 = np.sqrt(E0[k]/Tnorm)
#             ii = np.where(v > v0)[0] # ii = np.argwhere(v > v0).T[0] # argwhere has a weird output, but this should put it in a usable format - nh
#             if np.size(ii) > 0: # removed count variable
#                 v = np.concatenate([v[0 : ii[0]], [v0], v[ii[0] : ]])
#             else: 
#                 v = np.concatenate([v, v0]) # changed lines to np.concatenate - nh
        
#     vr = v[1 : ] # fixed indexing - nh
#     vx = np.concatenate([-reverse(vr), vr]) 

#     ixE0 = np.where(np.abs(vx) == v0)[0] # ixE0 = np.argwhere(abs(vx) == v0).T[0] # modified argwhere output - nh
#     if np.size(ixE0) == 1: # removed count variable; same changes made to following lines - nh
#         ixE0 = ixE0[0]
#     irE0 = np.where(vr == v0)[0] # irE0 = np.argwhere(vr == v0).T[0]
#     if np.size(irE0) == 1:
#         irE0 = irE0[0]
#     return vx,vr,Tnorm,ixE0,irE0 # changed to return outputs as a variables - GG



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
import copy
import numpy as np
from tqdm import tqdm
''' 
Version 01
Authors: Julio Balbin, Carlo Becerra
Date: August 20th, 2024
'''
# def create_vr_vx_mesh(nv: np.int32, Ti: np.ndarray, E0: np.ndarray = None, Tmax: np.float64 = 0.0):
#     if E0 is None:
#         E0 = np.array([0.0])
def create_vr_vx_mesh(nv:np.int32, Ti:np.array, E0=np.array([0.0]), Tmax: np.float64 = 0.0):

    ''' 
    nv & Ti are from:
    -> 1090904024_950to1050.sav

    Sets up optimum vr and vx velocity space mesh for Kinetic_Neutrals procedure 
    Input: 
    nv:     np.integer, number of elements desired in vr mesh
    Ti:     np.arrray,  Ti profile
    E0:     np.array,   energy where a velocity is desired ( optional )
    Tmax:   np.float64, ignore Ti above this value
    '''
    _Ti = copy.copy(Ti)  

    _Ti = np.concatenate((Ti, E0[E0 > 0]))  # Añadir valores de E0 a Ti
    if Tmax > 0:
        ii = np.where(_Ti < Tmax)
        _Ti = copy.copy(_Ti[ii])

    maxTi = np.nanmax(_Ti)
    minTi = np.nanmin(_Ti) 

    # Tnorm = np.nansum(_Ti)/len(_Ti) # np.count_nonzero(~np.isnan(data))
    Tnorm = np.nanmean(_Ti) # np.count_nonzero(~np.isnan(data))
    vmax  = 3.5

    if maxTi-minTi < 0.1*maxTi:
        v = np.arange(nv+1, dtype=np.float64)*(vmax/nv)
        print(v)
    else:
        G = 2*nv*np.sqrt(minTi/maxTi)/(1-np.sqrt(minTi/maxTi))
        b = vmax/(nv*(nv+G))
        a = G*b
        v = (a*np.arange(nv+1, dtype=np.float64)) + \
            (b*np.arange(nv+1, dtype=np.float64)**2)
    v0 = 0.0
    # Option: add velocity bins corresponding to E0
    for k in tqdm(range(0,len(E0)),desc=f'create_vrvxmesh.py: calc. v'):
        if E0[k]>0.0:
            v0 = np.sqrt(E0[k]/Tnorm)
            ii = np.where(v>v0)
            if len(ii[0]) > 0:
                v = np.concatenate((v[:ii[0][0]], [v0], v[ii[0][0]:]))
            else:
                v = np.concatenate((v, [v0]))

    vr = v[1:]

    vr_reversed = np.array(list(reversed(vr)))
    vx = np.concatenate((-vr_reversed, vr))

    ixE0 = np.where(np.abs(vx)==v0)
    if len(ixE0[0]==1):
        ixE0=ixE0[0]

    irE0 = np.where(vr==v0)
    print(irE0)
    if len(irE0[0]==1):
        irE0=irE0[0]

    return vx,vr,Tnorm,ixE0,irE0

if __name__ == "__main__":
    Ti = np.array([9.34678841, 9.34678841, 9.34678841, 9.34678841, 9.34678841,  
        9.34678841, 9.34678841, 9.34678841, 9.34678841, 9.34678841,  
        9.34678841, 9.34678841, 9.34678841, 9.34678841, 9.34678841,  
        9.34678841, 9.34678841, 9.34678841, 9.35320316, 9.37376491,  
        9.39664221, 9.55090158, 9.74462981, 10.01682774, 10.4194201,  
        10.85043168, 11.30267208, 11.76084605, 12.24974726, 12.78601165,  
        13.3353179, 13.89692424, 14.46181227, 15.00772385, 15.52499086,  
        16.03142479, 16.46155264, 16.86981647, 17.2232248, 17.51513957,  
        17.78057133, 17.92301564, 18.02953533, 18.08950149, 18.10810471,  
        18.10810471, 18.10810471, 18.10810471, 18.10810471, 18.10810471,  
        18.10917935, 18.12300874, 18.14222841, 18.20386174, 18.30656849,  
        18.44081541, 18.65152198, 18.90097004, 19.20510906, 19.55525068,  
        19.93960698, 20.36748855, 20.82293074, 21.29202885, 21.77093978,  
        22.25529679, 22.73730357, 23.21750178, 23.71196118, 24.21942508,  
        24.75868421, 25.38127615, 26.08562981, 26.98923976, 28.07260861,  
        29.32705803, 30.76285645, 32.39069096, 34.31278264, 36.50471401,  
        39.01950698, 41.95052863, 45.38358122, 48.89483031, 52.47814568,  
        56.13221024, 59.855224, 63.66947087, 67.30889627, 70.76738783,  
        74.12886042, 77.48437863, 80.83153765, 84.14353033, 87.4158932,  
        90.66862541, 93.91925609, 97.1668684, 100.40976548, 103.64686016,  
        106.88131942, 110.11531918, 113.34863026, 116.58109768, 119.81235607,  
        123.04318839, 126.27402108, 129.50485401, 132.73571541, 135.96659684,  
        139.19747161, 142.4283278, 145.65916564, 148.8900049, 152.12083297,  
        155.35165142, 158.58244833, 161.81320688, 165.04391078, 168.27451936,  
        171.50510443, 174.73568475, 177.9662626, 181.19682979, 184.42740785,  
        187.65798817, 190.88856844, 194.1191481, 197.34972896, 200.58029711,  
        203.81087269, 207.04146861, 210.27208269, 213.50268786, 216.73328038,  
        219.96386318, 223.19443731, 226.4250055, 229.65557301, 232.88615218,  
        236.11672325, 239.3472975, 242.57787344, 245.80844117, 249.03902203])
    nv = 20

    vx,vr,Tnorm,ixE0,irE0 = create_vr_vx_mesh(nv, Ti)
    print(vr)