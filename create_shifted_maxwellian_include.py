import numpy as np
from sval import sval
from global_vars import mH, q
from scipy.ndimage import shift
# #   This INCLUDE file is used by Kinetic_H2 and Kinetic_H
# #   The code is also written within the create_shifted_maxwellian function

# # # def create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
# # #                                       nx,nvx,nvr,vth,vth2,maxwell,vr2vx2_ran2,
# # #                                       Vr2pidVr,dVx,Vol,Vth_DVx,Vx_DVx,Vr_DVr,vr2vx2_2D,jpa,jpb,jna,jnb): # added vr,vx,Tnorm,mu to arguments; should not have been left out originally
def create_shifted_maxwellian_include(vr: np.ndarray,       # (40,) - Array of radial velocity values
                                      vx: np.ndarray,       # (80,) - Array of x-velocity values
                                      Tnorm: float,         # Scalar value for temperature normalization
                                      vx_shift: np.ndarray, # (80,) - Array of x-velocity shifts
                                      Tmaxwell: np.ndarray, # (80,) - Array of Maxwellian temperature values
                                      Shifted_Maxwellian_Debug: int, # Integer flag for debugging
                                      mu: float,            # Scalar value for mu (2.0 in this case)
                                      mol: int,             # Integer value for molecule type
                                      nx: int,              # Integer value for the mesh size in x
                                      nvx: int,             # Integer value for the mesh size in vx
                                      nvr: int,             # Integer value for the mesh size in vr
                                      vth: float,           # Scalar value for thermal velocity
                                      vth2: float,          # Scalar value for squared thermal velocity
                                      maxwell: np.ndarray,  # (80, 80, 40) - 3D array of Maxwellian distribution function
                                      vr2vx2_ran2: np.ndarray, # (80, 40) - 2D array of vr2vx2 relation
                                      Vr2pidVr: np.ndarray, # (40,) - Array of Vr2pidVr values
                                      dVx: np.ndarray,      # (80,) - Array of differences in vx
                                      Vol: np.ndarray,      # (80, 40) - 2D array of volumes
                                      Vth_DVx: np.ndarray,  # (82, 42) - 2D array for Vth_DVx relation
                                      Vx_DVx: np.ndarray,   # (82, 42) - 2D array for Vx_DVx relation
                                      Vr_DVr: np.ndarray,   # (82, 42) - 2D array for Vr_DVr relation
                                      vr2vx2_2D: np.ndarray, # (80, 40) - 2D array of vr2vx2 in 2D
                                      jpa: int,             # Integer value for index jpa
                                      jpb: int,             # Integer value for index jpb
                                      jna: int,             # Integer value for index jna
                                      jnb: int) -> np.ndarray: # Integer value for index jnb; output is maxwell

# Input Parameters:
#   - Vx_shift: np.ndarray (nx) - Array of velocity shifts in x (m/s)
#   - Tmaxwell: np.ndarray (nx) - Array of Maxwellian temperatures (eV)
#   - Shifted_Maxwellian_Debug: int - If set, debugging information will be printed
#   - mol: int - Molecule type (1=atom, 2=diatomic molecule)

# Output:
#   - maxwell: np.ndarray (nvr, nvx, nx) - Shifted Maxwellian distribution function with
#     a numerically evaluated vx moment close to Vx_shift and a temperature close to Tmaxwell

# Notes on Algorithm:
#   One might think that maxwell could be computed simply by directly evaluating the EXP function:
#
#       for i=0, nvr-1 do begin
#           arg = -(vr(i)^2 + (vx - Vx_shift/vth)^2) * mol * Tnorm / Tmaxwell
#           maxwell(i, *, k) = exp(arg) > (-80)
#       endfor
#
#   However, due to the discrete velocity space bins, this method does not necessarily lead to a digital
#   representation of a shifted Maxwellian (maxwell) that, when numerically integrated, has the desired vx
#   moment of Vx_shift and temperature Tmaxwell.
#
#   To ensure that maxwell has the desired vx and T moments when evaluated numerically, a compensation
#   scheme is employed, similar to that used in Interp_fVrVxX.

  # variables = {
  #     'vr': vr, 'vx': vx, 'Tnorm': Tnorm, 'vx_shift': vx_shift, 'Tmaxwell': Tmaxwell, 
  #     'Shifted_Maxwellian_Debug': Shifted_Maxwellian_Debug, 'mu': mu, 'mol': mol, 
  #     'nx': nx, 'nvx': nvx, 'nvr': nvr, 'vth': vth, 'vth2': vth2, 'maxwell': maxwell, 
  #     'vr2vx2_ran2': vr2vx2_ran2, 'Vr2pidVr': Vr2pidVr, 'dVx': dVx, 'Vol': Vol, 
  #     'Vth_DVx': Vth_DVx, 'Vx_DVx': Vx_DVx, 'Vr_DVr': Vr_DVr, 'vr2vx2_2D': vr2vx2_2D, 
  #     'jpa': jpa, 'jpb': jpb, 'jna': jna, 'jnb': jnb
  # }
  # # Abrir el archivo en modo de escritura
  # with open('variables_created_shifted_maxwellian_include.txt', 'w') as file:
  #     for name, value in variables.items():
  #         if hasattr(value, '__len__'):
  #             file.write(f'{name}: {value}, tipo: {type(value)}, len({name}): {value.shape}\n')
  #         else:
  #             file.write(f'{name}: {value}, tipo: {type(value)}\n')

  print('we are inside create_sihfted_maxwellian_include')
  print('maxwell.shape= ',maxwell.shape,type(maxwell))
  print('Tmaxwell.sjape= ',Tmaxwell.shape,type(Tmaxwell))
  print('dvx.shape= ',dVx.shape,type(dVx))




  AN=BN=np.zeros((2,nvx,nvr)) # fixed array creation - nh
  sgn=[-1,1]
  for k in range(nx):
    if Tmaxwell[k]>0:
      for i in range(nvr):
        arg=-(vr[i]**2+(vx-vx_shift[k]/vth)**2) * mol*Tnorm/Tmaxwell[k]
        arg=np.minimum(arg,0)
        # maxwell[i,:,k]=np.exp(arg)>(-80)
        maxwell[k,:,i]=np.e**(np.maximum(arg,-80)) # replaced < and > with np.minimum and np.maximum - nh

      maxwell[k,:,:]/=np.sum(Vr2pidVr*(np.matmul(dVx,maxwell[k,:,:]))) # fixed matmul argument order - nh

      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out1=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+(vx-vx_out1/vth)**2
        T_out1=(mol*mu*mH)*vth2*np.sum(Vr2pidVr*(np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:])))/(3*q) # fixed matmul argument order - nh
        vth_local=0.1*np.sqrt(2*Tmaxwell[k]*q/(mol*mu*mH))
        Terror=abs(Tmaxwell[k]-T_out1)/Tmaxwell[k]
        Verror=abs(vx_out1-vx_shift[k])/vth_local

      # Compute desired moments

      WxD=vx_shift[k]
      ED=WxD**2+3*q*Tmaxwell[k]/(mol*mu*mH)

      # Compute present moments of maxwell, WxMax, and EMax 

      WxMax=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
      EMax=vth2*np.sum(Vr2pidVr*np.matmul(dVx,(vr2vx2_2D*maxwell[k,:,:]))) # fixed matmul argument order - nh

      # Compute Nij from maxwell, padded with zeros

      nij=np.zeros((nvr+2,nvx+2)).T
      nij[1:nvx+1,1:nvr+1]=maxwell[k,:,:]*Vol

      Nijp1_vx_Dvx=np.roll(nij*Vx_DVx,-1,0)
      Nij_vx_Dvx  =nij*Vx_DVx
      Nijm1_vx_Dvx=np.roll(nij*Vx_DVx,1,0)
      Nip1j_vr_Dvr=np.roll(nij*Vr_DVr,-1,1)
      Nij_vr_Dvr  =nij*Vr_DVr
      Nim1j_vr_Dvr=np.roll(nij*Vr_DVr,1,1)

      # Compute Ap, Am, Bp, and Bm (0=p 1=m)

      _AN=np.roll(nij*Vth_DVx,1,0)-nij*Vth_DVx
      AN[0,:,:]=_AN[1:nvx+1,1:nvr+1]
      _AN=-np.roll(nij*Vth_DVx,-1,0)+nij*Vth_DVx
      AN[0,:,:]=_AN[1:nvx+1,1:nvr+1]

      BN[0,jpa+1:jpb+1,:]=Nijm1_vx_Dvx[jpa+2:jpb+2,1:nvr+1]-Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[0,jpa,:]=-Nij_vx_Dvx[jpa+2,1:nvr+1]
      BN[0,jnb,:]=Nij_vx_Dvx[jnb+2,1:nvr+1]
      BN[0,jna:jnb,:]=-Nijp1_vx_Dvx[jna+1:jnb+1,1:nvr+1]+Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[0,:,:]=BN[0,:,:]+ Nim1j_vr_Dvr[1:nvx+1,1:nvr+1]-Nij_vr_Dvr[1:nvx+1,1:nvr+1]

      BN[1,jpa+1:jpb+1,:]=-Nijp1_vx_Dvx[jpa+2:jpb+2,1:nvr+1]+Nij_vx_Dvx[jpa+2:jpb+2,1:nvr+1]
      BN[1,jpa,:]=-Nijp1_vx_Dvx[jpa+1,1:nvr+1]
      BN[1,jnb,:]=Nijm1_vx_Dvx[jnb+1,1:nvr+1]
      BN[1,jna:jnb,:]=Nijm1_vx_Dvx[jna+1:jnb+1,1:nvr+1]-Nij_vx_Dvx[jna+1:jnb+1,1:nvr+1]
      BN[1,:,1:nvr]=BN[1,:,1:nvr] - Nip1j_vr_Dvr[1:nvx+1,2:nvr+1]+Nij_vr_Dvr[1:nvx+1,2:nvr+1]
      BN[1,:,0]=BN[1,:,0] - Nip1j_vr_Dvr[1:nvx+1,1]

      # Remove padded zeros in Nij

      nij=nij[1:nvx+1,1:nvr+1]

      # Cycle through 4 possibilies of sign(a_Max),sign(b_Max)

      TB1=TB2=np.zeros(2)
      ia=0
      while ia<2:

        # Compute TA1, TA2

        TA1=vth*np.sum(np.matmul(vx,AN[ia,:,:])) # fixed matmul argument order - nh
        TA2=vth2*np.sum(vr2vx2_2D*AN[ia,:,:])
        ib=0
        while ib<2:

          # Compute TB1, TB2

          if TB1[ib]==0:
            TB1[ib]=vth*np.sum(np.matmul(vx,BN[ib,:,:])) # fixed matmul argument order - nh
          if TB2[ib]==0:
            TB2[ib]=vth2*np.sum(vr2vx2_2D*BN[ib,:,:])
          denom=TA2*TB1[ib]-TA1*TB2[ib]
          b_max=a_max=0
          if denom!=0 and TA1!=0:
            b_Max=(TA2*(WxD-WxMax)-TA1*(ED-EMax))/denom
            a_Max=(WxD-WxMax-TB1[ib]*b_Max)/TA1
          if a_max*sgn[ia]>0 and b_max*sgn[ib]>0:
            maxwell[k,:,:]=(nij+AN[ia,:,:]*a_Max+BN[ib,:,:]*b_Max)/Vol
            ia=ib=2
          ib+=1
        ia+=1
      maxwell[k,:,:]=maxwell[k,:,:]/np.sum(Vr2pidVr*np.matmul(dVx,maxwell[k,:,:])) # fixed matmul argument order - nh

      if Shifted_Maxwellian_Debug: # fixed capitalization
        vx_out2=vth*np.sum(Vr2pidVr*np.matmul((vx*dVx),maxwell[k,:,:])) # fixed matmul argument order - nh
        for i in range(nvr):
          vr2vx2_ran2[:,i]=vr[i]**2+((vx-vx_out2)/vth)**2
        T_out2=(mol*mu*mH)*vth2*np.sum(Vr2pidVr*np.matmul(dVx,vr2vx2_ran2*maxwell[k,:,:]))/(3*q) # fixed matmul argument order - nh
        Terror2=abs(Tmaxwell[k]-T_out2)/Tmaxwell[k]
        Verror2=abs(vx_shift[k]-vx_out2)/vth_local
        print('CREATE_SHIFTED_MAXWELLIAN=> Terror:'+sval(Terror)+'->'+sval(Terror2)+'  Verror:'+sval(Verror)+'->'+sval(Verror2))

  return maxwell

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

def create_shifted_maxwellian_include(vr: np.ndarray,       # (40,) - Array of radial velocity values
                                      vx: np.ndarray,       # (80,) - Array of x-velocity values
                                      Tnorm: float,         # Scalar value for temperature normalization
                                      vx_shift: np.ndarray, # (80,) - Array of x-velocity shifts
                                      Tmaxwell: np.ndarray, # (80,) - Array of Maxwellian temperature values
                                      Shifted_Maxwellian_Debug: int, # Integer flag for debugging
                                      mu: float,            # Scalar value for mu (2.0 in this case)
                                      mol: int,             # Integer value for molecule type
                                      nx: int,              # Integer value for the mesh size in x
                                      nvx: int,             # Integer value for the mesh size in vx
                                      nvr: int,             # Integer value for the mesh size in vr
                                      vth: float,           # Scalar value for thermal velocity
                                      vth2: float,          # Scalar value for squared thermal velocity
                                      maxwell: np.ndarray,  # (80, 80, 40) - 3D array of Maxwellian distribution function
                                      vr2vx2_ran2: np.ndarray, # (80, 40) - 2D array of vr2vx2 relation
                                      Vr2pidVr: np.ndarray, # (40,) - Array of Vr2pidVr values
                                      dVx: np.ndarray,      # (80,) - Array of differences in vx
                                      Vol: np.ndarray,      # (80, 40) - 2D array of volumes
                                      Vth_DVx: np.ndarray,  # (82, 42) - 2D array for Vth_DVx relation
                                      Vx_DVx: np.ndarray,   # (82, 42) - 2D array for Vx_DVx relation
                                      Vr_DVr: np.ndarray,   # (82, 42) - 2D array for Vr_DVr relation
                                      vr2vx2_2D: np.ndarray, # (80, 40) - 2D array of vr2vx2 in 2D
                                      jpa: int,             # Integer value for index jpa
                                      jpb: int,             # Integer value for index jpb
                                      jna: int,             # Integer value for index jna
                                      jnb: int) -> np.ndarray: # Integer value for index jnb; output is maxwell
    # Inicializar arrays
    AN = np.zeros((nvr, nvx, 2), dtype=np.float64)
    BN = np.zeros((nvr, nvx, 2), dtype=np.float64)
    # Definir el arreglo de signos
    sgn = np.array([1, -1], dtype=np.int32) 
    # Inicializar el array maxwell a cero
    maxwell = np.zeros((nvr, nvx, nx), dtype=np.float64)
    
    # Iterar sobre todos los valores de k
    for k in range(nx):
        # Si Tmaxwell[k] es mayor que 0.0, entonces hacer:
        if Tmaxwell[k] > 0.0:
            # Iterar sobre todos los valores de i
            for i in range(nvr):
                arg = -(vr[i]**2 + (vx - vx_shift[k]/vth)**2)*(mol*Tnorm/Tmaxwell[k])
                # Comparación como en IDL, los valores menores a cero, se convierten en 1, los mayores que cero, se convierten en 0
                arg = np.where(arg < 0.0, 1, 0)                
                # Calcular maxwell: Los valores mayyores a -80, se convierten en 1, los menores en 0, y se toma el e^{1 o 0}
                maxwell[i, :, k] = np.exp(np.where(arg > -80, 1, 0)) 
            print(maxwell[:,:,k].shape)
            print(maxwell[:,:,k])
            print(dVx.shape)
            maxwell[:, :, k] /= np.sum(Vr2pidVr * np.sum(maxwell[:, :, k] * dVx, axis=1))#(np.matmul(dVx,maxwell[:, :, k])))

            if Shifted_Maxwellian_Debug:
                # Calculate vx_out1
                print(vx*dVx)
                print((vx*dVx).shape)
                print(maxwell[:,:,k].shape)
                print(maxwell[:,:,k])
                vx_out1 = vth * np.sum(Vr2pidVr * ((np.matmul((vx * dVx),maxwell[:, :, k]))))   #### Vx, vx, vX y VX son la misma variable, IDL es insensible a estos cambios de mayúsculas y minúsculas.
                
                # Compute vr2vx2_ran2
                for i in range(len(vr)):
                    vr2vx2_ran2[i, :] = vr[i]**2 + (vx - vx_out1 / vth)**2
                
                # Calculate T_out1
                weighted_maxwell = np.matmul(dVx,maxwell[:, :, k])
                T_out1 = (mol * mu * mH) * vth2 * np.sum(Vr2pidVr * (vr2vx2_ran2 * weighted_maxwell)) / (3 * q)
                
                # Calculate vth_local
                vth_local = 0.1 * np.sqrt((2 * Tmaxwell[k] * q)/(mol * mu * mH))
                
                # Calculate errors
                Terror = abs(Tmaxwell[k] - T_out1) / Tmaxwell[k]
                Verror = abs(vx_out1 - vx_shift[k]) / vth_local
            
            # Calcular momentos deseados
            WxD = vx_shift[k]
            ED = WxD**2 + 3 * q * Tmaxwell[k] / (mol * mu * mH)

            # Calcular momentos actuales de maxwell
            WxMax = vth * np.sum(Vr2pidVr * (np.matmul((vx * dVx),maxwell[:, :, k])))
            EMax = vth2 * np.sum(Vr2pidVr * (np.matmul(dVx,(vr2vx2_2D * maxwell[:, :, k]))))
            
            # Calcular Nij y agregar ceros alrededor
            Nij = np.zeros((nvr + 2, nvx + 2))  # Crear un array de ceros con dimensiones aumentadas
            Nij[1:nvr+1, 1:nvx+1] = maxwell[:, :, k] * Vol  # Llenar el array interior con los valores de maxwell * vol
  
            # Calcular variaciones de Nij
            Nijp1_vx_Dvx = shift(Nij * Vx_DVx, shift=( 0,-1), mode='wrap')  # Desplazamiento a la izquierda en la segunda dimensión
            Nij_vx_Dvx = Nij * Vx_DVx
            Nijm1_vx_Dvx = shift(Nij * Vx_DVx, shift=( 0, 1), mode='wrap')  # Desplazamiento a la derecha en la segunda dimensión
            Nip1j_vr_Dvr = shift(Nij * Vr_DVr, shift=(-1, 0), mode='wrap')  # Desplazamiento hacia arriba en la primera dimensión
            Nij_vr_Dvr = Nij * Vr_DVr
            Nim1j_vr_Dvr = shift(Nij * Vr_DVr, shift=( 1, 0), mode='wrap')  # Desplazamiento hacia abajo en la primera dimensión

            _AN =  shift(Nij * Vth_DVx, shift=(0, 1),mode='wrap') - (Nij * Vth_DVx)
            AN[:,:,0] = _AN[1:nvr+1,1:nvx+1]

            _AN = -shift(Nij * Vth_DVx, shift=(0,-1),mode='wrap') + (Nij * Vth_DVx)
            AN[:,:,1] = _AN[1:nvr+1,1:nvx+1]


            BN[:, jpa+1:jpb+1   , 0] =  Nijm1_vx_Dvx[   1:nvr+1 , jpa+2:jpb+2   ] - Nij_vx_Dvx[1:nvr+1  , jpa+2:jpb+2]
            BN[:,   jpa         , 0] = -Nij_vx_Dvx[     1:nvr+1 , jpa+1         ]
            BN[:,   jnb         , 0] =  Nij_vx_Dvx[     1:nvr+1 , jnb+1         ]
            BN[:,   jna:jnb     , 0] = -Nijp1_vx_Dvx[   1:nvr+1 , jna+1:jnb+1   ] + Nij_vx_Dvx[1:nvr+1  , jna+1:jnb+1]
            BN[:,      :        , 0] += Nim1j_vr_Dvr[   1:nvr+1 ,     1:nvx+1   ] - Nij_vr_Dvr[1:nvr+1  ,     1:nvx+1]

            BN[:, jpa+1:jpb+1   , 1] = -Nijp1_vx_Dvx[1:nvr+1, jpa+2:jpb+2] + Nij_vx_Dvx[1:nvr+1,    jpa+2:jpb+2 ]
            BN[:,   jpa         , 1] = -Nijp1_vx_Dvx[1:nvr+1, jpa+1      ]
            BN[:,   jnb         , 1] =  Nijm1_vx_Dvx[1:nvr+1, jnb+1      ]
            BN[:,   jna:jnb     , 1] =  Nijm1_vx_Dvx[1:nvr+1, jna+1:jnb+1] - Nij_vx_Dvx[1:nvr+1,    jna+1:jnb+1 ]
            BN[1:nvr,  :        , 1] -= Nip1j_vr_Dvr[2:nvr+1,     1:nvx+1] - Nij_vr_Dvr[1:nvr+1,        1:nvx+1 ]
            BN[0,      :        , 1] -= Nip1j_vr_Dvr[1      ,     1:nvx+1]


            Nij = Nij[1:nvr+1, 1:nvx+1]

            # Inicializar los arrays de tipo flotante
            TB1 = np.zeros(2, dtype=np.float32)
            TB2 = np.zeros(2, dtype=np.float32)

            # Inicializar la variable
            ia = 0


            while ia < 2:
                TA1 = vth * np.sum(np.matmul(vx,AN[:, :, ia]))
                TA2 = vth2 * np.sum(vr2vx2_2D * AN[:, :, ia])
                ib = 0

                while ib < 2:
                    if TB1[ib] == 0:
                        TB1[ib] = vth * np.sum(np.matmul(vx,BN[:, :, ib]))
                    if TB2[ib] == 0:
                        TB2[ib] = vth2 * np.sum(vr2vx2_2D * BN[:, :, ib])
                    
                    denom = TA2 * TB1[ib] - TA1 * TB2[ib]
                    b_Max = 0.0
                    a_Max = 0.0
                    
                    if denom != 0.0 and TA1 != 0.0:
                        b_Max = (TA2 * (WxD - WxMax) - TA1 * (ED - EMax)) / denom
                        a_Max = (WxD - WxMax - TB1[ib] * b_Max) / TA1
                    
                    if a_Max * sgn[ia] > 0.0 and b_Max * sgn[ib] > 0.0:
                        maxwell[:, :, k] = (Nij + AN[:, :, ia] * a_Max + BN[:, :, ib] * b_Max) / Vol
                        ia = 2
                        ib = 2
                    ib += 1
                ia += 1
            maxwell[:, :, k] = maxwell[:, :, k] / (np.sum(Vr2pidVr * (np.matmul(dVx,maxwell[:, :, k]))))

            if Shifted_Maxwellian_Debug:
                # Calcular vx_out2
                vx_out2 = vth * np.sum(Vr2pidVr * (np.matmul((vx * dVx),maxwell[:, :, k])))
                
                # Calcular vr2vx2_ran2
                vr2vx2_ran2 = np.zeros((nvr, len(vx)))
                for i in range(nvr):
                    vr2vx2_ran2[i, :] = vr[i]**2 + (vx - vx_out2 / vth)**2
                
                # Calcular T_out2
                T_out2 = (mol * mu * mH) * vth2 * np.sum(Vr2pidVr * ((np.matmul(dVx,vr2vx2_ran2 * maxwell[:, :, k])))) / (3 * q)
                
                # Calcular errores
                Terror2 = np.abs(Tmaxwell[k] - T_out2) / Tmaxwell[k]
                Verror2 = np.abs(vx_shift[k] - vx_out2) / vth_local
                
                # Imprimir resultados
                print(f'CREATE_SHIFTED_MAXWELLIAN=> Terror: {Terror:.4f} -> {Terror2:.4f}  Verror: {Verror:.4f} -> {Verror2:.4f}')
    return maxwell
