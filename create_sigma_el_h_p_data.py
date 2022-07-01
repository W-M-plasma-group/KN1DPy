import numpy as np
from make_sigmav import make_sigmav
from sval import sval

def create_sigmav_el_p_h_data():

    #Creates a 2D SigmaV data table in particle energy and target temperature
    #   and saves it as a save set.
    
    sigma_function='Sigma_EL_P_H'

    mE=50
    Emin, Emax = 0.1, 2e4
    nT=50
    Tmin, Tmax = 0.1, 2e4
    
    E_Particle=10**(np.log10(Emin)+(np.log10(Emax)-np.log10(Emin))*np.arange(mE)/(mE-1))
    Mu_Particle=1
    T_Target=10**(np.log10(Tmin)+(np.log10(Tmax)-np.log10(Tmin))*np.arange(nT)/(nT-1))
    Mu_Target=1
    
    sigmav=np.zeros((mE,nT))
    
    for i in range(T_Target.size):
        print('Processing T='+sval(T_Target[i]))
        sigmav[i,:]=make_sigmav(E_Particle,Mu_Particle,T_Target[i],Mu_Target,sigma_function)

    Ln_E_Particle=np.log(E_Particle)
    Ln_T_Target=np.log(T_Target)

    print('Saving results in file: /home/labombard/edge/neutrals/sigmav_el_p_h_data.dat')

    #   SAVE FILE HERE