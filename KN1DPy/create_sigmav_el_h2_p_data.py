import numpy as np

from .make_sigmav import make_sigmav
from .sval import sval
from .sigma.sigma_el_p_hh import sigma_el_p_hh # added import for sigma_function - nh

def create_sigmav_el_h2_p_data():

    #Creates a 2D SigmaV data table in particle energy and target temperature
    #   and saves it as a save set.
    
    sigma_function=sigma_el_p_hh

    mE=50
    Emin, Emax = 0.1, 1e3
    nT=50
    Tmin, Tmax = 0.1, 1e3
    
    E_Particle=10**(np.log10(Emin)+(np.log10(Emax)-np.log10(Emin))*np.arange(mE)/(mE-1))
    Mu_Particle=2
    T_Target=10**(np.log10(Tmin)+(np.log10(Tmax)-np.log10(Tmin))*np.arange(nT)/(nT-1))
    Mu_Target=1
    
    sigmav=np.zeros((nT,mE))
    
    for i in range(T_Target.size):
        print('Processing T='+sval(T_Target[i]))
        sigmav[i,:]=make_sigmav(E_Particle,Mu_Particle,T_Target[i],Mu_Target,sigma_function) # fixed indexing - nh

    Ln_E_Particle=np.log(E_Particle)
    Ln_T_Target=np.log(T_Target)

    #print('Saving results in file: /home/labombard/edge/neutrals/sigmav_el_h2_p_data.dat') # orginal file location in IDL version
    print('Saving results in file: sigmav_el_h2_p_data.npz') # saves in current directory for now - nh

    np.savez('sigmav_el_h2_p_data',Ln_E_Particle=Ln_E_Particle,Ln_T_Target=Ln_T_Target,SigmaV=sigmav)
