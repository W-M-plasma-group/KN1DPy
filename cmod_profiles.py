import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
def read_sav(file_name):

# This function will print the data from an IDL .sav file
# Inputs:
#   file name - the name of the .sav file
# Outputs:
#   sav_data - a dictionary containing all of the data
    path = '//Users/Gwen/Desktop/Plasma_Physics/kn1d/' + file_name
    sav_data = readsav(path)
    return sav_data

data_file = read_sav('1090904016_950to1050.sav')

keys = data_file.keys() # gets the keys from the dictionary  
keys_list = list(keys)  # puts the keys into a list 
print(keys_list)
values = data_file.values() # gets the values from the dictionary 
values_list = list(values)  # puts the values into a list 

x = values_list[keys_list.index('x')]
print(keys_list.index('x'),x)
Ti = values_list[keys_list.index('t_i')]
Te = values_list[keys_list.index('t_e')]
ne = values_list[keys_list.index('n_e')]

plt.figure(constrained_layout=True)
plt.plot(x, ne, linestyle = 'solid', label = f'$n_e$')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title(f'CMOD Electron Density Profile')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/cmod_ne.png', dpi = 600)
plt.show() 

plt.figure(constrained_layout=True)
plt.plot(x, Te, linestyle = 'solid', label = f'$T_e$')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title(f'CMOD Electron Temperature Profile')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/cmod_te.png', dpi = 600)
plt.show() 

plt.figure(constrained_layout=True)
plt.plot(x, Ti, linestyle = 'solid', label = f'$T_i$')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title(f'CMOD Ion Temperature Profile')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/cmod_ti.png', dpi = 600)
plt.show() 