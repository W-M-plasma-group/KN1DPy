from KN1D import KN1D
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
    path = './sav_files/' + file_name
    sav_data = readsav(path)
    return sav_data

data_file = read_sav('1090904016_950to1050.sav')

keys = data_file.keys() # gets the keys from the dictionary  
keys_list = list(keys)  # puts the keys into a list 
print(keys_list)
#print(keys_list[0])

values = data_file.values() # gets the values from the dictionary 
values_list = list(values)  # puts the values into a list 
GuageH2 = 0.164
#print(values_list)
print('x', values_list[3].size, values_list[3])
print('xlim', values_list[4].size, values_list[4])
print('xsep', values_list[5].size, values_list[5])
print('GuageH2', GuageH2)
print('mu', values_list[7].size, values_list[7])
print('Ti', values_list[8].size, values_list[8] * 5)
print('Te', values_list[9].size, values_list[9] * 5)
print('ne', values_list[10].size, values_list[10])
print('vx', values_list[11].size, values_list[11])
print('LC', values_list[12].size, values_list[12])
print('PipeDia', values_list[13].size, values_list[13])

xH = values_list[keys_list.index('xh')]
nH = values_list[keys_list.index('nh')]
xH2 = values_list[keys_list.index('xh2')]
nH2 = values_list[keys_list.index('nh2')]
plt.figure(constrained_layout=True)
plt.plot(xH2, nH2, linestyle = 'solid', color = 'red', label = 'nH2')
plt.plot(np.full(10, 0.17), np.linspace(np.min(nH2), np.max(nH2), 10), linestyle = ':', color = 'green', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title('IDL CMOD Density Profile for H$_2$')
plt.legend()
#plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH2_cmod_idl.png', dpi = 1000)
plt.show()

plt.figure(constrained_layout=True)
plt.plot(xH, nH, linestyle = 'solid', color = 'blue', label = 'nH')
plt.plot(np.full(10, 0.17), np.linspace(np.min(nH), np.max(nH), 10), linestyle = ':', color = 'orange', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title('IDL CMOD Density Profile for H')
plt.legend()
#plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH_cmod_idl.png', dpi = 1000)
plt.show()

# KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
#         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
#         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
#         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
#         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0)

for i in range(1, 10):
    print(i)
    Te = values_list[keys_list.index('t_e')]
    Ti = values_list[keys_list.index('t_i')]
    ne = values_list[keys_list.index('n_e')]*i
    x = values_list[keys_list.index('x')]
    plt.figure(constrained_layout=True)
    plt.plot(x, ne, linestyle = 'solid', color = 'red', label = 'ne')
    plt.plot(np.full(10, 0.17), np.linspace(np.min(ne), np.max(ne), 10), linestyle = ':', color = 'green', label = 'sep')
    plt.xlabel('Distance from Wall (m)')
    plt.ylabel('Density (m$^{-3}$)')
    plt.yscale("log")
    plt.title('IDL CMOD Density Profile for ne')
    plt.legend()
    #plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH2_cmod_idl.png', dpi = 1000)
    plt.show()

    plt.figure(constrained_layout=True)
    plt.plot(x, Ti, linestyle = 'solid', color = 'blue', label = 'Ti')
    plt.plot(np.full(10, 0.17), np.linspace(np.min(Ti), np.max(Ti), 10), linestyle = ':', color = 'orange', label = 'sep')
    plt.xlabel('Distance from Wall (m)')
    plt.ylabel('Temp keV')
    plt.yscale("log")
    plt.title('IDL CMOD ion Temp Profile')
    plt.legend()
    #plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH_cmod_idl.png', dpi = 1000)
    plt.show()
    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
                xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D( values_list[3], values_list[4], values_list[5], GuageH2, values_list[7], Ti, Te,ne, values_list[11], values_list[14], values_list[15])#, \
                                                                                                            #debrief=True, Hdebrief=True, H2debrief=True, debug = True, Hdebug=True, H2debug=True)
    print('run_KN1D_success')
    plt.figure(constrained_layout=True)
    plt.plot(xH2, nH2, linestyle = 'solid', color = 'red', label = 'nH2')
    plt.plot(np.full(10, 0.17), np.linspace(np.min(nH2), np.max(nH2), 10), linestyle = ':', color = 'green', label = 'sep')
    plt.xlabel('Distance from Wall (m)')
    plt.ylabel('Density (m$^{-3}$)')
    plt.yscale("log")
    plt.ylim(10**16, 10**17)
    plt.title(f'Python CMOD Density Profile for H$_2$ (dens amp by {i})')
    plt.legend()
    plt.savefig(f'nH2_cmod_dens{i}_py.png', dpi = 800)
    plt.show()

    plt.figure(constrained_layout=True)
    plt.plot(xH, nH, linestyle = 'solid', color = 'blue', label = 'nH')
    plt.plot(np.full(10, 0.17), np.linspace(np.min(nH), np.max(nH), 10), linestyle = ':', color = 'orange', label = 'sep')
    plt.xlabel('Distance from Wall (m)')
    plt.ylabel('Density (m$^{-3}$)')
    plt.yscale("log")
    plt.xlim(0.14, 0.21)
    plt.ylim(10**12, 10**16)
    plt.title(f'Python CMOD Density Profile for H (dens amp by {i})')
    plt.legend()
    plt.savefig(f'nH_cmod_dens{i}_py.png', dpi = 800)
    plt.show()