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
    path = '//Users/Gwen/Desktop/Plasma_Physics/kn1d/' + file_name
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
print('Ti', values_list[8].size, values_list[8])
print('Te', values_list[9].size, values_list[9])
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
nH2_tot = [[], [], [], [], []]
nH_tot = [[], [], [], [], []]
xH2_tot = [[], [], [], [], []]
xH_tot = [[], [], [], [], []]
for i in range(1, 6):
    print(i)
    Te = values_list[keys_list.index('t_e')]*i
    Ti = values_list[keys_list.index('t_i')]*i
    ne = values_list[keys_list.index('n_e')]
    x = values_list[keys_list.index('x')]

    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
                xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D( values_list[3], values_list[4], values_list[5], GuageH2, values_list[7], Ti, Te,ne, values_list[11], values_list[14], values_list[15])
    
    nH2_tot[i - 1] = nH2.tolist()
    nH_tot[i - 1] = nH.tolist()
    xH2_tot[i - 1] = xH2.tolist()
    xH_tot[i - 1] = xH.tolist()                                          
    print('xH2', xH2)
    print('nH2', nH2)
    print('GammaxH2', GammaxH2)
    print('TH2', TH2)
    print('qxH2_total', qxH2_total)
    print('nHP', nHP)
    print('THP', THP)
    print('SH', SH)
    print('SP',SP)
    print('xH', xH)
    print('nH', nH)
    print('GammaxH', GammaxH)
    print('TH', TH)
    print('qxH_total', qxH_total)
    print('NetHSource', NetHSource)
    print('Sion', Sion)
    print('QH_total', QH_total)
    print('SideWallH', SideWallH)
    print('Lyman', Lyman)
    print('Balmer', Balmer)                                      
                                          
nH2_mins = []
nH2_maxs = []
xH2_mins = []
xH2_maxs = []
nH_mins = []
nH_maxs = []
xH_mins = []
xH_maxs = []

for i in range(0, 5):
    nH2_mins.append(min(nH2_tot[i]))
    nH2_maxs.append(max(nH2_tot[i]))
    xH2_mins.append(min(xH2_tot[i]))
    xH2_maxs.append(max(xH2_tot[i]))
    nH_mins.append(min(nH_tot[i]))
    nH_maxs.append(max(nH_tot[i]))
    xH_mins.append(min(xH_tot[i]))
    xH_maxs.append(max(xH_tot[i]))
plt.figure(constrained_layout=True)

for i in range(1, 6):
    plt.plot(xH2_tot[i - 1], nH2_tot[i - 1], linestyle = 'solid', label = f'input scaled {i}')
plt.plot(np.full(10, 0.17), np.linspace(min(nH2_mins), max(nH2_maxs), 10), linestyle = ':', color = 'green', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.ylim(min(nH2_mins), max(nH2_maxs))
plt.xlim(min(xH2_mins), 0.175)
plt.title(f'Python SPARC Density Profile for H$_2$ (temp amp)')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH2_cmod_temp_py.png', dpi = 800)
plt.show()
print('Press any key to continue')
input()
print('continuing')

plt.figure(constrained_layout=True)
for i in range(0, 6):
    plt.plot(xH_tot[i - 1], nH_tot[i - 1], linestyle = 'solid', label = f'input scaled {i}')
plt.plot(np.full(10, 0.17), np.linspace(min(nH_mins), max(nH_maxs), 10), linestyle = ':', color = 'orange', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.ylim(min(nH_mins), max(nH_maxs))
plt.xlim(min(xH_mins), max(xH_maxs))
plt.title(f'Python SPARC Density Profile for H (temp amp)')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH_cmod_temp_py.png', dpi = 800)
plt.show()      
                                                                                                       #debrief=True, Hdebrief=True, H2debrief=True, debug = True, Hdebug=True, H2debug=True)
