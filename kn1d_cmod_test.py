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
## .sav --> c-mod | 
data_file = read_sav('1090904018_950to1050.sav')

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
plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH2_cmod_idl.png', dpi = 1000)
plt.show()

plt.figure(constrained_layout=True)
plt.plot(xH, nH, linestyle = 'solid', color = 'blue', label = 'nH')
plt.plot(np.full(10, 0.17), np.linspace(np.min(nH), np.max(nH), 10), linestyle = ':', color = 'orange', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.title('IDL CMOD Density Profile for H')
plt.legend()
plt.savefig('/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH_cmod_idl.png', dpi = 1000)
plt.show()

# KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
#         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
#         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
#         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
#         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0)

xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
            xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D( values_list[3]\
                                                                                                        , values_list[4]\
                                                                                                            , values_list[5]\
                                                                                                                , GuageH2\
                                                                                                                    , values_list[7]\
                                                                                                                        , values_list[8]*10\
                                                                                                                            , values_list[9]*10\
                                                                                                                                , values_list[10]\
                                                                                                                                    , values_list[11]\
                                                                                                                                        , values_list[14]\
                                                                                                                                            , values_list[15])#, \
                                                                                                        #debrief=True, Hdebrief=True, H2debrief=True, debug = True, Hdebug=True, H2debug=True)

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
'''print('vxA', vxA)
print('vxM', vxM)'''

'''np.savetxt('myarray.txt', xH2)
np.savetxt('myarray.txt', nH2)
np.savetxt('myarray.txt', xH)
np.savetxt('myarray.txt', nH)'''

plt.plot(xH2, nH2, linestyle = 'dotted', color = 'red', label = 'Molecular Hydrogen')
plt.xlabel('Position')
plt.ylabel('Density')
plt.yscale("log")
plt.title('Density Profile')
plt.legend()
plt.show()
plt.plot(xH, nH, linestyle = 'dotted', color = 'blue', label = 'Atomic Hydrogen')
plt.yscale("log")
plt.xlabel('Position')
plt.ylabel('Density')
plt.title('Density Profile')
plt.legend()
plt.show()



# new data 
exit()
data_file = read_sav('1090904018_950to1050_towall.sav')

keys = data_file.keys() # gets the keys from the dictionary  
keys_list = list(keys)  # puts the keys into a list 
print(keys_list)
#print(keys_list[0])

values = data_file.values() # gets the values from the dictionary 
values_list = list(values)  # puts the values into a list 
GuageH2 = 0.164
#print(values_list)
print('x', values_list[3])
print('xlim', values_list[4])
print('xsep', values_list[5])
print('GuageH2', GuageH2)
print('mu', values_list[7])
print('Ti', values_list[8])
print('Te', values_list[9])
print('ne', values_list[10])
print('vx', values_list[11])
print('LC', values_list[12])
print('PipeDia', values_list[13])
#print(keys_list[3])

# x  = values_list[3] # this calls the x np.array from the list of values 
# print(x)
# KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
#         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
#         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
#         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
#         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0)

xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
            xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D( values_list[3], values_list[4], values_list[5], GuageH2, values_list[7], values_list[8], values_list[9], values_list[10], values_list[11], values_list[14], values_list[15])#, \
                                                                                                        #debrief=True, Hdebrief=True, H2debrief=True, debug = True, Hdebug=True, H2debug=True)

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
'''print('vxA', vxA)
print('vxM', vxM)'''

'''np.savetxt('myarray.txt', xH2)
np.savetxt('myarray.txt', nH2)
np.savetxt('myarray.txt', xH)
np.savetxt('myarray.txt', nH)'''


plt.plot(xH2, nH2, linestyle = 'dotted', color = 'red', label = 'Molecular Hydrogen')
plt.xlabel('Position')
plt.ylabel('Density')
plt.yscale("log")
plt.title('Density Profile')
plt.legend()
plt.show()
plt.plot(xH, nH, linestyle = 'dotted', color = 'blue', label = 'Atomic Hydrogen')
plt.yscale("log")
plt.xlabel('Position')
plt.ylabel('Density')
plt.title('Density Profile')
plt.legend()
plt.show()