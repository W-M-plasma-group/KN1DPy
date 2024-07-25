from KN1D import KN1D
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.io import readsav
from scipy import interpolate
def read_sav(file_name):

# This function will print the data from an IDL .sav file
# Inputs:
#   file name - the name of the .sav file
# Outputs:
#   sav_data - a dictionary containing all of the data
    path = '//Users/Gwen/Desktop/Plasma_Physics/kn1d/' + file_name
    sav_data = readsav(path)
    return sav_data
# Pull Cmod Data for the LC and PipeDia parameters 
data_file = read_sav('1090904018_950to1050.sav')

keys = data_file.keys() # gets the keys from the dictionary  
keys_list = list(keys)  # puts the keys into a list 
print(keys_list)


values = data_file.values() # gets the values from the dictionary 
values_list = list(values) 

LC = values_list[14]
PipeDia = values_list[15]

# generate empty lists to save data 
nH2_tot = [[], [], [], [], []]
nH_tot = [[], [], [], [], []]
xH2_tot = [[], [], [], [], []]
xH_tot = [[], [], [], [], []]

# generate text file to save data to
save_data = open("kn1d_sparc_temp_sweep_data.txt", "w")
save_data.write("Data for KN1DPy SPARC Temperature sweep: \n")
save_data.close()

for i in range(1, 6):

    # synthetic data for temperature 
    a = i
    b = a + a/20

    CORE = 0.58 - 0.4
    WALL = 0.58  - 0.58
    X_sep = 0.58  - 0.57 

    WIDTH = 0.01
    XSYM = X_sep + WIDTH/2 
    XKNEE = XSYM + WIDTH/2

    SLOPE = 5*a

    pts = 61
    x = np.linspace(WALL, CORE, pts)
    SP = round(((CORE - XKNEE) / (CORE - WALL) * pts)) 

    x_ = np.linspace(XKNEE, CORE, SP)
    y = a * np.tanh(-(2 * (XSYM - x) / WIDTH)) + b
    y_ = y - np.append( np.zeros(pts - SP), SLOPE * (XKNEE - x_))

    Ti = y_ * 1000
    Te = y_ * 1000

    # Synthetic data for density 
    a = 1
    b = a + a/20

    CORE = 0.58 - 0.4
    WALL = 0.58 - 0.58
    X_sep = 0.58 - 0.57 

    WIDTH = 0.01
    XSYM = X_sep + WIDTH/2 
    XKNEE = XSYM + WIDTH/2

    SLOPE = 5*a

    pts = 61
    x = np.linspace(WALL, CORE, pts)
    SP = round(((CORE - XKNEE) / (CORE - WALL) * pts)) 

    print(SP)
    x_ = np.linspace(XKNEE, CORE, SP)
    y = a * np.tanh(-(2 * (XSYM - x) / WIDTH)) + b
    y_ = y - np.append( np.zeros(pts - SP), SLOPE * (XKNEE - x_))

    ne = y_ * 10**20

    xlim = WALL
    xsep = X_sep
    GuageH2 = 0.100
    mu = 2
    vx = np.zeros(61)

    print('x', x.size, x)
    print('xlim', xlim)
    print('xsep', xsep)
    print('GuageH2', GuageH2)
    print('mu', mu)
    print('Ti', Ti.size, Ti)
    print('Te', Te.size, Te)
    print('ne', ne.size, ne)
    print('vx', vx.size, vx)
    print('LC', LC.size, LC)
    print('PipeDia', PipeDia.size, PipeDia)

    # KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
    #         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
    #         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
    #         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
    #         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0)

    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
                xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D( x, xlim, xsep, GuageH2, mu, Ti, Te, ne, vx, LC, PipeDia)#, \
    # save data to lists 
    nH2_tot[i - 1] = nH2.tolist()
    nH_tot[i - 1] = nH.tolist()
    xH2_tot[i - 1] = xH2.tolist()
    xH_tot[i - 1] = xH.tolist()

    # calculate gradient 
    f_H2 = np.array([xH2, nH2])
    f_H = np.array([xH, nH])

    grad_H2 = np.gradient(f_H2)
    grad_H = np.gradient(f_H)
    print(grad_H2)

    #normalized gradient 
    nH2_sep = interpolate.interp1d(xH2, nH2, fill_value='extrapolate')(xsep)
    ngrad_H2 = grad_H2/nH2_sep

    nH_sep = interpolate.interp1d(xH, nH, fill_value='extrapolate')(xsep)
    ngrad_H = grad_H/nH_sep

    # save data into txt file
    save_data = open("kn1d_sparc_temp_sweep_data.txt", "a")

    data = [f"xH2_{i}: {xH2} \n", f"nH2_{i}: {nH2} \n", f"xH_{i}: {xH} \n", f"nH_{i}: {nH} \n"]
    save_data.writelines(data)

    save_data.writelines("Gradient Data \n")

    grad_data = [f"grad_H2{i}: {grad_H2} \n", f"grad_H{i}: {grad_H} \n", f"ngrad_H2{i}: {ngrad_H2} \n", f"ngrad_H{i}: {ngrad_H} \n"]
    save_data.writelines(grad_data)
    save_data.close()
                                                                               
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

print('Press any key to continue')
input()
print('continuing')

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
    plt.plot(xH2_tot[i - 1], nH2_tot[i - 1], linestyle = 'solid', label = f'input scaled by {i}')
plt.plot(np.full(10, 0.01), np.linspace(min(nH2_mins), max(nH2_maxs), 10), linestyle = ':', color = 'green', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.ylim(min(nH2_mins), max(nH2_maxs))
plt.xlim(min(xH2_mins), max(xH2_maxs))
plt.title(f'Python SPARC Density Profile for H$_2$ (temp amp)')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH2_sparc_temp_py_HR.png', dpi = 800)
plt.show()

print('Press any key to continue')
input()
print('continuing')

plt.figure(constrained_layout=True)
for i in range(1, 6):# this was zero but i changed it make sure that nothing major changes 
    plt.plot(xH_tot[i - 1], nH_tot[i - 1], linestyle = 'solid', label = f'input scaled by {i}')
plt.plot(np.full(10, 0.01), np.linspace(min(nH_mins), max(nH_maxs), 10), linestyle = ':', color = 'orange', label = 'sep')
plt.xlabel('Distance from Wall (m)')
plt.ylabel('Density (m$^{-3}$)')
plt.yscale("log")
plt.ylim(min(nH_mins), max(nH_maxs))
plt.xlim(min(xH_mins), max(xH_maxs))
plt.title(f'Python SPARC Density Profile for H (temp amp)')
plt.legend()
plt.savefig(f'/Users/Gwen/Desktop/KN1DPy-Nick/Plots/nH_sparc_temp_py_HR.png', dpi = 800)
plt.show() 
