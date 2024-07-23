from KN1D_while_to_for import KN1D
import scipy.io as sio
from scipy.io import readsav
import numpy as np
print('test','\n')
import matplotlib.pyplot as plt
from scipy import interpolate
import time

#fileloc="C:\\Users\\nholl\\OneDrive\\Documents\\School Work\\Research Info\\ADAS Outputs\\"

testidl=True

def read_sav(file_name):

# This function will print the data from an IDL .sav file
# Inputs:
#   file name - the name of the .sav file
# Outputs:
#   sav_data - a dictionary containing all of the data
    path = file_name
    sav_data = readsav(path)
    return sav_data

fname='1090904024_950to1050'
fname2='test_qcx_h0'
#fname2='No_collisions'

if testidl:
    data_file = read_sav(fname+'.sav')

    keys = data_file.keys() # gets the keys from the dictionary  
    keys_list = list(keys)  # puts the keys into a list 
    #print(keys_list)
    #print(keys_list[0])

    values = data_file.values() # gets the values from the dictionary 
    values_list = list(values)  # puts the values into a list 
    #print(values_list)
    #print('x', values_list[3])
    #print('xlim', values_list[4])
    #print('xsept', values_list[5])
    #print('GuageH2', 50)
    #print('mu', values_list[7])
    #print('vx', values_list[11])
    #print('LC', values_list[12])
    #print('PipeDia', values_list[13])
    #print(keys_list[3])
    #raise Exception('check')
    for i in keys_list:
        #print(i,data_file[i])
        pass

    #interpfunc=interpolate.interp1d(data_file['x'],data_file['n_e'])
    #x=np.linspace(data_file['x'][0],data_file['x'][-1],1000)

    #plt.plot(x,interpfunc(x))
    #plt.xlabel(r'Distance (m)')
    #plt.ylabel(r'Density (m$^{-3}$)')

    #plt.show()
    # x  = values_list[3] # this calls the x np.array from the list of values 
    # print(x)
    # KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
    #         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
    #         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
    #         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
    #         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0)
    #print(data_file['nH2'])
    for i in ['x','x_lim','x_sep','guageH2','mu','t_i','t_e','n_e','vx','lc','d_pipe']:
        #print(i,':',data_file.get(i,50))
        pass
    #print(data_file.keys())
    #print(data_file['vx'])
    #print(data_file['vx1'])
    #print(data_file['vx2'])
    #print(data_file['vx1']-data_file['vx2'])

    t0=time.time()

    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
                xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer=\
                KN1D( data_file['x'], data_file['x_lim'], data_file['x_sep'], data_file['p_wall'], data_file['mu'], data_file['t_i'], data_file['t_e'], \
                data_file['n_e'], data_file['vx'], data_file['lc'], data_file['d_pipe'])

    tf=time.time()

    t=tf-t0

    print('Completed in '+str(t)+' seconds')

    plt.plot(xH2,nH2)
    plt.plot(data_file['xH2'],data_file['nH2'])
    plt.title(r'n$_{H_2}$ Comparison: Shot '+fname[:10])
    plt.yscale('log')
    plt.legend(['Python','IDL'])
    plt.xlabel('x (m)')
    plt.ylabel(r'Density (m$^{-3}$)')
    #plt.savefig(fileloc + fname+'_nH2_testplotb')
    plt.show()
    plt.clf()

    plt.plot(xH,nH)
    plt.plot(data_file['xH'],data_file['nH'])
    plt.yscale('log')
    plt.title(r'n$_{H}$ Comparison: Shot '+fname[:10])
    plt.legend(['Python','IDL'])
    plt.xlabel('x (m)')
    plt.ylabel(r'Density (m$^{-3}$)')
    #plt.savefig(fileloc + fname+'_nH_testplotb')
    plt.show()
    #raise Exception('check')

    #np.savez(fileloc + fname+'b',
    #xH2=xH2, nH2=nH2, GammaxH2=GammaxH2, TH2=TH2, qxH2_total=qxH2_total, nHP=nHP, THP=THP, SH=SH, SP=SP, \
    #        xH=xH, nH=nH, GammaxH=GammaxH, TH=TH, qxH_total=qxH_total, NetHSource=NetHSource, Sion=Sion, QH_total=QH_total, \
    #        SideWallH=SideWallH, Lyman=Lyman, Balmer=Balmer,t=t)

else:
    data_file = read_sav('1090904018_950to1050.sav')

    t0=time.time()
    #raise Exception('check')
    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
                xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer=\
                KN1D( data_file['x'], data_file['x_lim'], data_file['x_sep'], data_file['p_wall']*50, data_file['mu'], data_file['t_i'], data_file['t_e'], \
                data_file['n_e'], data_file['vx'], data_file['lc'], data_file['d_pipe'],\
                adas_rec_h1s=None, adas_ion_h0=None, adas_qcx_h0='qcx#h0_ory#h1')

    tf=time.time()

    t=tf-t0

    print('Completed in '+str(t)+' seconds')

    #np.savez(fileloc + fname2,\
    #xH2=xH2, nH2=nH2, GammaxH2=GammaxH2, TH2=TH2, qxH2_total=qxH2_total, nHP=nHP, THP=THP, SH=SH, SP=SP, \
    #        xH=xH, nH=nH, GammaxH=GammaxH, TH=TH, qxH_total=qxH_total, NetHSource=NetHSource, Sion=Sion, QH_total=QH_total, \
    #        SideWallH=SideWallH, Lyman=Lyman, Balmer=Balmer,t=t)



# runs in 2:46
