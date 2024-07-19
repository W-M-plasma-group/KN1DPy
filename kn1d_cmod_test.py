from KN1D import KN1D
import scipy.io as sio
from scipy.io import readsav
import numpy as np
print('test')
def read_sav(file_name):

# This function will print the data from an IDL .sav file
# Inputs:
#   file name - the name of the .sav file
# Outputs:
#   sav_data - a dictionary containing all of the data
    path = 'C:\Users\julio\projects_git\kn1d_gl2\kn1d_python\sav_files_2' + file_name
    sav_data = readsav(path)
    return sav_data

data_file = read_sav('1090904018_950to1050.sav')

keys = data_file.keys() # gets the keys from the dictionary  
keys_list = list(keys)  # puts the keys into a list 
print(keys_list)
#print(keys_list[0])

values = data_file.values() # gets the values from the dictionary 
values_list = list(values)  # puts the values into a list 
#print(values_list)
print('x', values_list[3])
print('xlim', values_list[4])
print('xsept', values_list[5])
print('GuageH2', 50)
print('mu', values_list[7])
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

KN1D( values_list[3], values_list[4], values_list[5], 50, values_list[7], values_list[8], values_list[9], values_list[10], values_list[11], values_list[13], values_list[14])



