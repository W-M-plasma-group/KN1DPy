from KN1DPy.kn1d import kn1d
import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
import sys

import time


standard_out = sys.stdout

##Input

data_file = './dev_files/kn1d_test_inputs.sav'
#data_file = './dev_files/1090904018_950to1050.sav'
print("Loading file: "  + data_file)
sav_data = readsav(data_file)

##Output

print("Beginning KN1D")
start = time.time()
results = kn1d(sav_data['x'], sav_data['x_lim'], sav_data['x_sep'], sav_data['p_wall'], sav_data['mu'], sav_data['t_i'], 
               sav_data['t_e'], sav_data['n_e'], sav_data['vx'], sav_data['lc'], sav_data['d_pipe'])
end = time.time()

print("Elapsed Time: ", end-start)
print()

#print result data
output = open('Results/output.txt', 'w')
sys.stdout = output

print("xH2")
print(results[0])
print()

print("nH2")
print(results[1])
print()

print("xH")
print(results[9])
print()

print("nH")
print(results[10])
print()

print("Full Results")
print(results)

output.close()
sys.stdout = standard_out

#Create Results Plots
plt.plot(np.arange(0, len(results[0])), results[0])
plt.title("xH2")
plt.savefig('Results/xH2.png')
plt.clf()

plt.plot(np.arange(0, len(results[1])), results[1])
plt.title("nH2")
plt.savefig('Results/nH2.png')
plt.clf()

plt.plot(np.arange(0, len(results[9])), results[9])
plt.title("xH")
plt.savefig('Results/xH.png')
plt.clf()

plt.plot(np.arange(0, len(results[10])), results[10])
plt.title("nH")
plt.savefig('Results/nH.png')
plt.clf()

plt.plot(results[0], results[1])
plt.title("nH2 vs xH2")
plt.savefig('Results/nH2_xH2.png')
plt.clf()

plt.plot(results[9], results[10])
plt.title("nH vs xH")
plt.savefig('Results/nH_xH.png')
plt.clf()

plt.plot(results[9], results[14])
plt.title("Net H Source")
plt.savefig('Results/NetHSource.png')
plt.clf()