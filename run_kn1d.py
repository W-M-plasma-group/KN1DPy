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
               sav_data['t_e'], sav_data['n_e'], sav_data['vx'], sav_data['lc'], sav_data['d_pipe'], debug = 1, debrief = 1)
end = time.time()

print("Elapsed Time: ", end-start)
print()

#print result data
output = open('Results/output.txt', 'w')
sys.stdout = output

for key, value in results.items():
    print(key)
    print(value)
    print()

output.close()
sys.stdout = standard_out

#Create Results Plots
plt.plot(np.arange(0, len(results["xH2"])), results["xH2"])
plt.title("xH2")
plt.savefig('Results/xH2.png')
plt.clf()

plt.plot(np.arange(0, len(results["nH2"])), results["nH2"])
plt.title("nH2")
plt.savefig('Results/nH2.png')
plt.clf()

plt.plot(np.arange(0, len(results["xH"])), results["xH"])
plt.title("xH")
plt.savefig('Results/xH.png')
plt.clf()

plt.plot(np.arange(0, len(results["nH"])), results["nH"])
plt.title("nH")
plt.savefig('Results/nH.png')
plt.clf()

plt.plot(results["nH2"], results["xH2"])
plt.title("nH2 vs xH2")
plt.savefig('Results/nH2_xH2.png')
plt.clf()

plt.plot(results["nH"], results["xH"])
plt.title("nH vs xH")
plt.savefig('Results/nH_xH.png')
plt.clf()

plt.plot(results["xH"], results["NetHSource"])
plt.title("Net H Source")
plt.savefig('Results/NetHSource.png')
plt.clf()