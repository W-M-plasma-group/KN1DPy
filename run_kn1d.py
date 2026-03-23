from dataclasses import dataclass, asdict

from KN1DPy.kn1d import kn1d
from scipy.io import readsav
import numpy as np
import sys

import time

import sys
np.set_printoptions(linewidth=225)
np.set_printoptions(threshold=sys.maxsize)

standard_out = sys.stdout

##Input

data_file = './sav_files/kn1d_test_inputs.sav'
# data_file = './sav_files/1090904018_950to1050.sav'
# data_file = './sav_files/1090904029_950to1050_towall.sav'
print("Loading file: "  + data_file)
sav_data = readsav(data_file)

##Output

print("Beginning KN1D")
start = time.time()
results = kn1d(sav_data['x'], sav_data['x_lim'], sav_data['x_sep'], sav_data['p_wall'], sav_data['mu'], sav_data['t_i'], 
               sav_data['t_e'], sav_data['n_e'], sav_data['vx'], sav_data['lc'], sav_data['d_pipe'], 
               max_gen=100, Hdebug=0, H2debug=0, debrief = 1, Hdebrief = 1, H2debrief = 1, compute_errors = 1)
end = time.time()

print("Elapsed Time: ", end-start)
print()

#print result data
output = open('Results/output.txt', 'w')
sys.stdout = output

for key, value in asdict(results).items():
    print(key)
    print(value)
    print()

output.close()
sys.stdout = standard_out