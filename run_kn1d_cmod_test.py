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
data_file = './sav_files/cmod_test_in.sav'

# data_file = './sav_files/kn1d_test_inputs.sav'
# data_file = './sav_files/1090904018_950to1050.sav'
# data_file = './sav_files/1090904029_950to1050_towall.sav'
print("Loading file: "  + data_file)
sav_data = readsav(data_file)

sav_data['Ti'] = sav_data['Ti'] * 1e3  # convert from keV to eV
sav_data['Te'] = sav_data['Te'] * 1e3  # convert from keV to eV
sav_data['n'] = sav_data['n'] * 1e20  # convert from 1e20 m^-3 to

##Output
file = 'cmod_test'

print("Beginning KN1D")
start = time.time()
results = kn1d(sav_data['x'], sav_data['xlimiter'], sav_data['xsep'], sav_data['GaugeH2'], sav_data['mu'], sav_data['Ti'], 
               sav_data['Te'], sav_data['n'], sav_data['vxi'], sav_data['LC'], sav_data['PipeDia'], 
               max_gen=100, Hdebug=0, H2debug=0, debrief = 1, Hdebrief = 1, H2debrief = 1, compute_errors = 1)
end = time.time()

print("Elapsed Time: ", end-start)
print()

sys.stdout = standard_out