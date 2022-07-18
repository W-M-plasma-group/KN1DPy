import numpy as np
from scipy.io import readsav
import netCDF4 as nc

#   Bellow are two functions to read and save data from an input file into a dictionary. 
# The first file sav_read is used to read and save .sav files and the second file
# nc_read is to read and save .nc files (netCDF).
#
# Gwendolyn Galleher 

def sav_read(sav_path, nc_path):
    # Inputs:
    #   sav_path - the path to the .sav input file
    #   nc_path  - the path to the .nc file you are creating 
    #Ouputs:
    #    input_dict - a dictionary of all inputs from the input file
    sav_data = readsav(sav_path)
    fn = nc_path
    ds = nc.Dataset(fn, 'w', format = 'NETCDF4') 
    for k,v in sav_data.items(): 
        setattr(ds, k, v)
    input_dict = ds.__dict__
    return input_dict

def nc_read(nc_path):
    # Inputs:
    #   nc_path  - the path to the .nc file you are creating 
    #Ouputs:
    #    input_dict - a dictionary of all inputs from the input file
    fn = nc_path
    ds = nc.Dataset(fn, 'w', format = 'NETCDF4') 
    input_dict = ds.__dict__
    return input_dict
