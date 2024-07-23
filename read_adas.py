import numpy as np
from scipy import optimize
from scipy import interpolate
from string import whitespace

# ADAS has procedures for reading from their data files
# This is an abridged code to cut down on external libraries

settings={}

def is_empty(s):
    for i in s:
        if not i in whitespace:
            return False
    return True

def readfile(name,reaction):
    #with open(settings[name]) as file:
    with open(name+'.dat') as file:
        rawdat=file.read()
    lines=[i for i in rawdat.split('\n') if not is_empty(i)]

    if reaction[:3]=='rec': # log fit?
        n0=int(lines[1].split()[-1])
        nlvl=int(lines[5+n0].split()[-1])
        nid=0
        for i in range(nlvl):
            if lines[9+n0+i].split()[1][:2].lower() in name[4:]:
                nid=i
        Te=lines[11+n0+nlvl].replace('D','e')
        Te=np.array([float(i) for i in Te.split()[2:]])
        sigv=lines[13+n0+nlvl+nid].replace('D','e')
        sigv=np.array([float(i) for i in sigv.split()[1:]])
        

    elif reaction[:3]=='ion': 
        Te=lines[2]+lines[3]+lines[4]+lines[5]
        Te=np.array([float(i) for i in Te.replace('D','e').split()])
        sigv=lines[6]+lines[7]+lines[8]+lines[9]
        sigv=np.array([float(i) for i in sigv.replace('D','e').split()])

    elif reaction[:3]=='qcx': # tanh fit?
        Te=lines[4]
        Te=np.array([float(i) for i in Te.split()[:-3]])
        sig=lines[10].replace('D','e') # total sigma
        sig=np.array([float(i) for i in sig.split()[:-4]])
        return Te,sig

    return Te,sigv

def interp_sig(e,ref_file,reaction):
    eref,sig=readfile(ref_file,reaction)
    interpfunc= interpolate.PchipInterpolator(eref,sig)
    return interpfunc(e)
