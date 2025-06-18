import numpy as np

def comb(data, npts=2000):
  if data.size>npts:
    return np.arange(npts+1)*(data.size-1)/npts
  return np.arange(data.size)
