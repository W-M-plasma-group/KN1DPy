from scipy import interpolate
import numpy as np

#Interpolates an inputted 2D data array along a specified path or trajectory.

def path_interp_2d(p,px,py,x,y):

#   p - Inputted 2D data np.array
#   px - independent coordinate cooresponding to first dimension of p (np.array)
#   py - independent coordinate cooresponding to second dimension of p (np.array)
#   x - inputted 1st coordinate of trajectory (vector)
#   y - inputted 2nd coordinate of trajectory (vector)

#   Check for non-monotonic px, py

  dpx=px-np.roll(px,1)
  if np.min(dpx[1:])<0:
    raise Exception('ERROR in PATH_INTERP_2D => PX is non-monotonic!')

  dpy=py-np.roll(py,1)
  if np.min(dpy[1:])<0:
    raise Exception('ERROR in PATH_INTERP_2D => PY is non-monotonic!')

#   compute coordinates normalized to array indices

  ipx=np.arange(px.size)
  ipy=np.arange(py.size)
  ix=interpolate.interp1d(px,ipx)(x)
  iy=interpolate.interp1d(py,ipy)(y)

  return interpolate.interp1d(p,ix)(iy)