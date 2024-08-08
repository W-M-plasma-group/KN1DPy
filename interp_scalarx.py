import numpy as np
from scipy import interpolate
from warnings import warn

# def interp_scalarx(fa,xa,xb,do_warn=0,debug=0):

#     # Input function 'a'
# 	#   fa	=  density np.array
#     #   Xa	=  spatial coordinate np.array

#     # Desired space coordinates of Output function 'b'
#     #   Xb	=  spatial coordinate np.array()

#     # do_warn can be set to issue the warnings at the end of the file when needed

#     # debug wasn't used in the IDL code, but it was still listed in the inputs
#     #   if it really doesn't do anything, it can be removed later

#   nxa=xa.size
  
#   if fa.size!=nxa:
#     raise Exception('Number of elements in fa and Xa do not agree!')
#   okk=np.array(np.nonzero(np.logical_and(xb<max(xa), xb>min(xa)))[0])
#   if len(okk)==0:
#     raise Exception('No values of Xb are within range of Xa')
  
#   k0,k1=okk[0],okk[-1]
#   nxb=xb.size
#   fb=np.zeros(nxb)

#   fb[k0:k1+1]=interpolate.interp1d(xa,fa)(xb[okk])

#   # warnings only issued if 'warn' parameter has been set
#   # for warnings, a message is shown but the code continues
#   # Need to check the IDL output to see if a full exception is better
#   if do_warn!=0: 
#     big=max(abs(fb))
#     if k0>0 or k0<nxb-1:
#       if k0>0 and abs(fb[k0])>do_warn*big:
#         warn('Non-zero value of fb detected at min(Xa) boundary')
#       if k1<nxb-1 and abs(fb[k1])>do_warn*big:
#         warn('Non-zero value of fb detected at max(Xa) boundary')

#   return fb




import numpy as np
from scipy.interpolate import interp1d
import warnings

def interp_scalarx(fa, xa, xb, do_warn=0, debug=0):
    """
    Interpolates a scalar function from one set of spatial coordinates to another.
    
    Parameters:
    fa (np.array): Density array.
    xa (np.array): Input spatial coordinates.
    xb (np.array): Desired output spatial coordinates.
    do_warn (float, optional): Threshold for issuing warnings. Default is 0.
    debug (int, optional): Debug flag (not used). Default is 0.
    
    Returns:
    np.array: Interpolated density array at the new spatial coordinates.
    """
    if fa.size != xa.size:
        raise ValueError('Number of elements in fa and xa do not agree!')
    
    in_range_mask = np.logical_and(xb < max(xa), xb > min(xa))
    if not np.any(in_range_mask):
        raise ValueError('No values of xb are within the range of xa')
    
    ok_indices = np.nonzero(in_range_mask)[0]
    start_idx, end_idx = ok_indices[0], ok_indices[-1]
    
    interpolator = interp1d(xa, fa, bounds_error=False, fill_value=0)
    fb = interpolator(xb)
    
    if do_warn != 0:
        max_fb_abs = np.max(np.abs(fb))
        if start_idx > 0 and np.abs(fb[start_idx]) > do_warn * max_fb_abs:
            warnings.warn('Non-zero value of fb detected at min(xa) boundary')
        if end_idx < xb.size - 1 and np.abs(fb[end_idx]) > do_warn * max_fb_abs:
            warnings.warn('Non-zero value of fb detected at max(xa) boundary')
    
    return fb

# # Example usage:
# fa = np.array([1, 2, 3])
# xa = np.array([0, 1, 2])
# xb = np.array([-1, 0.5, 1.5, 3])
# fb = interp_scalar(fa, xa, xb, do_warn=0.1)
# print(fb)
