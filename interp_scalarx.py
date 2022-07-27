import numpy as np
from scipy import interpolate
from warnings import warn

def interp_scalarx(fa,xa,xb,do_warn=0,debug=0):

    # Input function 'a'
	#   fa	=  density np.array
    #   Xa	=  spatial coordinate np.array

    # Desired space coordinates of Output function 'b'
    #   Xb	=  spatial coordinate np.array()

    # do_warn can be set to issue the warnings at the end of the file when needed

    # debug wasn't used in the IDL code, but it was still listed in the inputs
    #   if it really doesn't do anything, it can be removed later

  nxa=xa.size
  
  if fa.size!=nxa:
    raise Exception('Number of elements in fa and Xa do not agree!')
  okk=np.array(np.nonzero(np.logical_and(xb<max(xa), xb>min(xa)))[0])
  if len(okk)==0:
    raise Exception('No values of Xb are within range of Xa')
  
  k0,k1=okk[0],okk[-1]
  nxb=xb.size
  fb=np.zeros(nxb)

  fb[k0:k1+1]=interpolate.interp1d(xa,fa)(xb[okk])

  # warnings only issued if 'warn' parameter has been set
  # for warnings, a message is shown but the code continues
  # Need to check the IDL output to see if a full exception is better
  if do_warn!=0: 
    big=max(abs(xb))
    if k0>0 or k0<nxb-1:
      if k0>0 and abs(fb[k0])>warn*big:
        warn('Non-zero value of fb detected at min(Xa) boundary')
      if k1<nxb-1 and abs(fb[k1])>warn*big:
        warn('Non-zero value of fb detected at max(Xa) boundary')

  return fb
