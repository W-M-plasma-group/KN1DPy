import numpy as np

def poly(x, c):

    if type(x) == list:
	    x = np.array(x)
    n=len(c)-1
    y=c[n]
    for i in range( n-1 ,-1, -1 ):
	    y = y * x + c[i]
    return y
