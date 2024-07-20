import numpy as np

# This simply defines several variables regarding the size of an input.
# It is slighlty different from the IDL version simply because python and IDL have different functions.
# Instead of returning a number based on the type it returns the type, this means we will have to ammend some if statements 
# in other functions.
# Gwendolyn Galleher

def type_of(arg):
    n = np.size(arg)
    id = type(arg)
    nDim = np.ndim(arg)
    shape = np.shape(arg)
    return id
