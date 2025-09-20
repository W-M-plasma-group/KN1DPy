# searches a sorted, non-repeating list
# if list ascends (ex: 1,2,3,4)- 
#    returns index of last item less than input value
#    returns -1 if value less than first item
# if list descends (ex: 4,3,2,1)-
#    retruns index of last item greater than input value
#    returns -1 if value greater than first item
# table must be a np.array, 
# value can be scalar, list, or np.array
# returns a list of equal length to len(value)

import numpy as np

def locate(table, value):
  if table.size == 0:
    return -1
  out = [] # this is the output
  if type(value)!=list and type(value)!=np.ndarray:
    value=[value] # converts value to list if it is not already a list or np.array
  if table[0]<table[1]: # if table is ascending
    for i in value:
      if i>=table[-1]: # case that value is greater than the last item
        out+= [table.size-1]
      else:
        out+= [np.where(table>i)[0][0]-1] # adds index to out
  else: # if table is descending
    for i in value:
      if i<=table[-1]: # case that value is less than the last item
        out+= [table.size-1]
      else:
        out+= [np.where(table<i)[0][0]-1] # adds index to out
  return out
