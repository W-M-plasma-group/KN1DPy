import numpy as np
from reverse import reverse

def integ_bl(X,Y,value_only=None,rev=None):
# X, Y are input np.arrays
# value_only returns only the final calculated value if set.
#   -default returns list
# unsure what rev is used for, but if set it reverses X, Y, and the output

  if rev!=None: #reverses arrays if specified
    Y=reverse(Y)
    X=-reverse(X)

  ans=np.zeros(Y.size) #will be the output
  for i in range(1,Y.size):
    ans[i]=ans[i-1]+.5*(X[i]-X[i-1])*(Y[i]+Y[i-1])

  if value_only!=None: #returns only last element
    return ans[-1]
  if rev!=None: #reverses output
    return reverse(ans)
  return ans
