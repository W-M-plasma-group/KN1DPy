from .revoff import revoff
from .revon import revon

def printcr(string, rev=None):
  if rev!=None:
    revon()
    print(string)
    revoff()
  else:
    print(string)