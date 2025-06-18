from .revon import revon
from .revoff import revoff

def printncr(string, rev=None):
  if rev!=None:
    revon()
    print(string, end="")
    revoff()
  else:
    print(string, end="")