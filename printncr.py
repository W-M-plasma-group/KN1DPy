from Revon import Revon
from Revoff import Revoff

def printncr(string, rev=None):
  if rev!=None:
    Revon()
    print(string, end="")
    Revoff()
  else:
    print(string, end="")