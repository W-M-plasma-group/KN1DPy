from Revoff import Revoff
from Revon import Revon

def printcr(string, rev=None):
  if rev!=None:
    Revon()
    print(string)
    Revoff()
  else:
    print(string)