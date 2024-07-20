# removes leading / trailing spaces
# truncates string to a specified length

def sval(s,length=None):
  return str(s).strip()[:length]
