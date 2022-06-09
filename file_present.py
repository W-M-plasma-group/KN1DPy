#determines whether a file exists

from os.path import exists

def file_present(file):
  return exists(file)
