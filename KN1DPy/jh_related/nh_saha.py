import numpy as np

#   Evaluates the Saha equilibrium population density (m^-3)
# for atomic hydrogen 
#
# Inputs: 
#   Density     - array, electron density (=hydrogen ion density) (m^-3)
#   Te          - array, electron temperature (eV)
#   p           - array, hydrogen energy level, p=1 is ground state

#   These constants were commented in the original code but they aren't called, 
# but I am including them anyway 
# a=h^3/(2*!pi*m*k)^1.5
#  h=6.626e-34 !j-s
#  k=1.603e-19 !j/eV
#  m=0.91094e-30 kg
#  a=(h/(sqrt(2*!pi*m)*sqrt(k)))^3
#  a = {j^3 s^3 kg^-1.5 j^-1.5 eV^1.5)
#  a = {j^1.5 s^3 kg^-1.5 eV^1.5)
#  j=kg m^2 s^-2
#  a = {kg^1.5 m^3 s^-3 s^3 kg^-1.5 eV^1.5)
#  a = {m^3 eV^1.5)
#  a=3.310e-28

def nh_saha(Density, Te, p):
    if len(Density) != len(Te):
        raise Exception('Number of Elements of Density and Te are different!')
    if hasattr(p, "__len__"):
        raise Exception('‘p’ must be a scalar')
    if p<0:
        raise Exception('“p” must greater than 0')
    result = [1.0e32] * len(Density)

    ok = np.array([]) # updated how ok is defined to resolve errors - GG
    for i in range(0, len(Density)):
        if 0.0 < Density[i] < 1.0e32 and 0.0 < Te[i] < 1.e32:
            ok = np.append(ok, i)
    # converts array from a float array to an int array
    ok = ok.astype(int) 

    if len(ok) > 0:
        for i in ok:
            result[i] = Density[i] * (3.310E-28 * Density[i]) * p * p * np.exp(13.6057 / (p * p * Te[i])) / (Te[i] ** 1.5)
            # this returns many infinite values and 
            # I can't tell if that is an issue with the code inputs or if they are supposed to be that big - GG
    return result
