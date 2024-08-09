# import numpy as np
# #   Evaluates the Saha equilibrium population density (m^-3)
# # for atomic hydrogen 
# #
# # Inputs: 
# #   Density     - array, electron density (=hydrogen ion density) (m^-3)
# #   Te          - array, electron temperature (eV)
# #   p           - array, hydrogen energy level, p=1 is ground state

# #   These constants were commented in the original code but they aren't called, 
# # but I am including them anyway 
# # a=h^3/(2*!pi*m*k)^1.5
# #  h=6.626e-34 !j-s
# #  k=1.603e-19 !j/eV
# #  m=0.91094e-30 kg
# #  a=(h/(sqrt(2*!pi*m)*sqrt(k)))^3
# #  a = {j^3 s^3 kg^-1.5 j^-1.5 eV^1.5)
# #  a = {j^1.5 s^3 kg^-1.5 eV^1.5)
# #  j=kg m^2 s^-2
# #  a = {kg^1.5 m^3 s^-3 s^3 kg^-1.5 eV^1.5)
# #  a = {m^3 eV^1.5)
# #  a=3.310e-28

# def NHSaha(Density, Te, p):
#     if len(Density) != len(Te):
#         raise Exception('Number of Elements of Density and Te are different!')
#     if hasattr(p, "__len__"):
#         raise Exception('‘p’ must be a scalar')
#     if p<0:
#         raise Exception('“p” must greater than 0')
#     result = [1.0e32] * len(Density) #----

#     ok = np.array([]) # updated how ok is defined to resolve errors - GG
#     for i in range(0, len(Density)):
#         if 0.0 < Density[i] < 1.0e32 and 0.0 < Te[i] < 1.e32:
#             ok = np.append(ok, i)
#     # converts array from a float array to an int array
#     ok = ok.astype(int) 

#     if len(ok) > 0:
#         for i in ok:
#             result[i] = Density[i] * (3.310E-28 * Density[i]) * p * p * np.exp(13.6057 / (p * p * Te[i])) / (Te[i] ** 1.5)
#             # this returns many infinite values and 
#             # I can't tell if that is an issue with the code inputs or if they are supposed to be that big - GG
#     return result


import numpy as np
from typing import Union

def NHSaha(Density: np.ndarray, Te: np.ndarray, p: int) -> np.ndarray:
    """
    Computes the Saha population density for a given energy level p.
    
    Parameters:
    -----------
    Density : np.ndarray
        Electron density (=hydrogen ion density) in m^-3.
    Te : np.ndarray
        Electron temperature in eV.
    p : int
        Hydrogen energy level, where p=1 corresponds to the ground state.
        
    Returns:
    --------
    np.ndarray
        The Saha population density for the specified level p.

    Raises:
    -------
    ValueError:
        If the number of elements in Density and Te are different, 
        if p is not a scalar, or if p is less than 0.
    
    Notes:
    ------
    This function was adapted from IDL code originally written by B. LaBombard on 6/29/99.
    The Saha population density is calculated based on the formula:
    
    a = h^3 / (2 * π * m * k)^1.5
    where:
        h = Planck's constant = 6.626e-34 J-s
        k = Boltzmann constant = 1.603e-19 J/eV
        m = Electron mass = 0.91094e-30 kg
        a = 3.310e-28 m^3 eV^1.5

    The population density is then calculated as:
    result = Density * (3.310E-28 * Density) * p^2 * exp(13.6057 / (p^2 * Te)) / Te^1.5
    """
    
    if len(Density) != len(Te):
        raise ValueError('Number of elements of Density and Te are different!')
    if not np.isscalar(p):
        raise ValueError('‘p’ must be a scalar')
    if p < 0:
        raise ValueError('“p” must be greater than 0')

    # Initialize the result array with a default value of 1.0e32
    result = np.full(len(Density), 1.0e32)

    # Determine the indices where Density and Te are within the valid range
    conditions = np.logical_and.reduce([
        (0.0 < Density),
        (Density < 1.0e32),
        (0.0 < Te),
        (Te < 1.0e32)
    ])
    
    # Get the indices where all conditions are true
    ok = np.where(conditions)[0]

    if len(ok) > 0:
        # Compute the Saha population density for the valid indices
        result[ok] = Density[ok] * (3.310E-28 * Density[ok]) * p**2 * np.exp(13.6057 / (p**2 * Te[ok])) / (Te[ok]**1.5)

        # Optional: Handle infinite values in the result
        if np.any(np.isinf(result)):
            print("Warning: Infinite values generated in the result")

    return result
