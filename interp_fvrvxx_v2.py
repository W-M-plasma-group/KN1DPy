import numpy as np
from warnings import warn
from scipy.interpolate import RegularGridInterpolator

from Make_dVr_dVx import Make_dVr_dVx
from locate import locate
from sval import sval

from global_vars import mH, q

def compute_thermal_velocities(Tnorma, Tnormb, mu, mH, q):
    Vtha = np.sqrt(2 * q * Tnorma / (mu * mH))
    Vtha2 = Vtha ** 2
    Vthb = np.sqrt(2 * q * Tnormb / (mu * mH))
    Vthb2 = Vthb ** 2
    return Vtha, Vtha2, Vthb, Vthb2

def validate_input_dimensions(fa, Vra, Vxa, Xa):
    print(f'Debug: fa.shape = {fa.shape}, Vra.size = {Vra.size}, Vxa.size = {Vxa.size}, Xa.size = {Xa.size}')
    if fa.shape[0] != Xa.size:
        raise ValueError(f'Number of elements in fa[:,0,0] ({fa.shape[0]}) and Xa ({Xa.size}) do not agree!')
    if fa.shape[1] != Vxa.size:
        raise ValueError(f'Number of elements in fa[0,:,0] ({fa.shape[1]}) and Vxa ({Vxa.size}) do not agree!')
    if fa.shape[2] != Vra.size:
        raise ValueError(f'Number of elements in fa[0,0,:] ({fa.shape[2]}) and Vra ({Vra.size}) do not agree!')

def interp_fvrvxx(fa, Vra, Vxa, Xa, Tnorma, Vrb, Vxb, Xb, Tnormb, do_warn=None, debug=False, correct=True, g=None):
    """
    Interpolates distribution functions used by Kinetic_Neutrals.py,
    Kinetic_H.py, Kinetic_H2.py, and other related procedures.
    
    Parameters:
    fa - np.ndarray (nXa, nVxa, nVra) distribution function
    Vra - np.ndarray (nVra) - radial velocity
    Vxa - np.ndarray (nVxa) - axial velocity
    Xa - np.ndarray (nXa) - spatial coordinate
    Tnorma - float, Normalization temperature for Vra & Vxa
    Vrb - np.ndarray (nVrb) - radial velocity
    Vxb - np.ndarray (nVxb) - axial velocity
    Xb - np.ndarray (nXb) - spatial coordinate
    Tnormb - float, Normalization temperature for Vrb & Vxb
    
    Returns:
    fb - np.ndarray (nVrb, nVxb, nXb) distribution function
    fb is scaled if necessary to make its digital integral over all velocity space equal to that of fa.
    """
    prompt = 'INTERP_FVRVXX => '
    mu = 1

    # Compute thermal velocities
    Vtha, Vtha2, Vthb, Vthb2 = compute_thermal_velocities(Tnorma, Tnormb, mu, mH, q)
    
    # Validate input dimensions
    validate_input_dimensions(fa, Vra, Vxa, Xa)
    
    fV = np.sqrt(Tnormb / Tnorma)
    
    # Find valid indices for interpolation
    oki = np.where((fV * Vrb <= Vra.max()) & (fV * Vrb >= Vra.min()))[0]
    if oki.size < 1:
        raise ValueError('No values of Vrb are within range of Vra')
    i0, i1 = oki[0], oki[-1]

    okj = np.where((fV * Vxb <= Vxa.max()) & (fV * Vxb >= Vxa.min()))[0]
    if okj.size < 1:
        raise ValueError('No values of Vxb are within range of Vxa')
    j0, j1 = okj[0], okj[-1]

    okk = np.where((Xb <= Xa.max()) & (Xb >= Xa.min()))[0]
    if okk.size < 1:
        raise ValueError('No values of Xb are within range of Xa')
    k0, k1 = okk[0], okk[-1]

    # Perform multidimensional interpolation
    interpolator = RegularGridInterpolator((Xa, Vxa, Vra), fa, bounds_error=False, fill_value=0)
    fb_points = np.array(np.meshgrid(Xb, Vxb, Vrb, indexing='ij')).T.reshape(-1, 3)
    fb = interpolator(fb_points).reshape((len(Xb), len(Vxb), len(Vrb)))

    # Normalize fb
    if correct:
        fa_sum = np.sum(fa)
        fb_sum = np.sum(fb)
        if fb_sum != 0:
            fb = fb * (fa_sum / fb_sum)
    
    # Check boundaries
    if do_warn is not None:
        big = np.max(fb)
        for i in range(min(i0, fb.shape[0])):
            if np.any(fb[i, :, :] > do_warn * big):
                warn('Non-zero value of fb detected at min(Vra) boundary')
        for i in range(i1 + 1, min(len(Vrb), fb.shape[0])):
            if np.any(fb[i, :, :] > do_warn * big):
                warn('Non-zero value of fb detected at max(Vra) boundary')
        for j in range(min(j0, fb.shape[1])):
            if np.any(fb[:, j, :] > do_warn * big):
                warn('Non-zero value of fb detected at min(Vxa) boundary')
        for j in range(j1 + 1, min(len(Vxb), fb.shape[1])):
            if np.any(fb[:, j, :] > do_warn * big):
                warn('Non-zero value of fb detected at max(Vxa) boundary')
        for k in range(min(k0, fb.shape[2])):
            if np.any(fb[:, :, k] > do_warn * big):
                warn('Non-zero value of fb detected at min(Xa) boundary')
        for k in range(k1 + 1, min(len(Xb), fb.shape[2])):
            if np.any(fb[:, :, k] > do_warn * big):
                warn('Non-zero value of fb detected at max(Xa) boundary')
    
    return fb
