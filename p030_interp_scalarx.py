import numpy as np
from scipy.interpolate import interp1d
def interp_scalar_x( f_a: np.ndarray
                    ,x_a: np.ndarray
                    # ,f_b: np.ndarray
                    ,x_b: np.ndarray
                    ,debug   = 0    
                    ,correct = 1
                    ,warn: float = None # float, acceptable truncation level.
                    ):
    nx_a = len(x_a)

    if len(f_a) != nx_a:
        raise ValueError('Number of elements in fa and x_a do not agree!')

    okk = np.where((x_b <= np.max(x_a)) & (x_b >= np.min(x_a)))[0]
    if len(okk) < 1:
        raise ValueError('No values of x_b are within range of x_a')
    k0 = okk[ 0]
    k1 = okk[-1]

    nx_b = len(x_b)
    f_b  = np.zeros(nx_b)
    interpolator = interp1d(x_a, f_a, kind='linear', fill_value="extrapolate")
    f_b[k0:k1 + 1] = interpolator(x_b[okk])

    if warn is not None:
        big = np.max(np.abs(f_b))
        if (k0 > 0) or (k1 < nx_b-1):
            if (k0 > 0)      and (np.abs(f_b[k0]) > warn*big):
                print('Non-zero value of fb detected at min(x_a) boundary')
            if (k1 < nx_b-1) and (np.abs(f_b[k1]) > warn*big):
                print('Non-zero value of fb detected at max(x_a) boundary')
    
    return f_b
