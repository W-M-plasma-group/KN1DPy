import numpy as np
from warnings import warn

from .make_dvr_dvx import VSpace_Differentials
from .create_shifted_maxwellian import compensate_distribution
from .utils import sval, locate, Bound, interp_1d
from .kinetic_mesh import KineticMesh

def _get_interpolation_bounds(a, b, a_name="a", b_name="b"):
    '''
    Internal method for usage in interp_fvrvxx
    Generates bounds for where b is within the range of a

    Parameters
    ----------
        a : ndarray
            array that determines bounds of b
        b : ndarray
            array being bounded
        a_name, b_name : str (optional)
            variable names for exception message

    Returns
    -------
        tuple : (int, int)
            start, end of the interpolation range
    
    Raises
    -------
        Exception if interpolation range is 0
    '''

    ii = np.where((min(a) <= b) & (b <= max(a)))[0]
    if len(ii) < 1:
        raise Exception(f"No values of {b_name} are within range of {a_name}")
    return Bound(ii[0], ii[-1])


def _test_bounds(fb, test_bound : Bound, var_len, test_axis, iter_bound1 : Bound, iter_bound2 : Bound, do_warn, var_name="a"):
    '''
    Internal method for usage in interp_fvrvxx
    Tests boundaries for vr, vx, x

    Parameters
    ----------
        fb : ndarray
            3D array, distribution function
        test_bound : Bound
            boundary being tested
        var_len : int
            Length of the variable distribution whose bound is being tested
        test_axis : int
            Determines axis for slicing fb, 0, 1, or 2
        iter_bound1, iter_bound2 : Bound
            Boundaries being tested over
        do_warn : float
            Acceptable truncation level
        var_name : str (optional)
            variable name for warning message
    
    Issues
    -------
        Warning if 0 found at boundary edge
    '''

    big = np.max(fb)
    start_error = 0
    end_error = 0
    if (test_bound.start > 0) or (test_bound.end < var_len-1):
        iter_slice1 = iter_bound1.slice(0,1)
        iter_slice2 = iter_bound2.slice(0,1)
        if test_axis == 0:
            min_slice = fb[test_bound.start, iter_slice1, iter_slice2]
            max_slice = fb[test_bound.end, iter_slice1, iter_slice2]
        elif test_axis == 1:
            min_slice = fb[iter_slice1, test_bound.start, iter_slice2]
            max_slice = fb[iter_slice1, test_bound.end, iter_slice2]
        elif test_axis == 2:
            min_slice = fb[iter_slice1, iter_slice2, test_bound.start]
            max_slice = fb[iter_slice1, iter_slice2, test_bound.end]
        else:
            raise Exception("Invalid test axis")
        
        if (start_error == 0) and (test_bound.start > 0) and np.any(min_slice > do_warn*big):
            warn(f"Non-zero value of fb detected at min({var_name}) boundary")
        if (end_error == 0) and (test_bound.end < var_len-1) and np.any(max_slice > do_warn*big):
            warn(f"Non-zero value of fb detected at max({var_name}) boundary")


def interp_fvrvxx(fa: np.ndarray, mesh_a : KineticMesh, mesh_b : KineticMesh, do_warn=None, debug=False, correct=1):
    '''
    Interpolates distribution functions used by kinetic neutral procedures

    Parameters
    ----------
        fa : ndarray
            Input distribution function, 3D array of shape (vra, vxa, xa)
        mesh_a : KineticMesh
            Mesh information for input distribution
        mesh_b : KineticMesh
            Mesh information for desired output distribution
        do_warn : float or None (optional)
            Accebtable truncation level. If None, warnings are not generated
                For interpolations outside the phase space set by
                (Vra, Vxa, Xa), the values of fb are set to zero.
                This may not be acceptable. A test is performed on
                fb at the boundaries. If fb at the boundaries is greater
                than do_warn times the maximum value of fb,
                a warning message is generated.
        debug : bool
            If True, generate debug statements
            
    Returns
    -------
        fb : ndarray
            Interpolated distribution function, scaled if necessary to make its 
            digital integral over all velocity space equal to that of fa
            3D array of of shape (vrb, vxb, xb)
    '''

    prompt = 'INTERP_FVRVXX => '

    v_scale = np.sqrt(mesh_b.Tnorm / mesh_a.Tnorm) # velocity ratio (scales velocities from mesh_a to mesh_b)
    
    # Check shape agreement for fa
    if fa.shape != (mesh_a.vr.size, mesh_a.vx.size, mesh_a.x.size):
        raise Exception('fa (' + str(fa.shape) + ') does not have shape (vra, vxa, xa)' + str((mesh_a.vr.size, mesh_a.vx.size, mesh_a.x.size)))


    # --- Get interpolation Bounds ---

    get_range = lambda a, b : np.where((min(a) <= b) & (b <= max(a)))[0]

    vr_bound = _get_interpolation_bounds(mesh_a.vr, v_scale*mesh_b.vr, "Vra", "Vrb")
    vx_bound = _get_interpolation_bounds(mesh_a.vx, v_scale*mesh_b.vx, "Vxa", "Vxb")
    x_bound = _get_interpolation_bounds(mesh_a.x, mesh_b.x, "Xa", "Xb")

    fb = np.zeros((mesh_b.vr.size, mesh_b.vx.size, mesh_b.x.size))


    # --- Generate differentials ---

    vdiff_a = VSpace_Differentials(mesh_a.vr, mesh_a.vx)
    vdiff_b = VSpace_Differentials(mesh_b.vr, mesh_b.vx)
    

    # --- Compute Weights ---

    if debug:
        print(prompt+'computing new weight')

    # NOTE Removed saving weights temporarily, re-implement later

    # NOTE This is slightly more confusing than the original method, but should be more efficient
    # Set area contributions to Weight array
    # Get arrays of element-wise min/max values for vr and vx, comparing mesh_a and mesh_b
    vr_min = np.maximum(v_scale*vdiff_b.vr_left_bound[:, np.newaxis, np.newaxis, np.newaxis],
                                vdiff_a.vr_left_bound[np.newaxis, np.newaxis, :, np.newaxis])
    vr_max = np.minimum(v_scale*vdiff_b.vr_right_bound[:, np.newaxis, np.newaxis, np.newaxis],
                                vdiff_a.vr_right_bound[np.newaxis, np.newaxis, :, np.newaxis])
    
    vx_min = np.maximum(v_scale*vdiff_b.vx_left_bound[np.newaxis, :, np.newaxis, np.newaxis],
                                vdiff_a.vx_left_bound[np.newaxis, np.newaxis, np.newaxis, :])
    vx_max = np.minimum(v_scale*vdiff_b.vx_right_bound[np.newaxis, :, np.newaxis, np.newaxis],
                                vdiff_a.vx_right_bound[np.newaxis, np.newaxis, np.newaxis, :])

    # Calculate weights
    condition = (vr_max > vr_min) & (vx_max > vx_min)
    weight_value = 2*np.pi*(vr_max**2 - vr_min**2)*(vx_max - vx_min) / (vdiff_b.dvr_vol[:, np.newaxis, np.newaxis, np.newaxis]*vdiff_b.dvx[np.newaxis, :, np.newaxis, np.newaxis])
    weight = np.where(condition, weight_value, 0)

    # Convert to 2D
    weight = np.reshape(weight, (mesh_b.vr.size*mesh_b.vx.size, mesh_a.vr.size*mesh_a.vx.size), order = 'F')


    # --- Correct fb so that it has the same Wx and E moments as fa ---

    if correct:

        # --- Compute Desired Moments ---

        # Determine fb distribution on mesh_a.x grid from weight array
        fa_reshaped = np.reshape(fa, (mesh_a.vr.size*mesh_a.vx.size, mesh_a.x.size), order = 'F')
        fb_on_xa = np.matmul(weight, fa_reshaped)

        #   Compute desired vx_moment and energy_moments of fb, but on the xa grid
        vx_moment_on_xa = np.zeros(mesh_a.x.size)
        energy_moment_on_xa = np.zeros(mesh_a.x.size)

        for k in range(mesh_a.x.size):
            density_a = np.sum(vdiff_a.dvr_vol*(np.matmul(fa[:,:,k], vdiff_a.dvx)))
            if density_a > 0:
                vx_moment_on_xa[k] = np.sqrt(mesh_a.Tnorm)*np.sum(vdiff_a.dvr_vol*(np.matmul(fa[:,:,k], (mesh_a.vx*vdiff_a.dvx)))) / density_a
                energy_moment_on_xa[k] = mesh_a.Tnorm*np.sum(vdiff_a.dvr_vol*(np.matmul((vdiff_a.vmag_squared*fa[:,:,k]), vdiff_a.dvx))) / density_a

        # Compute desired moments on xb grid
        target_vx = np.zeros(mesh_b.x.size)
        target_energy = np.zeros(mesh_b.x.size)

        for k in range(x_bound.start, x_bound.end+1):
            position = np.maximum(locate(mesh_a.x, mesh_b.x[k]), 0)
            kr = np.minimum(position+1, mesh_a.x.size-1)
            kl = np.minimum(position, kr-1)

            interp_fraction = (mesh_b.x[k] - mesh_a.x[kl]) / (mesh_a.x[kr] - mesh_a.x[kl])
            fb[:,:,k] = np.reshape((fb_on_xa[:,kl] + interp_fraction*(fb_on_xa[:,kr] - fb_on_xa[:,kl])), fb[:,:,k].shape, order='F')
            target_vx[k] = vx_moment_on_xa[kl] + interp_fraction*(vx_moment_on_xa[kr] - vx_moment_on_xa[kl])
            target_energy[k] = energy_moment_on_xa[kl] + interp_fraction*(energy_moment_on_xa[kr] - energy_moment_on_xa[kl])

        #   Process each spatial location
        for k in range(mesh_b.x.size):
            if target_energy[k] is None:
                continue

            #   Compute nb, Wxb, and Eb - these are the current moments of fb

            nb = np.sum(vdiff_b.dvr_vol*(np.matmul(fb[:,:,k], vdiff_b.dvx)))
            if nb <= 0:
                continue

            while True:
                
                # --- Adjust fb for desired weights ---

                fb[:,:,k], s = compensate_distribution(fb[:,:,k], vdiff_b, mesh_b.vr, mesh_b.vx, np.sqrt(mesh_b.Tnorm), target_vx[k], target_energy[k], nb=nb, assume_pos=False)
                if s >= 1:
                    break


    # --- Test Boundaries ---

    if do_warn != None:
        # vr_bound
        _test_bounds(fb, vr_bound, mesh_b.vr.size, 0, vx_bound, x_bound, do_warn, var_name="Vra")
        # vx_bound
        _test_bounds(fb, vx_bound, mesh_b.vx.size, 1, vr_bound, x_bound, do_warn, var_name="Vxa")
        # x_bound
        _test_bounds(fb, x_bound, mesh_b.x.size, 2, vr_bound, vx_bound, do_warn, var_name="Xa")


    # --- Rescale ---

    tot_a = np.zeros(mesh_a.x.size)
    for k in range(mesh_a.x.size):
        tot_a[k] = np.sum(vdiff_a.dvr_vol*(np.matmul(fa[:,:,k], vdiff_a.dvx)))
        
    tot_b = np.zeros(mesh_b.x.size)
    tot_b[x_bound.slice(0,1)] = interp_1d(mesh_a.x, tot_a, mesh_b.x[x_bound.slice(0,1)], fill_value="extrapolate")

    ii = np.where(fb > 0)
    if ii[0].size > 0:
        min_tot = np.min(np.array(fb[ii]))
        for k in x_bound.range():
            tot = np.sum(vdiff_b.dvr_vol*(np.matmul(fb[:,:,k], vdiff_b.dvx)))
            if tot > min_tot:
                if debug:
                    print(prompt + 'Density renormalization factor =' + sval(tot_b[k] / tot))
                fb[:,:,k] = fb[:,:,k]*tot_b[k]/tot


    # NOTE Plotting/Debugging was here in the original code
    # May be added later, but has been left out for now

    return fb
