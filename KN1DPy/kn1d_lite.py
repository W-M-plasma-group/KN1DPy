import numpy as np
from numpy.typing import NDArray
from scipy import interpolate
from dataclasses import dataclass

from .make_dvr_dvx import VSpace_Differentials
from .rates.johnson_hinnov.johnson_hinnov import Johnson_Hinnov
from .kinetic_mesh import KineticMesh
from .kinetic_h import KineticH

from .common import constants as CONST


@dataclass
class KN1DLiteResults:
    '''
    Data structure to hold and pass the results for kn1d_lite.
    Includes mesh information so that runs can be chained.
    '''
    xH: NDArray
    vr: NDArray
    vx: NDArray
    Tnorm: float
    fH: NDArray
    nH: NDArray
    GammaxH: NDArray
    VxH: NDArray
    TH: NDArray
    qxH_total: NDArray
    Sion: NDArray
    fHBC: NDArray
    GammaxHBC: float


def kn1d_lite(
    x, mu, Ti, Te, n, vxi,
    incident_n0,

    # Simple mode
    energies_eV=None,
    velocities_ms=None,
    fractions=None,

    # Advanced mode
    fH_BC=None,

    # Other
    truncate=1e-3,
    max_gen=50,
    compute_errors=False,
    debrief=False,
    debug=False,
    config_path='./config.json',
) -> KN1DLiteResults:
    '''
    Run KN1D atomic neutral transport with a user-specified incident neutral
    boundary condition, bypassing the H2 molecular stage entirely.

    Parameters
    ----------
    x : ndarray(nx)
        Cross-field coordinate (meters).
    mu : float
        1 = hydrogen, 2 = deuterium.
    Ti : ndarray(nx)
        Ion temperature profile (eV).
    Te : ndarray(nx)
        Electron temperature profile (eV).
    n : ndarray(nx)
        Density profile (m^-3).
    vxi : ndarray(nx)
        Plasma velocity profile (m/s, negative is towards wall).
    incident_n0 : float
        Incident neutral density at the boundary (m^-3).  This is the
        *incident* flux only — output nH will be higher due to reflections
        and charge exchange.

    Simple mode (assumes one or more mono-energetic components for incoming neutrals)
    ---------------------------------------------------
    energies_eV : array-like, optional
        Energy of each component in eV.  Default [3.0].
    velocities_ms : array-like, optional
        Speed of each component in m/s.  Overrides energies_eV if given.
    fractions : array-like, optional
        Fraction of incident_n0 in each component (must sum to 1).
        Default [1.0].

    # NOTE
    The default is that all the neutrals enter with 3eV. 
    An input of energies_eV=[3, 10] and fractions=[0.7, 0.3] would mean that 70% of the incident neutrals have 3eV and 30% have 10eV.  
    If velocities_ms is given instead, those speeds are used directly (and energies_eV is ignored).

    Advanced mode (arbitrary distribution function)
    ------------------------------------------------
    fH_BC : ndarray(nvr, nvx), optional
        Boundary distribution function with arbitrary normalisation.
        Scaled internally so that integral = incident_n0.  Negative-vx
        entries are zeroed out before use.

    # NOTE
    The simple mode should suffice for most users. Only use advanced if you want to explore the
    shape of the distribution function on the neutral penetration.

    Other
    -----
    truncate : float
        Convergence threshold.  Default 1e-3.
    max_gen : int
        Maximum collision generations.  Default 50.
    compute_errors : bool
    debrief : bool
    debug : bool
    config_path : str

    Returns
    -------
    KN1DLiteResults
    '''

    prompt = 'KN1D_lite => '

    # ------------------------------------------------------------------
    # Validate / resolve BC inputs
    # ------------------------------------------------------------------

    advanced_mode = fH_BC is not None
    simple_mode = not advanced_mode

    if advanced_mode and (energies_eV is not None or velocities_ms is not None or fractions is not None):
        print(prompt + 'Warning: fH_BC supplied; energies_eV/velocities_ms/fractions are ignored.')

    if simple_mode:
        if velocities_ms is not None and energies_eV is not None:
            print(prompt + 'Warning: both velocities_ms and energies_eV supplied; velocities_ms takes precedence.')

        if velocities_ms is not None:
            component_vs = np.atleast_1d(np.asarray(velocities_ms, dtype=float))
        else:
            if energies_eV is None:
                energies_eV = [3.0]
            energies_eV = np.atleast_1d(np.asarray(energies_eV, dtype=float))
            component_vs = np.sqrt(2.0 * CONST.Q * energies_eV / (mu * CONST.H_MASS))

        if fractions is None:
            fractions = np.ones(len(component_vs)) / len(component_vs)
        fractions = np.atleast_1d(np.asarray(fractions, dtype=float))

        if len(fractions) != len(component_vs):
            raise ValueError(prompt + 'fractions and energies/velocities must have the same length.')
        if not np.isclose(np.sum(fractions), 1.0):
            raise ValueError(prompt + f'fractions must sum to 1.0 (got {np.sum(fractions):.4f}).')

    # ------------------------------------------------------------------
    # Build velocity-space mesh (E0 targets energy bins on component energies)
    # ------------------------------------------------------------------

    if simple_mode:
        if velocities_ms is not None:
            E0 = 0.5 * mu * CONST.H_MASS * component_vs**2 / CONST.Q
        else:
            E0 = energies_eV
    else:
        E0 = np.array([0.0])

    jh = Johnson_Hinnov()

    kh_mesh = KineticMesh('h', mu, x, Ti, Te, n, np.zeros_like(x), jh=jh,
                          E0=E0, config_path=config_path)

    vth = np.sqrt(2.0 * CONST.Q * kh_mesh.Tnorm / (mu * CONST.H_MASS))
    kh_differentials = VSpace_Differentials(kh_mesh.vr, kh_mesh.vx)

    # ------------------------------------------------------------------
    # Build boundary distribution function fHBC and GammaxHBC
    # ------------------------------------------------------------------

    fHBC = np.zeros((kh_mesh.vr.size, kh_mesh.vx.size))

    if simple_mode:
        GammaxHBC = 0.0
        for frac, v_ms in zip(fractions, component_vs):
            v_norm = v_ms / vth
            ix = int(np.argmin(np.abs(kh_mesh.vx - v_norm)))
            fHBC[0, ix] += (frac * incident_n0) / (kh_differentials.dvr_vol[0] * kh_differentials.dvx[ix])
            GammaxHBC += frac * incident_n0 * v_ms
    else:
        # Zero out negative-vx entries
        fH_BC = np.array(fH_BC, dtype=float)
        neg_vx = kh_mesh.vx < 0
        fH_BC[:, neg_vx] = 0.0

        # Scale so density integral equals incident_n0
        current_n = np.sum(kh_differentials.dvr_vol * (fH_BC @ kh_differentials.dvx))
        if current_n <= 0:
            raise ValueError(prompt + 'fH_BC integrates to zero or negative density after zeroing negative vx.')
        scale = incident_n0 / current_n
        fHBC = fH_BC * scale

        # GammaxHBC from the scaled distribution (positive-vx flux)
        pos_vx_flux = np.maximum(kh_mesh.vx, 0.0) * kh_differentials.dvx
        GammaxHBC = vth * np.sum(kh_differentials.dvr_vol * (fHBC @ pos_vx_flux))

    # ------------------------------------------------------------------
    # Interpolate vxi onto atomic mesh
    # ------------------------------------------------------------------

    vxiA = interpolate.interp1d(x, vxi, fill_value='extrapolate')(kh_mesh.x)

    # ------------------------------------------------------------------
    # Set up zero molecular inputs
    # ------------------------------------------------------------------

    fH2A = np.zeros((kh_mesh.vr.size, kh_mesh.vx.size, kh_mesh.x.size))
    fSHA = np.zeros_like(fH2A)
    nHPA = np.zeros(kh_mesh.x.size)
    THPA = np.ones(kh_mesh.x.size)    # avoid divide-by-zero in rate calculations
    fH_init = np.zeros_like(fH2A)

    # ------------------------------------------------------------------
    # Run KineticH
    # ------------------------------------------------------------------

    kinetic_h = KineticH(kh_mesh, mu, vxiA, fHBC, GammaxHBC, jh=jh,
                         ni_correct=True, truncate=truncate, max_gen=max_gen,
                         compute_errors=compute_errors, debrief=debrief, debug=debug,
                         config_path=config_path)

    kh_results = kinetic_h.run_procedure(fH2A, fSHA, fH_init, nHPA, THPA)

    # ------------------------------------------------------------------
    # Pack results
    # ------------------------------------------------------------------

    return KN1DLiteResults(
        xH=kh_mesh.x,
        vr=kh_mesh.vr,
        vx=kh_mesh.vx,
        Tnorm=kh_mesh.Tnorm,
        fH=kh_results.fH,
        nH=kh_results.nH,
        GammaxH=kh_results.GammaxH,
        VxH=kh_results.VxH,
        TH=kh_results.TH,
        qxH_total=kh_results.qxH_total,
        Sion=kh_results.Sion,
        fHBC=fHBC,
        GammaxHBC=GammaxHBC,
    )
