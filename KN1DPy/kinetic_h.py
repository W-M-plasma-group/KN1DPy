from dataclasses import dataclass
from typing import Optional
import numpy as np
from numpy.typing import NDArray

from .utils import sval, get_config
from .make_dvr_dvx import VSpace_Differentials
from .create_shifted_maxwellian import create_shifted_maxwellian
from .kinetic_mesh import KineticMesh

from .rates.collrad.collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from .rates.johnson_hinnov.johnson_hinnov import Johnson_Hinnov
from .rates.janev.sigmav_ion_h0 import sigmav_ion_h0
from .rates.janev.sigmav_rec_h1s import sigmav_rec_h1s
from .rates.adas.adas_ionisation import scd_adas, acd_adas
from .rates.janev.sigma_cx_h0 import sigma_cx_h0
from .rates.janev.sigma_el_h_h import sigma_el_h_h
from .rates.janev.sigma_el_h_hh import sigma_el_h_hh
from .rates.janev.sigma_el_p_h import sigma_el_p_h
from .rates.janev.sigmav_cx_h0 import sigmav_cx_h0

from .common import constants as CONST
from .common.Kinetic_H import *


# Dataclasses for use in kinetic_h

@dataclass
class KHCollisions():
    '''
    Collision settings for Kinetic H procedure
    '''
    H2_H_EL: bool = False
    H_H_EL: bool = False
    H_P_EL: bool = False
    H_P_CX: bool = False
    SIMPLE_CX: bool = False

@dataclass
class MeshEqCoefficients():
    '''
    Mesh Equation values used in kinetic_h iteration
    Eqs. (3.22), (3.25), (3.30), (3.33)
    '''
    A: NDArray
    B: NDArray
    C: NDArray
    D: NDArray
    F: NDArray
    G: NDArray

@dataclass
class CollisionType():
    '''
    Data class for grouping H_H, H_P, and H_H2 elastic collision data
    '''
    H_H: NDArray
    H_P: NDArray
    H_H2: NDArray

@dataclass
class KHResults():
    '''
    Variables for results of KineticH.run_procedure()
    See run_procedure for more detail on individual variables
    '''
    fH: NDArray
    nH: NDArray
    GammaxH: NDArray
    VxH: NDArray
    pH: NDArray
    TH: NDArray
    qxH: NDArray
    qxH_total: NDArray
    NetHSource: NDArray
    Sion: NDArray
    QH: NDArray
    RxH: NDArray
    QH_total: NDArray
    AlbedoH: float
    SideWallH: NDArray
    nH_gen0: Optional[NDArray] = None
    nH_generations: Optional[NDArray] = None


class KineticH():
    '''
    This class is part of the "KN1D" atomic and molecular neutral transport code.

    This class contains the data and methods to solve a 1-D spatial, 2-D velocity kinetic neutral transport 
    problem for atomic hydrogen (H) or deuterium by computing successive generations of 
    charge exchange and elastic scattered neutrals. The routine handles electron-impact 
    ionization, proton-atom charge exchange, radiative recombination, and elastic
    collisions with hydrogenic ions, neutral atoms, and molecules.

    The positive vx half of the atomic neutral distribution function is inputted at x(0) 
    (with arbitrary normalization) and the desired flux of hydrogen atoms entering the slab,
    at x(0) is specified. Background profiles of plasma ions, (e.g., Ti(x), Te(x), n(x), vxi(x),...)
    molecular ions, (nHP(x), THP(x)), and molecular distribution function (fH) are inputted.

    Optionally, the hydrogen source velocity distribution function is also inputted.
    (The H source and fH2 distribution functions can be computed using procedure 
    "Kinetic_H2.pro".) The code returns the atomic hydrogen distribution function, fH(vr,vx,x) 
    for all vx, vr, and x of the specified vr,vx,x grid.

    Since the problem involves only the x spatial dimension, all distribution functions
    are assumed to have rotational symmetry about the vx axis. Consequently, the distributions
    only depend on x, vx and vr where vr =sqrt(vy^2+vz^2)

    History:

        B. LaBombard   First coding based on Kinetic_Neutrals.pro 		22-Dec-2000

        For more information, see write-up: "A 1-D Space, 2-D Velocity, Kinetic 
        Neutral Transport Algorithm for Hydrogen Atoms in an Ionizing Plasma", B. LaBombard

    Variable names contain characters to help designate species -
    atomic neutral (H), molecular neutral (H2), molecular ion (HP), proton (i) or (P)
    '''

    # Theta-prime Coordinate
    ntheta = 5 # use 5 theta mesh points for theta integration
    dtheta = np.ones(ntheta) / ntheta
    cos_theta = np.cos(np.pi*(np.arange(ntheta) + 0.5) / ntheta)

    # Internal Print Formatting
    prompt = 'Kinetic_H => '


    def __init__(self, mesh: KineticMesh, mu: int, vxi: NDArray, fHBC: NDArray, GammaxHBC: float, jh: Johnson_Hinnov = None,
                 recomb: bool = True, ni_correct: bool = False, truncate: float = 1e-4, max_gen: int = 100,
                 compute_errors: bool = False, debrief: int = 0, debug: int = 0, config_path: str = './config.toml',
                 return_gen0: bool = False, return_all_generations: bool = False):
        '''
        Parameters
        ----------
            mesh : KineticMesh
                Mesh data for h kinetic procedure, must be of type 'h'
                Includes coordinate data and temperature/density profiles
            mu : int
                1=hydrogen, 2=deuterium
            vxi : ndarray
                flow speed profile (m/s)
            fHBC : ndarray
                2D array, input boundary condition. Specifies shape of atom velocity distribution (fH) at x=0
            GammaxHBC : float
                Desired neutral atom flux density at x=0 (m^-2 s^-1)
            jh_coeffs : JH_Coef, defualt=None
                Common blocks used to pass data for JH methods
            recomb : bool, default=True
                If true, includes recombination as a source of atomic neutrals in the algorithm
            ni_correct : bool, default=False
                If true, Corrects hydrogen ion density according to quasineutrality: ni=ne-nHp
            truncate : float, default=1.0e-4
                Convergence threshold for generations
            max_gen : int, default=50
                Max number of generations
            compute_errors : bool, default=False
                If true, compute error estimates
            debug : int, default=0
                - 0=do not execute debug code
                - 1=summary debug
                - 2=detail debug
                - 3=very detailed debug
            debrief : int, default=0
                - 0=do not print
                - 1=print summary information
                - 2=print detailed information
        '''

        # --- Settings ---

        # Configuration Options
        self.config = get_config(config_path)
        
        col = self.config['collisions']
        self.COLLISIONS = KHCollisions(col['H2_H_EL'], col['H_H_EL'], col['H_P_EL'], col['H_P_CX'], col['SIMPLE_CX'])
        
        self.ion_rate_option = self.config['kinetic_h']['ion_rate']

        # Small numerical tolerances to avoid divide-by-zero in velocity grid
        # spacing and wall pressure calculations respectively.
        self.DeltaVx_tol = 0.01
        self.Wpp_tol = 0.001

        # Internal Debug switches
        self.CI_Test = self.config['kinetic_h']['ci_test']
        self.Do_Alpha_CX_Test = self.config['kinetic_h']['alpha_cx_test']

        # Run Settings
        self.debrief = debrief
        self.truncate = truncate
        self.max_gen = max_gen
        self.ni_correct = ni_correct
        self.compute_errors = (compute_errors and debrief)
        self.recomb = recomb
        self.debug = debug
        self.return_gen0 = return_gen0
        self.return_all_generations = return_all_generations

        # Override settings for debug
        if self.debug > 0:
            self.debrief = np.maximum(self.debrief, 1)

        
        # --- Main Attributes ---

        # Main attributes
        self.mesh = mesh
        self.mu = mu
        self.vxi = vxi
        self.fHBC = fHBC
        self.GammaxHBC = GammaxHBC

        # Shorthand sizes for main mesh variables
        self.nvr = mesh.vr.size
        self.nvx = self.mesh.vx.size
        self.nx = self.mesh.x.size

        self.vx_neg = np.nonzero(self.mesh.vx < 0)[0]
        self.vx_pos = np.nonzero(self.mesh.vx > 0)[0]
        self.vx_zero = np.nonzero(self.mesh.vx == 0)[0]


        # --- Internal Variables ---

        self.vth = np.sqrt((2*CONST.Q*self.mesh.Tnorm) / (self.mu*CONST.H_MASS))

        # Vr^2-2*Vx^2
        self.vr2_2vx2_2D = np.asarray([(vr**2) - 2*(self.mesh.vx**2) for vr in self.mesh.vr])

        # Differential Values
        differential = VSpace_Differentials(self.mesh.vr, self.mesh.vx)
        self.dvr_vol = differential.dvr_vol
        self.dvx = differential.dvx

        # FHBC_Input
        self._init_fhbc_input()


        # Common Blocks
        self.Input = Kinetic_H_Input()
        self.Internal = Kinetic_H_Internal()
        self.Output = Kinetic_H_Output(self.nx)
        self.H2_Moments = Kinetic_H_H2_Moments()
        self.Errors = Kinetic_H_Errors()

        # Setup Johnson-Hinov
        self.jh = jh
        if (jh == None) and (self.ion_rate_option == "jh"):
            self.jh = Johnson_Hinnov()


        self._test_init_parameters()

        # Initial Computations
        # Some may not be used depending on inputs
        self._init_static_internals()

        if self.compute_errors:
            self._compute_vbar_error()

        return
    
    
    def run_procedure(self, fH2: NDArray = None, fSH: NDArray = None, fH: NDArray = None, nHP: NDArray = None, THP: NDArray = None) -> KHResults:
        '''
        Solves a 1-D spatial, 2-D velocity kinetic neutral transport 
        problem for atomic hydrogen or deuterium (H)

        Parameters
        ----------
            fH2 : ndarray, default=None
                3D array, molecular distribution function. If None, H-H2 collisions are not computed
            fSH : ndarray, defualt=None
                Source velocity distribution function. If None, zero array is used
            fH : ndarray, default=None
                3D array, atomic distribution function. If None, zero array is used
            nHP : ndarray, defualt=None
                Molecular ion density profile (m^-3). If None, zero array is used
            THP : ndarray, defualt=None
                Molecular ion temperature profile (m^-3). If None, array of 3.0 used

        Returns
        -------
        KHResults

            fH : ndarray
                3D array, atomic distribution function.
            nH : ndarray, defualt=None
                Neutral atom density profile (m^-3).
            GammaxH : ndarray
                Neutral flux profile (m^-2 s^-1)
            VxH : ndarray
                Neutral velocity profile (m s^-1)
            pH : ndarray
                Neutral pressure (eV m^-2)
            TH : ndarray
                Neutral temperature profile (eV)
            qxH : ndarray
                Neutral random heat flux profile (watts m^-2)
            qxH_total : ndarray
                Total neutral heat flux profile (watts m^-2)
            NetHSource : ndarray
                Net H0 source (m^-3 s^-1)
            Sion : ndarray
                H ionization rate (m^-3 s^-1)
            QH : ndarray
                Rate of net thermal energy transfer into neutral atoms (watts m^-3)
            RxH : ndarray
                Rate of x momentum transfer to neutral atoms (N m^-2)
            QH_total : ndarray
                Net rate of total energy transfer into neutral atoms (watts m^-3)
            AlbedoH : float
                Ratio of atomic particle flux with Vx < 0 divided by particle flux with Vx > 0 at x=0
            SideWallH : ndarray
                Atomic sink rate from interation with 'side walls' (m^-3 s^-1)
        '''

        nvr, nvx, nx = self.nvr, self.nvx, self.nx

        # --- Initialize Inputs ---

        if fH2 is None:
            fH2 = np.zeros((nvr,nvx,nx))
        if fSH is None:
            fSH = np.zeros((nvr,nvx,nx))
        if nHP is None:
            nHP = np.zeros(nx)
        if THP is None:
            THP = np.full(nx, 1.0)
        if fH is None:
            fH = np.zeros((nvr,nvx,nx))
        self._test_input_parameters(fH2, fSH, nHP, THP, fH)

        # If fH2 is zero, then turn off elastic H2 <-> H collisions
        self.COLLISIONS.H2_H_EL = self.config['collisions']['H2_H_EL']
        if np.sum(fH2) <= 0:
            self.COLLISIONS.H2_H_EL = False

        # Scale input molecular distribution function to agree with desired flux
        fH[:,self.vx_pos,0] = self.fHBC_input[:,self.vx_pos]


        # --- Compute Variables---

        self._compute_dynamic_internals(fH, fH2, nHP, THP, fSH)

        # Compute nH
        nH = np.zeros(nx)
        for k in range(nx):
            nH[k] = np.sum(self.dvr_vol*(np.matmul(fH[:,:,k], self.dvx)))

        # Compute Side-Wall collision rate
        gamma_wall = np.zeros((nvr,nvx,nx))
        for k in range(nx):
            if self.mesh.PipeDia[k] > 0:
                gamma_wall[:,:,k] = 2*self.mesh.vr / self.mesh.PipeDia[k]


        # --- Iteration ---
        
        fH, nH, alpha_c, Beta_CX_sum, collision_freqs, m_sums, NHG, igen = self._run_iteration_scheme(fH, nH, gamma_wall)


        # --- Compute Results ---

        results = self._compile_results(fH, nH, fSH, gamma_wall, alpha_c, Beta_CX_sum, collision_freqs, m_sums)

        if self.return_gen0:
            results.nH_gen0 = NHG[:, 0].copy()
        if self.return_all_generations:
            results.nH_generations = NHG[:, :igen + 1].copy()

        if self.compute_errors:
            self._compute_final_errors(results, Beta_CX_sum, m_sums, alpha_c, collision_freqs)

        
        # --- Save input Variables ---

        self.Input.fH2_s = fH2
        self.Input.fSH_s = fSH
        self.Input.nHP_s = nHP
        self.Input.THP_s = THP
        self.Input.fH_s = fH


        self._debrief_msg("Finished", 0)

        return results
        


    # ------ Main Procedure Functions ------

    def _run_iteration_scheme(self, fH, nH, gamma_wall):
        '''
        Implements the main procdure iteration, iterates until fH converges.
        Convergence is determined by evaluating changes in ion density, with iterations
        terminating when density change is low enough.
        '''

        nvr, nvx, nx = self.nvr, self.nvx, self.nx

        #	Set iteration scheme
        fH_iterate = False
        if self.COLLISIONS.H_H_EL or self.COLLISIONS.H_P_EL or self.COLLISIONS.H2_H_EL: 
            fH_iterate = True

        # Begin Iteration
        fHG = np.zeros((nvr,nvx,nx))
        NHG = np.zeros((nx,self.max_gen+1))
        while True:

            nH_input = np.copy(nH)


            # --- Compute Collision Frequency ---

            # Omega Values (collision frequencies)
            collision_freqs = self._compute_omega_values(fH, nH)
            # Total Collision Frequency
            alpha_c = self._compute_collision_frequency(collision_freqs, gamma_wall)

            # Generate Coefficients
            meq_coeffs = self._compute_mesh_equation_coefficients(alpha_c)


            # --- 0th Generation ---
                            
            # Compute first-flight (0th generation) neutral distribution function
            self._debrief_msg('Computing atomic neutral generation#0', 0)

            fHG[:,self.vx_pos,0] = fH[:,self.vx_pos,0]
            for k in range(nx-1):
                fHG[:,self.vx_pos,k+1] = fHG[:,self.vx_pos,k]*meq_coeffs.A[:,self.vx_pos,k] + meq_coeffs.F[:,self.vx_pos,k]
            for k in range(nx-1,0,-1):
                fHG[:,self.vx_neg,k-1] = fHG[:,self.vx_neg,k]*meq_coeffs.C[:,self.vx_neg,k] + meq_coeffs.G[:,self.vx_neg,k]
                    
            # Compute first-flight neutral density profile
            for k in range(nx):
                NHG[k,0] = np.sum(self.dvr_vol*(fHG[:,:,k] @ self.dvx))

            # Set total atomic neutral distribution function to first flight generation
            fH = fHG.copy()
            nH = NHG[:,0].copy()


            # --- Iterative Generations ---

            fH, nH, fHG, NHG, Beta_CX_sum, m_sums, igen = self._run_generations(fH, nH, fHG, NHG, meq_coeffs, collision_freqs, fH_iterate)
            self.Internal.MH_H_sum = m_sums.H_H

            # Compute H density profile
            for k in range(nx):
                nH[k] = np.sum(self.dvr_vol*(fH[:,:,k] @ self.dvx))


            # --- End Iteration ---

            if fH_iterate:
                # Compute 'seed error': Delta_nHs=(|nHs-nH|)/max(nH)
                # If Delta_nHs is less than 10*truncate then stop iterating fH
                self.Internal.Delta_nHs = np.max(np.abs(nH_input - nH)) / np.max(nH)
                if self.Internal.Delta_nHs <= 10*self.truncate:
                    break
            else:
                break  # No outer iteration needed when fH_iterate=False


        # --- Update Last Generation ---

        # Update Beta_CX_sum using last generation
        Beta_CX = self._compute_beta_cx(fHG)
        Beta_CX_sum += Beta_CX
        
        # Update MH_*_sum using last generation
        m_vals = self._compute_mh_values(fHG, NHG[:,igen])
        m_sums.H_H += m_vals.H_H
        m_sums.H_P += m_vals.H_P
        m_sums.H_H2 += m_vals.H_H2

        return fH, nH, alpha_c, Beta_CX_sum, collision_freqs, m_sums, NHG, igen


    def _run_generations(self, fH, nH, fHG, NHG, meq_coeffs, collision_freqs, fH_iterate):
        '''
        Iterate through and compute generations of collision
        '''

        nvr, nvx, nx = self.nvr, self.nvx, self.nx
        vxp, vxn = self.vx_pos, self.vx_neg

        Beta_CX_sum = np.zeros((nvr,nvx,nx))
        m_sums = CollisionType(np.zeros((nvr,nvx,nx)), np.zeros((nvr,nvx,nx)), np.zeros((nvr,nvx,nx)))

        fH_generations = False
        if fH_iterate or self.COLLISIONS.H_P_CX: 
            fH_generations = True

        igen = 0
        while True:

            if igen >= self.max_gen:
                raise Exception(
                    f'Kinetic_H: failed to converge after {self.max_gen} generations. '
                    f'The {self.max_gen}th generation is still contributing a non-negligible amount '
                    f'to the total neutral density. This means there are neutrals undergoing '
                    f'{self.max_gen} charge exchange or scattering events before ionisation, which '
                    f'is unlikely in typical tokamak conditions and probably indicates a problem '
                    f'with the input profiles.'
                )
            if not fH_generations:
                break
            igen += 1
            self._debrief_msg('Computing atomic neutral generation#'+sval(igen), 0)


            # Compute Beta_CX from previous generation
            Beta_CX = self._compute_beta_cx(fHG)
            # Sum charge exchange source over all generations
            Beta_CX_sum += Beta_CX

            # Elastic collision maxwellians
            m_vals = self._compute_mh_values(fHG, NHG[:,igen-1])
            m_sums.H_H += m_vals.H_H
            m_sums.H_P += m_vals.H_P
            m_sums.H_H2 += m_vals.H_H2

            # Compute next generation molecular distribution
            OmegaM = collision_freqs.H_H*m_vals.H_H + collision_freqs.H_P*m_vals.H_P + collision_freqs.H_H2*m_vals.H_H2
            fHG[:] = 0
            for k in range(nx-1):
                fHG[:,vxp,k+1] = meq_coeffs.A[:,vxp,k]*fHG[:,vxp,k] + meq_coeffs.B[:,vxp,k]*(Beta_CX[:,vxp,k+1] + OmegaM[:,vxp,k+1] + Beta_CX[:,vxp,k] + OmegaM[:,vxp,k])
            for k in range(nx-1, 0, -1):
                fHG[:,vxn,k-1] = meq_coeffs.C[:,vxn,k]*fHG[:,vxn,k] + meq_coeffs.D[:,vxn,k]*(Beta_CX[:,vxn,k-1] + OmegaM[:,vxn,k-1] + Beta_CX[:,vxn,k] + OmegaM[:,vxn,k])
            for k in range(nx):
                NHG[k,igen] = np.sum(self.dvr_vol*(fHG[:,:,k] @ self.dvx))

            # Add result to total neutral distribution function
            fH += fHG
            nH += NHG[:,igen]


            # Compute 'generation error': Delta_nHG=max(NHG(*,igen)/max(nH))
            # and decide if another generation should be computed
            Delta_nHG = np.max(NHG[:,igen] / np.max(nH))
            if (Delta_nHG < self.truncate) or (fH_iterate and (Delta_nHG < 0.003*self.Internal.Delta_nHs)):
                # If fH 'seed' is being iterated, then do another generation until the 'generation error'
                # is less than 0.003 times the 'seed error' or is less than TRUNCATE
                break

        return fH, nH, fHG, NHG, Beta_CX_sum, m_sums, igen
    

    def _compile_results(self, fH, nH, fSH, gamma_wall, alpha_c, Beta_CX_sum, collision_freqs, m_sums):
        '''
        Computes final results of kinetic_h2 procedure, compiles into KH2Results dataclass
        '''

        dvr_vol, dvx = self.dvr_vol, self.dvx
        vx, vr = self.mesh.vx, self.mesh.vr

        # GammaxH - particle flux in x direction
        GammaxH = np.zeros(self.nx)
        for k in range(self.nx):
            GammaxH[k] = self.vth*np.sum(dvr_vol*(fH[:,:,k] @ (vx*dvx)))

        # VxH - x velocity
        VxH = GammaxH / nH
        VxH_vth = VxH / self.vth

        # Magnitude of random velocity at each mesh point 
        # vr2vx2_ran.shape = (self.nvr, self.nvx, self.nx)
        # vr2vx2_ran[i,j,k] = vr[i]**2 + (vx[j] - VxH2[k])**2
        vr2vx2_ran = vr[:,None,None]**2 + (vx[None,:,None] - VxH_vth[None,None,:])**2

        # pH - pressure 
        pH = np.zeros(self.nx)
        pH_coef = (self.mu*CONST.H_MASS)*(self.vth**2) / (3*CONST.Q)
        for k in range(self.nx):
            pH[k] = np.sum(dvr_vol*((vr2vx2_ran[:,:,k]*fH[:,:,k]) @ dvx))
        pH *= pH_coef

        # TH - temperature
        TH = pH/nH

        # piH_xx, piH_yy, piH_zz
        piH_coef = (self.mu*CONST.H_MASS)*(self.vth**2) / CONST.Q
        for k in range(self.nx):
            self.Output.piH_xx[k] = np.sum(dvr_vol*(fH[:,:,k] @ (dvx*(vx - VxH_vth[k])**2)))
            self.Output.piH_yy[k] = np.sum((dvr_vol*(vr**2))*(fH[:,:,k] @ dvx))
        self.Output.piH_xx = piH_coef*self.Output.piH_xx - pH
        self.Output.piH_yy = 0.5*piH_coef*self.Output.piH_yy - pH
        self.Output.piH_zz = np.copy(self.Output.piH_yy)

        # qxH
        qxH = np.zeros(self.nx)
        qxH_coef = 0.5*(self.mu*CONST.H_MASS)*(self.vth**3)
        for k in range(self.nx):
            qxH[k] = np.sum(dvr_vol*((vr2vx2_ran[:,:,k]*fH[:,:,k]) @ (dvx*(vx - VxH_vth[k]))))
        qxH *= qxH_coef

        if self.recomb:
            self.Output.SRecomb = self.vth*self.Internal.ni*self.Internal.Rec
        else:
            self.Output.SRecomb[:] = 0

        QH = np.zeros(self.nx)
        RxH = np.zeros(self.nx)
        NetHSource = np.zeros(self.nx)
        Sion = np.zeros(self.nx)
        SideWallH = np.zeros(self.nx)

        E_coef = 0.5*(self.mu*CONST.H_MASS)*(self.vth**2)
        Rx_coef = (self.mu*CONST.H_MASS)*self.vth
        for k in range(self.nx):
            # C = RHS of Boltzman equation for total fH
            C = self.vth*(self.Internal.Sn[:,:,k] + Beta_CX_sum[:,:,k] - alpha_c[:,:,k]*fH[:,:,k]
                            + collision_freqs.H_P[k]*m_sums.H_P[:,:,k] + collision_freqs.H_H2[k]*m_sums.H_H2[:,:,k] + collision_freqs.H_H[k]*m_sums.H_H[:,:,k])
            
            QH[k] = E_coef*np.sum(dvr_vol*((vr2vx2_ran[:,:,k]*C) @ dvx))
            RxH[k] = Rx_coef*np.sum(dvr_vol*(C @ (dvx*(vx - VxH_vth[k]))))
            NetHSource[k] = np.sum(dvr_vol*(C @ dvx))
            Sion[k] = self.vth*nH[k]*self.Internal.alpha_ion[k]
            self.Output.SourceH[k] = np.sum(dvr_vol*(fSH[:,:,k] @ dvx))
            SideWallH[k] = np.sum(dvr_vol*((gamma_wall[:,:,k]*fH[:,:,k]) @ dvx))

            if self.COLLISIONS.H_P_CX:
                CCX = self.vth*(Beta_CX_sum[:,:,k] - self.Internal.Alpha_CX[:,:,k]*fH[:,:,k])
                self.Output.RxHCX[k] = Rx_coef*np.sum(dvr_vol*(CCX @ (dvx*(vx - VxH_vth[k]))))
                self.Output.EHCX[k] = E_coef*np.sum(dvr_vol*((self.Internal.vr2vx2[:,:,k]*CCX) @ dvx))

            if self.COLLISIONS.H2_H_EL:
                CH_H2 = self.vth*collision_freqs.H_H2[k]*(m_sums.H_H2[:,:,k] - fH[:,:,k])
                self.Output.RxH2_H[k] = Rx_coef*np.sum(dvr_vol*(CH_H2 @ (dvx*(vx - VxH_vth[k]))))
                self.Output.EH2_H[k] = E_coef*np.sum(dvr_vol*((self.Internal.vr2vx2[:,:,k]*CH_H2) @ dvx))

            if self.COLLISIONS.H_P_EL:
                CH_P = self.vth*collision_freqs.H_P[k]*(m_sums.H_P[:,:,k] - fH[:,:,k])
                self.Output.RxP_H[k] = Rx_coef*np.sum(dvr_vol*(CH_P @ (dvx*(vx - VxH_vth[k]))))
                self.Output.EP_H[k] = E_coef*np.sum(dvr_vol*((self.Internal.vr2vx2[:,:,k]*CH_P) @ dvx))

            CW_H = -self.vth*(gamma_wall[:,:,k]*fH[:,:,k])
            self.Output.RxW_H[k] = Rx_coef*np.sum(dvr_vol*(CW_H @ (dvx*(vx - VxH_vth[k]))))
            self.Output.EW_H[k] = E_coef*np.sum(dvr_vol*((self.Internal.vr2vx2[:,:,k]*CW_H) @ dvx))
            
            if self.COLLISIONS.H_H_EL:
                CH_H = self.vth*collision_freqs.H_H[k]*(m_sums.H_H[:,:,k] - fH[:,:,k])
                # vr2vx2_ran2[i,j] = vr[i]**2 + 2*(vx[j] - VxH2_vth[k])**2
                vr2_2vx_ran2 = vr[:,None]**2 - 2*(vx[None,:] - VxH_vth[k])**2
                self.Output.Epara_PerpH_H[k] = -E_coef*np.sum(dvr_vol*((vr2_2vx_ran2*CH_H) @ dvx))

        #	qxH_total
        qxH_total = (0.5*nH*(self.mu*CONST.H_MASS)*VxH*VxH + 2.5*pH*CONST.Q)*VxH + CONST.Q*self.Output.piH_xx*VxH + qxH

        #	QH_total
        QH_total = QH + RxH*VxH + 0.5*(self.mu*CONST.H_MASS)*NetHSource*VxH*VxH

        #	Albedo
        gammax_plus = self.vth*np.sum(dvr_vol*(fH[:,self.vx_pos,0] @ (vx[self.vx_pos]*dvx[self.vx_pos]))) 
        gammax_minus = self.vth*np.sum(dvr_vol*(fH[:,self.vx_neg,0] @ (vx[self.vx_neg]*dvx[self.vx_neg])))
        AlbedoH = 0.0
        if np.abs(gammax_plus) > 0:
            AlbedoH = -gammax_minus/gammax_plus


        results = KHResults(fH, nH, GammaxH, VxH, pH, TH, qxH, qxH_total, NetHSource, Sion, QH, RxH, QH_total, AlbedoH, SideWallH)
        return results




    # ------ Computational Functions ------

    def _compute_omega_values(self, fH, nH):
        '''
        Compute elastic momentum transfer frequencies (omega) using Eqs.(3.12-3.14)
        '''
        
        nvr, nvx, nx = self.nvr, self.nvx, self.nx

        Omega_H_P = np.zeros(nx)
        Omega_H_H2 = np.zeros(nx)
        Omega_H_H = np.zeros(nx)

        # Compute Omega values if nH is non-zero
        if np.any(nH <= 0):
            return CollisionType(Omega_H_H, Omega_H_P, Omega_H_H2)
            
        # Compute VxH
        VxH = np.zeros(nx)
        if self.COLLISIONS.H_P_EL or self.COLLISIONS.H2_H_EL or self.COLLISIONS.H_H_EL:
            for k in range(nx):
                VxH[k] = self.vth*np.sum(self.dvr_vol*(fH[:,:,k] @ (self.mesh.vx*self.dvx))) / nH[k]

        #	Compute Omega_H_P for present fH and Alpha_H_P if H_P elastic collisions are included
        if self.COLLISIONS.H_P_EL:
            self._debrief_msg('Computing Omega_H_P', 1)
            for k in range(nx):
                DeltaVx = (VxH[k] - self.vxi[k]) / self.vth
                MagDeltaVx = np.maximum(abs(DeltaVx), self.DeltaVx_tol)
                DeltaVx = np.sign(DeltaVx)*MagDeltaVx
                Omega_H_P[k] = np.sum(self.dvr_vol*((self.Internal.Alpha_H_P[:,:,k]*fH[:,:,k]) @ self.dvx)) / (nH[k]*DeltaVx)
            Omega_H_P = np.maximum(Omega_H_P, 0)

        #	Compute Omega_H_H2 for present fH and Alpha_H_H2 if H_H2 elastic collisions are included
        if self.COLLISIONS.H2_H_EL:
            self._debrief_msg('Computing Omega_H_H2', 1)
            for k in range(nx):
                DeltaVx = (VxH[k] - self.H2_Moments.VxH2[k]) / self.vth
                MagDeltaVx = np.maximum(abs(DeltaVx), self.DeltaVx_tol)
                DeltaVx = np.sign(DeltaVx)*MagDeltaVx
                Omega_H_H2[k] = np.sum(self.dvr_vol*((self.Internal.Alpha_H_H2[:,:,k]*fH[:,:,k]) @ self.dvx)) / (nH[k]*DeltaVx)
            Omega_H_H2 = np.maximum(Omega_H_H2, 0)

        #	Compute Omega_H_H for present fH if H_H elastic collisions are included
        if self.COLLISIONS.H_H_EL:
            self._debrief_msg('Computing Omega_H_H', 1)

            Wperp_paraH = np.zeros(nx)
            vr2_2vx_ran2 = np.zeros((nvr,nvx))
            if np.sum(self.Internal.MH_H_sum) <= 0:
                for k in range(nx):
                    # vr2vx2_ran2[i,j] = vr[i]**2 + 2*(vx[j] - VxH[k])**2
                    vr2_2vx_ran2 = self.mesh.vr[:,None]**2 - 2*(self.mesh.vx[None,:] - VxH[k])**2
                    Wperp_paraH[k] = np.sum(self.dvr_vol*((vr2_2vx_ran2*fH[:,:,k]) @ self.dvx)) / nH[k]
            else:
                for k in range(nx):
                    M_fH = self.Internal.MH_H_sum[:,:,k] - fH[:,:,k]
                    Wperp_paraH[k] = -np.sum(self.dvr_vol*((self.vr2_2vx2_2D*M_fH) @ self.dvx)) / nH[k]
            
            for k in range(nx):
                Work = fH[:,:,k].reshape((nvr*nvx), order='F')
                Alpha_H_H = (self.Internal.SIG_H_H @ Work).reshape((nvr,nvx), order='F')
                MagWpp = np.maximum(np.abs(Wperp_paraH[k]), self.Wpp_tol)
                Wpp = np.sign(Wperp_paraH[k])*MagWpp
                Omega_H_H[k] = np.sum(self.dvr_vol*((Alpha_H_H*Work.reshape((nvr,nvx), order='F')) @ self.dvx)) / (nH[k]*Wpp)
            Omega_H_H = np.maximum(Omega_H_H, 0)

        return CollisionType(Omega_H_H, Omega_H_P, Omega_H_H2)
    

    def _compute_collision_frequency(self, collision_freqs: CollisionType, gamma_wall):
        '''
        Computes total elastic scattering frequency (Eq. 3.15)
        and total collision frequency (Eq. 3.16)
        '''

        # Total Elastic scattering frequency
        Omega_EL = collision_freqs.H_P + collision_freqs.H_H2 + collision_freqs.H_H

        # Total collision frequency
        alpha_c = np.zeros((self.nvr,self.nvx,self.nx))
        if self.COLLISIONS.H_P_CX:
            for k in range(self.nx):
                alpha_c[:,:,k] = self.Internal.Alpha_CX[:,:,k] + self.Internal.alpha_ion[k] + Omega_EL[k] + gamma_wall[:,:,k]
        else:
            for k in range(self.nx):
                alpha_c[:,:,k] = self.Internal.alpha_ion[k] + Omega_EL[k] + gamma_wall[:,:,k]

        self._test_grid_spacing(alpha_c)

        return alpha_c
    

    def _compute_mesh_equation_coefficients(self, alpha_c):
        '''
        Define parameters Ak, Bk, Ck, Dk, Fk, Gk using Eqs. (3.22), (3.25), (3.30), (3.33)
        '''

        Ak = np.zeros((self.nvr,self.nvx,self.nx))
        Bk = np.zeros((self.nvr,self.nvx,self.nx))
        Ck = np.zeros((self.nvr,self.nvx,self.nx))
        Dk = np.zeros((self.nvr,self.nvx,self.nx))
        Fk = np.zeros((self.nvr,self.nvx,self.nx))
        Gk = np.zeros((self.nvr,self.nvx,self.nx))

        for k in range(0, self.nx-1):
            x_diffs = self.mesh.x[k+1] - self.mesh.x[k]
            for j in self.vx_pos:
                denom = 2*self.mesh.vx[j] + x_diffs*alpha_c[:,j,k+1]
                Ak[:,j,k] = (2*self.mesh.vx[j] - x_diffs*alpha_c[:,j,k]) / denom
                Bk[:,j,k] = x_diffs / denom
                Fk[:,j,k] = x_diffs*(self.Internal.Sn[:,j,k+1]+self.Internal.Sn[:,j,k]) / denom
            for j in self.vx_neg:
                denom = -2*self.mesh.vx[j] + x_diffs*alpha_c[:,j,k]
                Ck[:,j,k+1] = (-2*self.mesh.vx[j] - x_diffs*alpha_c[:,j,k+1]) / denom
                Dk[:,j,k+1] = x_diffs / denom
                Gk[:,j,k+1] = x_diffs*(self.Internal.Sn[:,j,k+1]+self.Internal.Sn[:,j,k]) / denom

        return MeshEqCoefficients(Ak, Bk, Ck, Dk, Fk, Gk)
    

    def _compute_beta_cx(self, fH):
        '''
        Compute charge exchange source (beta_cx) with Eq. (3.11a) or (3.11b)
        '''

        Beta_CX = np.zeros((self.nvr,self.nvx,self.nx))
        if self.COLLISIONS.H_P_CX:
            
            self._debrief_msg('Computing Beta_CX', 1)

            if self.COLLISIONS.SIMPLE_CX:
                # Option (B): Compute charge exchange source with assumption that CX source neutrals have ion distribution function
                for k in range(self.nx):
                    Beta_CX[:,:,k] = self.Internal.fi_hat[:,:,k]*np.sum(self.dvr_vol*((self.Internal.Alpha_CX[:,:,k]*fH[:,:,k]) @ self.dvx))
            else:
                # Option (A): Compute charge exchange source using fH and vr x sigma x v_v at each velocity mesh point
                for k in range(self.nx):
                    Work = fH[:,:,k].reshape((self.nvr*self.nvx), order='F')
                    Beta_CX[:,:,k] = self.Internal.ni[k]*self.Internal.fi_hat[:,:,k]*((self.Internal.SIG_CX @ Work).reshape((self.nvr,self.nvx), order='F'))

        return Beta_CX
    

    def _compute_mh_values(self, fH, nH):
        '''
        Compute collision distributions using Eqs. (3.6)-(3.8)
        '''
        
        MH_H = np.zeros((self.nvr,self.nvx,self.nx))
        MH_P = np.zeros((self.nvr,self.nvx,self.nx))
        MH_H2 = np.zeros((self.nvr,self.nvx,self.nx))
        VxHG = np.zeros(self.nx)
        THG = np.zeros(self.nx)
        if self.COLLISIONS.H_H_EL or self.COLLISIONS.H_P_EL or self.COLLISIONS.H2_H_EL:

            # Compute VxHG, THG
            for k in range(0, self.nx):
                VxHG[k] = self.vth*np.sum(self.dvr_vol*(fH[:,:,k] @ (self.mesh.vx*self.dvx))) / nH[k]
                vr2vx2_ran2 = (self.mesh.vr[:, None]**2 + (self.mesh.vx[None, :] - VxHG[k]/self.vth)**2)
                THG[k] = (self.mu*CONST.H_MASS)*(self.vth**2)*np.sum(self.dvr_vol*((vr2vx2_ran2*fH[:,:,k]) @ self.dvx)) / (3*CONST.Q*nH[k])

            if self.COLLISIONS.H_H_EL:

                self._debrief_msg('Computing MH_H', 1)

                # Compute MH_H
                Maxwell = create_shifted_maxwellian(self.mesh.vr, self.mesh.vx, THG, VxHG, self.mu, 1, self.mesh.Tnorm)
                MH_H = Maxwell*nH

            if self.COLLISIONS.H_P_EL:

                self._debrief_msg('Computing MH_P', 1)

                # Compute MH_P 
                vx_shift = (VxHG + self.vxi) / 2
                Tmaxwell = THG + (2/4)*(self.mesh.Ti - THG + self.mu*CONST.H_MASS*((self.vxi - VxHG)**2) / (6*CONST.Q))
                Maxwell = create_shifted_maxwellian(self.mesh.vr, self.mesh.vx, Tmaxwell, vx_shift, self.mu, 1, self.mesh.Tnorm)
                MH_P = Maxwell*nH

            if self.COLLISIONS.H2_H_EL:

                self._debrief_msg('Computing MH_H2', 1)

                # Compute MH_H2
                vx_shift = (VxHG + 2*self.H2_Moments.VxH2) / 3
                Tmaxwell = THG + (4/9)*(self.H2_Moments.TH2 - THG + 2*self.mu*CONST.H_MASS*((self.H2_Moments.VxH2 - VxHG)**2) / (6*CONST.Q))
                Maxwell = create_shifted_maxwellian(self.mesh.vr, self.mesh.vx, Tmaxwell, vx_shift, self.mu, 1, self.mesh.Tnorm)
                MH_H2 = Maxwell*nH
        
        return CollisionType(MH_H, MH_P, MH_H2)




    # ------ Variable Functions ------

    # --- Initialization ---


    def _init_fhbc_input(self):
        '''
        Computes fH2BC_input, used to scale molecular distribution function (fH) to desired flux
        '''

        self.fHBC_input = np.zeros(self.fHBC.shape)
        self.fHBC_input[:,self.vx_pos] = self.fHBC[:,self.vx_pos]
        gamma_input = 1.0
        if abs(self.GammaxHBC) > 0:
            gamma_input = self.vth*np.sum(self.dvr_vol*(self.fHBC_input @ (self.mesh.vx*self.dvx)))
        ratio = abs(self.GammaxHBC) / gamma_input
        self.fHBC_input = self.fHBC_input*ratio
        if abs(ratio - 1) > 0.01*self.truncate:
            self.fHBC = self.fHBC_input
    

    def _init_static_internals(self):
        '''
        Computes various internal variables based on constant mesh data
        '''

        self._init_grid()
        self._init_protons()
        self._init_sigv()
        self._init_v_v2()
        self._init_sig_cx()
        self._init_sig_h_h()
        self._init_sig_h_h2()
        self._init_sig_h_p()

        return


    def _init_grid(self):
        '''
        Computes internal vr2vx2, vr2vx_vxi2, ErelH_P
        '''

        self._debrief_msg('Computing vr2vx2, vr2vx_vxi2, ErelH_P', 1)

        # Magnitude of total normalized v^2 at each mesh point
        self.Internal.vr2vx2 = np.zeros((self.nvr,self.nvx,self.nx))
        for i in range(self.nvr):
            for k in range(self.nx):
                self.Internal.vr2vx2[i,:,k] = self.mesh.vr[i]**2 + self.mesh.vx**2

        # Magnitude of total normalized (v-vxi)^2 at each mesh point
        self.Internal.vr2vx_vxi2 = np.zeros((self.nvr,self.nvx,self.nx))
        for i in range(self.nvr):
            for k in range(self.nx):
                self.Internal.vr2vx_vxi2[i,:,k] = self.mesh.vr[i]**2 + (self.mesh.vx - self.vxi[k]/self.vth)**2

        # Atomic hydrogen ion energy in local rest frame of plasma at each mesh point
        self.Internal.ErelH_P = (0.5*CONST.H_MASS*self.Internal.vr2vx_vxi2*(self.vth**2)) / CONST.Q
        # sigmav_cx does not handle neutral energies below 0.1 eV or above above 20 keV
        self.Internal.ErelH_P = np.clip(self.Internal.ErelH_P, 0.1, 2.0e4)

        return


    def _init_protons(self):

        # Ti/mu at each mesh point
        self._debrief_msg('Computing Ti/mu at each mesh point', 1)
        self.Internal.Ti_mu = np.zeros((self.nvr,self.nvx,self.nx))
        for k in range(self.nx):
            self.Internal.Ti_mu[:,:,k] = self.mesh.Ti[k] / self.mu

        # Compute Fi_hat
        self._debrief_msg('Computing fi_hat', 1)
        self.Internal.fi_hat = create_shifted_maxwellian(self.mesh.vr, self.mesh.vx, self.mesh.Ti, self.vxi, self.mu, 1, self.mesh.Tnorm)

        return
    

    def _init_sigv(self):

        self._debrief_msg('Computing sigv', 1)

        # Compute sigmav rates for each reaction with option to use rates
        # from CR model of Johnson-Hinnov

        self.Internal.sigv = np.zeros((self.nx,3))

        # Reaction R1:  e + H -> e + H(+) + e   (ionization)
        if self.ion_rate_option == "collrad":
            self.Internal.sigv[:,1] = collrad_sigmav_ion_h0(self.mesh.ne, self.mesh.Te) # from COLLRAD code (DEGAS-2)
        elif self.ion_rate_option == "jh":
            self.Internal.sigv[:,1] = self.jh.jhs_coef(self.mesh.ne, self.mesh.Te, no_null=True) # Johnson-Hinnov, limited Te range
        elif self.ion_rate_option == "adas":
            self.Internal.sigv[:,1] = scd_adas(self.mesh.ne, self.mesh.Te) # ADAS SCD, density-dependent
        else:
            self.Internal.sigv[:,1] = sigmav_ion_h0(self.mesh.Te) # from Janev et al., up to 20keV

        # Reaction R2:  e + H(+) -> H(1s) + hv  (radiative recombination)
        if self.ion_rate_option == "jh":
            self.Internal.sigv[:,2] = self.jh.jhalpha_coef(self.mesh.ne, self.mesh.Te, no_null=True)
        elif self.ion_rate_option == "adas":
            self.Internal.sigv[:,2] = acd_adas(self.mesh.ne, self.mesh.Te) # ADAS ACD, density-dependent
        else:
            self.Internal.sigv[:,2] = sigmav_rec_h1s(self.mesh.Te)

        # H ionization rate (normalized by vth) = reaction 1
        self.Internal.alpha_ion = (self.mesh.ne*self.Internal.sigv[:,1]) / self.vth

        # Recombination rate (normalized by vth) = reaction 2
        self.Internal.Rec = (self.mesh.ne*self.Internal.sigv[:,2]) / self.vth

        return
    

    def _init_v_v2(self):
        '''
        Compute arrays for charge exchange and elastic collision computations
        Computes v_v2, v_v, vr2_vx2, vx_vx, and vr2pidvrdvx
        '''

        self._debrief_msg('Computing v_v2, v_v, vr2_vx2, and vx_vx', 1)

        vr = self.mesh.vr
        vx = self.mesh.vx
        cos_theta = self.cos_theta

        v_starter = vr[:,None,None]**2 + vr[None,:,None]**2 - 2*vr[:,None,None]*vr[None,:,None]*cos_theta[None,None,:]
        vx_diff = (vx[:, None] - vx[None, :])
        vx_diff2 = vx_diff**2

        # v_v2=(v-v_prime)^2 at each double velocity space mesh point, including theta angle
        #   self.Internal.v_v2.shape = (nvr,nvx,nvr,nvx,ntheta))
        #   v_v2[i,j,k,l,m] = v_starter[i,k,m] + (vx[j] - vx[l])**2
        self.Internal.v_v2 = v_starter[:,None,:,None,:] + vx_diff2[None,:,None,:,None]

        # vr2_vx2=0.125* [ vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2 ]
        # at each double velocity space mesh point, including theta angle
        #   self.Internal.vr2_vx2.shape = (nvr,nvx,nvr,nvx,ntheta)
        #   vr2_vx2[i,j,k,l,m] = v_starter[i,k,m] + 2*(vx[j] - vx[l])**2
        self.Internal.vr2_vx2 = v_starter[:,None,:,None,:] - 2*vx_diff2[None,:,None,:,None]

        #	v_v=|v-v_prime| at each double velocity space mesh point, including theta angle
        self.Internal.v_v = np.sqrt(self.Internal.v_v2)

        # vx_vx=(vx-vx_prime) at each double velocity space mesh point
        #   self.Internal.vx_vx.shape = (nvr,nvx,nvr,nvx))
        #   vr_vx[:,j,:,l] = (vx[j] - vx[l])
        self.Internal.vx_vx = np.tile(vx_diff[None,:,None,:], (self.nvr,1,self.nvr,1))

        # Set Vr'2pidVr'*dVx' for each double velocity space mesh point
        #   self.Internal.Vr2pidVrdVx.shape = (nvr,nvx,nvr,nvx)
        #   Vr2pidVrdVx[i,j,k,l] = dvr_vol[k] * dvx[l], repeated over i,j
        self.Internal.Vr2pidVrdVx = np.tile(self.dvr_vol[None,None,:,None]*self.dvx[None,None,None,:], (self.nvr,self.nvx,1,1))

        return
    

    def _init_sig_cx(self):
        '''
        Compute SigmaV_CX from sigma directly for present velocity space grid.
        Charge Exchange Option A
        '''

        self._debrief_msg('Computing SIG_CX', 1)

        #	Compute sigma_cx * v_v at all possible relative velocities
        _Sig = np.zeros((self.nvr*self.nvx*self.nvr*self.nvx, self.ntheta))
        _Sig[:] = (self.Internal.v_v*sigma_cx_h0(self.Internal.v_v2*(0.5*CONST.H_MASS*(self.vth**2)/CONST.Q))).reshape(_Sig.shape, order='F')

        #	Set SIG_CX = vr' x Integral{v_v*sigma_cx} 
        #		over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        self.Internal.SIG_CX = np.zeros((self.nvr*self.nvx, self.nvr*self.nvx))
        self.Internal.SIG_CX[:] = (self.Internal.Vr2pidVrdVx*((_Sig @ self.dtheta).reshape(self.Internal.Vr2pidVrdVx.shape, order='F'))).reshape(self.Internal.SIG_CX.shape, order='F')

        #	SIG_CX is now vr' * sigma_cx(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

        return
    

    def _init_sig_h_h(self):
        '''
        Compute SIG_H_H for present velocity space grid
        '''

        self._debrief_msg('Computing SIG_H_H', 1)

        #	Compute sigma_H_H * vr2_vx2 * v_v at all possible relative velocities
        _Sig = np.zeros((self.nvr*self.nvx*self.nvr*self.nvx,self.ntheta))
        _Sig[:] = (self.Internal.vr2_vx2*self.Internal.v_v*sigma_el_h_h(self.Internal.v_v2*(0.5*CONST.H_MASS*self.mu*(self.vth**2)/CONST.Q), vis=True) / 8).reshape(_Sig.shape, order='F')

        #	Note: For viscosity, the cross section for D -> D is the same function of center of mass energy as H -> H.

        #	Set SIG_H_H = vr' x Integral{vr2_vx2*v_v*sigma_H_H} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        self.Internal.SIG_H_H = np.zeros((self.nvr*self.nvx,self.nvr*self.nvx))
        self.Internal.SIG_H_H[:] = (self.Internal.Vr2pidVrdVx*(_Sig @ self.dtheta).reshape(self.Internal.Vr2pidVrdVx.shape, order='F')).reshape(self.Internal.SIG_H_H.shape, order='F')
        #	SIG_H_H is now vr' * sigma_H_H(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

        return
    

    def _init_sig_h_h2(self):
        '''
        Compute SIG_H_H2 for present velocity space grid
        '''

        self._debrief_msg('Computing SIG_H_H2', 1)

        # Compute sigma_H_H2 * v_v at all possible relative velocities
        _Sig = np.zeros((self.nvr*self.nvx*self.nvr*self.nvx,self.ntheta))
        _Sig[:] = (self.Internal.v_v*sigma_el_h_hh(self.Internal.v_v2*(0.5*CONST.H_MASS*(self.vth**2)/CONST.Q))).reshape(_Sig.shape, order='F')

        # Note: using H energy here for cross-sections tabulated as H->H2

        # Set SIG_H_H2 = vr' x vx_vx x Integral{v_v*sigma_H_H2} over theta=0,
        #   2pi times differential velocity space element Vr'2pidVr'*dVx'
        self.Internal.SIG_H_H2 = np.zeros((self.nvr*self.nvx,self.nvr*self.nvx))
        self.Internal.SIG_H_H2[:] = (self.Internal.Vr2pidVrdVx*self.Internal.vx_vx*(_Sig @ self.dtheta).reshape(self.Internal.Vr2pidVrdVx.shape, order='F')).reshape(self.Internal.SIG_H_H2.shape, order='F')

        # SIG_H_H2 is now vr' *vx_vx * sigma_H_H2(v_v) * v_v 
        #   (intergated over theta) for all possible ([vr,vx],[vr',vx'])
        
        return
    
    
    def _init_sig_h_p(self):
        '''
        Compute SIG_H_P for present velocity space grid
        '''

        self._debrief_msg('Computing SIG_H_P', 1)

        # Compute sigma_H_P * v_v at all possible relative velocities
        _Sig = np.zeros((self.nvr*self.nvx*self.nvr*self.nvx,self.ntheta))
        _Sig[:] = (self.Internal.v_v*sigma_el_p_h(self.Internal.v_v2*(0.5*CONST.H_MASS*(self.vth**2)/CONST.Q))).reshape(_Sig.shape, order='F')

        # Set SIG_H_P = vr' x vx_vx x Integral{v_v*sigma_H_P} over theta=0,
        #   2pi times differential velocity space element Vr'2pidVr'*dVx'
        self.Internal.SIG_H_P = np.zeros((self.nvr*self.nvx,self.nvr*self.nvx))
        self.Internal.SIG_H_P[:] = (self.Internal.Vr2pidVrdVx*self.Internal.vx_vx*(_Sig @ self.dtheta).reshape(self.Internal.Vr2pidVrdVx.shape, order='F')).reshape(self.Internal.SIG_H_P.shape, order='F')

        # SIG_H_P is now vr' *vx_vx * sigma_H_P(v_v) * v_v (intergated over theta) 
        #   for all possible ([vr,vx],[vr',vx'])

        return
    


    # --- Procedural ---

    def _compute_dynamic_internals(self, fH, fH2, nHP, THP, fSH):
        '''
        Determines which internal variables need to be recomputed based on changes in input across iterations
        '''

        New_Molecular_Ions = True
        if (self.Input.nHP_s is not None) and np.array_equal(self.Input.nHP_s, nHP) and np.array_equal(self.Input.THP_s, THP):
            New_Molecular_Ions = False

        New_fH2 = True
        if (self.Input.fH2_s is not None) and np.array_equal(self.Input.fH2_s, fH2):
            New_fH2 = False

        New_H_Seed = True
        if (self.Input.fH_s is not None) and  np.array_equal(self.Input.fH_s, fH):
            New_H_Seed = False

        # Reset H2 Moments
        self.H2_Moments.nH2 = np.zeros(self.nx)
        self.H2_Moments.VxH2 = np.zeros(self.nx)
        self.H2_Moments.TH2 = np.full(self.nx, 1.0)

        if New_H_Seed:
            self.Internal.MH_H_sum = np.zeros((self.nvr,self.nvx,self.nx))
            self.Internal.Delta_nHs = 1
        if New_fH2 and (np.sum(fH2) > 0.0):
            self._compute_fh2_moments(fH2)
        if New_Molecular_Ions:
            self._compute_ni(nHP)
        self._compute_sn(fSH)

        # Set up arrays for charge exchange and elastic collision computations, if needed
        if ((self.Internal.Alpha_CX is None) | New_Molecular_Ions) and self.COLLISIONS.H_P_CX:
            self._compute_alpha_cx()
        if ((self.Internal.Alpha_H_H2 is None) | New_fH2) and self.COLLISIONS.H2_H_EL:
            self._compute_alpha_h_h2(fH2)
        if ((self.Internal.Alpha_H_P is None) | New_Molecular_Ions) and self.COLLISIONS.H_P_EL:
            self._compute_alpha_h_p()

        return


    def _compute_fh2_moments(self, fH2):
        '''
        Computes moments from molecular hydrogen distribution functions
        '''

        self._debrief_msg('Computing vx and T moments of fH2', 1)
        
        # Compute x flow velocity and temperature of molecular species
        for k in range(self.nx):
            self.H2_Moments.nH2[k] = np.sum(self.dvr_vol*(fH2[:,:,k] @ self.dvx))
            if self.H2_Moments.nH2[k] <= 0:
                continue
            self.H2_Moments.VxH2[k] = self.vth*np.sum(self.dvr_vol*(fH2[:,:,k] @ (self.mesh.vx*self.dvx))) / self.H2_Moments.nH2[k]
            vr2vx2_ran2 = self.mesh.vr[:, None]**2 + (self.mesh.vx[None, :] - self.H2_Moments.VxH2[k]/self.vth)**2
            self.H2_Moments.TH2[k] = (2*self.mu*CONST.H_MASS)*(self.vth**2)*np.sum(self.dvr_vol*((vr2vx2_ran2*fH2[:,:,k]) @ self.dvx)) / (3*CONST.Q*self.H2_Moments.nH2[k])


    def _compute_ni(self, nHP):
        '''
        Computes ni profile
        '''

        self._debrief_msg('Computing ni profile', 1)

        self.Internal.ni = self.mesh.ne
        if self.ni_correct:
            self.Internal.ni = self.mesh.ne - nHP
        self.Internal.ni = np.maximum(self.Internal.ni, 0.01*self.mesh.ne)


    def _compute_sn(self, fSH):
        '''
        Compute Total Atomic Hydrogen Source using Eq. (3.18)
        '''

        self.Internal.Sn = np.zeros((self.nvr,self.nvx,self.nx))

        # Add Recombination (optionally) and User-Supplied Hydrogen Source (velocity space distribution)
        for k in range(self.nx):
            self.Internal.Sn[:,:,k] = fSH[:,:,k]/self.vth
            if self.recomb:
                self.Internal.Sn[:,:,k] = self.Internal.Sn[:,:,k] + self.Internal.fi_hat[:,:,k]*self.Internal.ni[k]*self.Internal.Rec[k]

    
    def _compute_alpha_cx(self):
        '''
        Compute charge exchange collision frequency (alpha_cx) using Eq.(2.10a) or (2.10b)
        '''
 
        self._debrief_msg('Computing Alpha_CX', 1)

        if self.COLLISIONS.SIMPLE_CX:
            # Option (B): Use maxwellian weighted <sigma v>

            # Charge Exchange sink rate
            self.Internal.Alpha_CX = sigmav_cx_h0(self.Internal.Ti_mu, self.Internal.ErelH_P) / self.vth
            for k in range(self.nx):
                self.Internal.Alpha_CX[:,:,k] = self.Internal.Alpha_CX[:,:,k]*self.Internal.ni[k]

        else:
            # Option (A): Compute SigmaV_CX from sigma directly via SIG_CX
            self.Internal.Alpha_CX = np.zeros((self.nvr,self.nvx,self.nx))
            for k in range(self.nx):
                Work = (self.Internal.fi_hat[:,:,k]*self.Internal.ni[k]).reshape((self.nvr*self.nvx), order='F')
                self.Internal.Alpha_CX[:,:,k] = (self.Internal.SIG_CX @ Work).reshape(self.Internal.Alpha_CX[:,:,k].shape, order='F')
            
            if self.Do_Alpha_CX_Test: # NOTE Not tested/implemented
                Alpha_CX_Test = sigmav_cx_h0(self.Internal.Ti_mu, self.Internal.ErelH_P) / self.vth
                for k in range(self.nx):
                    Alpha_CX_Test[:,:,k] = Alpha_CX_Test[:,:,k]*self.Internal.ni[k]
                print('Compare alpha_cx and alpha_cx_test')

        return


    def _compute_alpha_h_h2(self, fH2):
        '''
        Compute H:H2 Elastic momentum transfer frequency using Eq.(3.13)
        '''
        
        self._debrief_msg('Computing Alpha_H_H2', 1)

        # Compute Alpha_H_H2 for inputted fH, if it is needed and has not
        #   already been computed with the present input parameters
        
        self.Internal.Alpha_H_H2 = np.zeros((self.nvr,self.nvx,self.nx))
        for k in range(self.nx):
            Work = fH2[:,:,k].reshape((self.nvr*self.nvx), order='F')
            self.Internal.Alpha_H_H2[:,:,k] = (self.Internal.SIG_H_H2 @ Work).reshape(self.Internal.Alpha_H_H2[:,:,k].shape, order='F')

        return


    def _compute_alpha_h_p(self):
        '''
        Compute H:P Elastic momentum transfer frequency using Eq.(3.12)
        '''

        self._debrief_msg('Computing Alpha_H_P', 1)

        # Compute Alpha_H_P for present Ti and ni 
        #   if it is needed and has not already been computed with the present parameters

        self.Internal.Alpha_H_P = np.zeros((self.nvr,self.nvx,self.nx))
        for k in range(self.nx):
            Work = (self.Internal.fi_hat[:,:,k]*self.Internal.ni[k]).reshape((self.nvr*self.nvx), order='F')
            self.Internal.Alpha_H_P[:,:,k] = (self.Internal.SIG_H_P @ Work).reshape(self.Internal.Alpha_H_P[:,:,k].shape, order='F')    

        return
    



    # ------ Error Computation ------

    def _compute_vbar_error(self):
        if self.debrief>1:
            print(self.prompt+'Computing Vbar_Error')

        #	Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti[x]/Tnorm)/sqrt(pi)

        vx_shift = np.zeros(self.nx)
        Maxwell = create_shifted_maxwellian(self.mesh.vr, self.mesh.vx, self.mesh.Ti, vx_shift, self.mu, 1, self.mesh.Tnorm)
        
        vbar_test = np.zeros((self.nvr,self.nvx,self.ntheta))
        self.Errors.vbar_error = np.zeros(self.nx)
        for m in range(self.ntheta):
            vbar_test[:,:,m] = self.Internal.vr2vx2[:,:,0]
        _vbar_test = np.zeros((self.nvr*self.nvx,self.ntheta))
        _vbar_test[:] = (self.vth*np.sqrt(vbar_test)).reshape(_vbar_test.shape, order='F')
        vbar_test = np.zeros((self.nvr,self.nvx))
        vbar_test[:] = (_vbar_test @ self.dtheta).reshape(vbar_test.shape, order='F')
        for k in range(self.nx):
            vbar = np.sum(self.dvr_vol*((vbar_test*Maxwell[:,:,k]) @ self.dvx))
            vbar_exact = 2*self.vth*np.sqrt(self.mesh.Ti[k]/self.mesh.Tnorm) / np.sqrt(np.pi)
            self.Errors.vbar_error[k] = abs(vbar-vbar_exact) / vbar_exact
        if self.debrief > 0:
            print(self.prompt+'Maximum Vbar error = ', sval(max(self.Errors.vbar_error)))


    def _compute_final_errors(self, results, Beta_CX_sum, m_sums, alpha_c, collision_freqs):
        if self.debrief > 1:
            print(self.prompt+'Computing Collision Operator, Mesh, and Moment Normalized Errors')

        #	Compute Mesh Errors

        self.Errors.mesh_error = np.zeros((self.nvr,self.nvx,self.nx))
        max_mesh_error = 0.0
        min_mesh_error = 0.0
        mtest = 5
        self.Errors.moment_error = np.zeros((self.nx,mtest))
        max_moment_error = np.zeros(mtest)
        self.Errors.C_error = np.zeros(self.nx)
        self.Errors.CX_error = np.zeros(self.nx)
        self.Errors.H_H_error = np.zeros((self.nx, 3))
        H_H2_error = np.zeros((self.nx, 3))
        H_P_error = np.zeros((self.nx, 3))
        max_H_H_error = np.zeros(3)
        max_H_H2_error = np.zeros(3)
        max_H_P_error = np.zeros(3)

        NetHSource2 = self.Output.SourceH + self.Output.SRecomb - results.Sion - results.SideWallH
        for k in range(self.nx):
            self.Errors.C_error[k] = abs(results.NetHSource[k] - NetHSource2[k]) / max(abs(np.array([results.NetHSource[k], NetHSource2[k]])))

        #	Test conservation of particles for charge exchange operator
        if self.COLLISIONS.H_P_CX:
            for k in range(self.nx):
                CX_A = np.sum(self.dvr_vol*((self.Internal.Alpha_CX[:,:,k]*results.fH[:,:,k]) @ self.dvx))
                CX_B = np.sum(self.dvr_vol*(Beta_CX_sum[:,:,k] @ self.dvx))
                self.Errors.CX_error[k] = np.abs(CX_A - CX_B) / np.max(np.abs(np.array([CX_A, CX_B])))

        #	Test conservation of particles, x momentum, and total energy of elastic collision operators
        for m in range(0, 3):
            for k in range(0, self.nx):
                if m < 2:
                    TfH = np.sum(self.dvr_vol*(results.fH[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                else:
                    TfH = np.sum(self.dvr_vol*((self.Internal.vr2vx2[:,:,k]*results.fH[:,:,k]) @ self.dvx))

                if self.COLLISIONS.H_H_EL:
                    if m < 2:
                        TH_H = np.sum(self.dvr_vol*(m_sums.H_H[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                    else:
                        TH_H = np.sum(self.dvr_vol*((self.Internal.vr2vx2[:,:,k]*m_sums.H_H[:,:,k]) @ self.dvx))
                    self.Errors.H_H_error[k,m] = np.abs(TfH - TH_H) / np.max(np.abs(np.array([TfH, TH_H])))
                
                if self.COLLISIONS.H2_H_EL:
                    if m < 2:
                        TH_H2 = np.sum(self.dvr_vol*(m_sums.H_H2[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                    else:
                        TH_H2 = np.sum(self.dvr_vol*((self.Internal.vr2vx2[:,:,k]*m_sums.H_H2[:,:,k]) @ self.dvx))
                    H_H2_error[k,m] = np.abs(TfH - TH_H2) / np.max(np.abs(np.array([TfH, TH_H2])))

                if self.COLLISIONS.H_P_EL:
                    if m < 2:
                        TH_P = np.sum(self.dvr_vol*(m_sums.H_P[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                    else:
                        TH_P = np.sum(self.dvr_vol*((self.Internal.vr2vx2[:,:,k]*m_sums.H_P[:,:,k]) @ self.dvx))
                    H_P_error[k,m] = np.abs(TfH - TH_P) / np.max(np.abs(np.array([TfH, TH_P])))

            max_H_H_error[m] = np.max(self.Errors.H_H_error[:,m])
            max_H_H2_error[m] = np.max(H_H2_error[:,m])
            max_H_P_error[m] = np.max(H_P_error[:,m])

        if self.CI_Test:
            #	Compute Momentum transfer rate via full collision integrals for charge exchange and 
            #		mixed elastic scattering.
            #		Then compute error between this and actual momentum transfer 
            #		resulting from CX and BKG (elastic) models.

            if self.COLLISIONS.H_P_CX: # P -> H charge exchange momentum transfer via full collision integral
                print(self.prompt, 'Computing P -> H2 Charge Exchange Momentum Transfer')
                _Sig = np.zeros((self.nvr*self.nvx*self.nvr*self.nvx,self.ntheta))
                _Sig[:] = (self.Internal.v_v*sigma_cx_h0(self.Internal.v_v2*(0.5*CONST.H_MASS*(self.vth**2) / CONST.Q))).reshape(_Sig.shape, order='F')
                SIG_VX_CX = np.zeros((self.nvr*self.nvx,self.nvr*self.nvx))
                SIG_VX_CX[:] = (self.Internal.Vr2pidVrdVx*self.Internal.vx_vx*((_Sig @ self.dtheta).reshape(self.Internal.vx_vx.shape, order='F'))).reshape(SIG_VX_CX.shape, order='F')
                alpha_vx_cx = np.zeros((self.nvr,self.nvx,self.nx))

                for k in range(0, self.nx):
                    Work = (self.Internal.fi_hat[:,:,k]*self.Internal.ni[k]).reshape((self.nvr*self.nvx), order='F')
                    alpha_vx_cx[:,:,k] = (SIG_VX_CX @ Work).reshape(alpha_vx_cx[:,:,k].shape, order='F')

                RxCI_CX = np.zeros(self.nx)
                for k in range(0, self.nx):
                    RxCI_CX[k] = -(self.mu*CONST.H_MASS)*(self.vth**2)*np.sum(self.dvr_vol*((alpha_vx_cx[:,:,k]*results.fH[:,:,k]) @ self.dvx))

                norm = np.max(np.abs(np.array([self.Output.RxHCX, RxCI_CX])))
                CI_CX_error = np.zeros(self.nx)
                for k in range(0, self.nx):
                    CI_CX_error[k] = np.abs(self.Output.RxHCX[k] - RxCI_CX[k]) / norm

                print(self.prompt,'Maximum normalized momentum transfer error in CX collision operator: ', sval(np.max(CI_CX_error)))

            if self.COLLISIONS.H_P_EL: # P -> H momentum transfer via full collision integral
                RxCI_P_H = np.zeros(self.nx)
                for k in range(0, self.nx):
                    RxCI_P_H[k] = -(1/2)*(self.mu*CONST.H_MASS)*(self.vth**2)*np.sum(self.dvr_vol*((self.Internal.Alpha_H_P[:,:,k]*results.fH[:,:,k]) @ self.dvx))

                norm = np.max(np.abs(np.array([self.Output.RxP_H, RxCI_P_H])))
                CI_P_H_error = np.zeros(self.nx)
                for k in range(0, self.nx):
                    CI_P_H_error[k] = np.abs(self.Output.RxP_H[k] - RxCI_P_H[k]) / norm 

                print(self.prompt, 'Maximum normalized momentum transfer error in P -> H elastic BKG collision operator: ', sval(np.max(CI_P_H_error)))

            if self.COLLISIONS.H2_H_EL: # H2 -> H momentum transfer via full collision integral
                RxCI_H2_H = np.zeros(self.nx)
                for k in range(0, self.nx):
                    RxCI_H2_H[k] = -(2/3)*(self.mu*CONST.H_MASS)*(self.vth**2)*np.sum(self.dvr_vol*((self.Internal.Alpha_H_H2[:,:,k]*results.fH[:,:,k]) @ self.dvx))
                
                norm = np.max(np.abs(np.array([self.Output.RxH2_H, RxCI_H2_H])))
                CI_H2_H_error = np.zeros(self.nx)
                for k in range(0, self.nx):
                    CI_H2_H_error[k] = np.abs(self.Output.RxH2_H[k] - RxCI_H2_H[k])/norm
                
                print(self.prompt, 'Maximum normalized momentum transfer error in H2 -> H elastic BKG collision operator: ', sval(np.max(CI_H2_H_error)))

            if self.COLLISIONS.H_H_EL: # H -> H perp/parallel energy transfer via full collision integral
                Epara_Perp_CI = np.zeros(self.nx)
                for k in range(0, self.nx):
                    Work = results.fH[:,:,k].reshape((self.nvr*self.nvx), order='F')
                    Alpha_H_H = (self.Internal.SIG_H_H @ Work).reshape((self.nvr,self.nvx), order='F')
                    Epara_Perp_CI[k] = 0.5*(self.mu*CONST.H_MASS)*(self.vth**3)*np.sum(self.dvr_vol*((Alpha_H_H*results.fH[:,:,k]) @ self.dvx)) 
                
                norm = np.max(np.abs(np.array([self.Output.Epara_PerpH_H, Epara_Perp_CI])))
                CI_H_H_error = np.zeros(self.nx)
                for k in range(0, self.nx):
                    CI_H_H_error[k] = np.abs(self.Output.Epara_PerpH_H[k] - Epara_Perp_CI[k]) / norm 
                
                print(self.prompt, 'Maximum normalized perp/parallel energy transfer error in H -> H elastic BKG collision operator: ', sval(np.max(CI_H_H_error)))

        #	Mesh Point Error based on fH satisfying Boltzmann equation

        T1 = np.zeros((self.nvr,self.nvx,self.nx))
        T2 = np.zeros((self.nvr,self.nvx,self.nx))
        T3 = np.zeros((self.nvr,self.nvx,self.nx))
        T4 = np.zeros((self.nvr,self.nvx,self.nx))
        T5 = np.zeros((self.nvr,self.nvx,self.nx))
        for k in range(0, self.nx-1):
            for j in range(0, self.nvx):
                T1[:,j,k] = 2*self.mesh.vx[j]*(results.fH[:,j,k+1] - results.fH[:,j,k]) / (self.mesh.x[k+1] - self.mesh.x[k]) 
            T2[:,:,k] = (self.Internal.Sn[:,:,k+1] + self.Internal.Sn[:,:,k])
            T3[:,:,k] = Beta_CX_sum[:,:,k+1] + Beta_CX_sum[:,:,k]
            T4[:,:,k] = alpha_c[:,:,k+1]*results.fH[:,:,k+1] + alpha_c[:,:,k]*results.fH[:,:,k]
            T5[:,:,k] = collision_freqs.H_P[k+1]*m_sums.H_P[:,:,k+1] + collision_freqs.H_H2[k+1]*m_sums.H_H2[:,:,k+1] + collision_freqs.H_H[k+1]*m_sums.H_H[:,:,k+1] + \
                        collision_freqs.H_P[k]*m_sums.H_P[:,:,k] + collision_freqs.H_H2[k]*m_sums.H_H2[:,:,k] + collision_freqs.H_H[k]*m_sums.H_H[:,:,k]
            self.Errors.mesh_error[:,:,k] = np.abs(T1[:,:,k] - T2[:,:,k] - T3[:,:,k] + T4[:,:,k] - T5[:,:,k]) / \
                                np.max(np.abs(np.array([T1[:,:,k], T2[:,:,k], T3[:,:,k], T4[:,:,k], T5[:,:,k]])))
        ave_mesh_error = np.sum(self.Errors.mesh_error) / np.size(self.Errors.mesh_error)
        max_mesh_error = np.max(self.Errors.mesh_error)
        min_mesh_error = np.min(self.Errors.mesh_error[:,:,0:self.nx-1])

        #	Moment Error
        for m in range(0, mtest):
            for k in range(0, self.nx-1):
                MT1 = np.sum(self.dvr_vol*(T1[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                MT2 = np.sum(self.dvr_vol*(T2[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                MT3 = np.sum(self.dvr_vol*(T3[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                MT4 = np.sum(self.dvr_vol*(T4[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                MT5 = np.sum(self.dvr_vol*(T5[:,:,k] @ (self.dvx*(self.mesh.vx**m))))
                #NOTE This is correct for the original code, but is it correct mathematically?
                self.Errors.moment_error[k,m] = np.abs(MT1 - MT2 - MT3 + MT4 - MT5) / np.max(np.abs(np.array([MT1, MT2, MT3, MT4, MT5])))
            max_moment_error[m] = np.max(self.Errors.moment_error[:,m])

        #	Compute error in qxH_total

        #		qxH_total2 total neutral heat flux profile (watts m^-2)
        #			This is the total heat flux transported by the neutrals
        #			computed in a different way from:

        #			qxH_total2(k)=Vth**3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#(Vx*dVx)))*0.5*(mu*mH)

        #			This should agree with qxH_total if the definitions of nH, pH, piH_xx,
        #			TH, VxH, and qxH are coded correctly.
        qxH_total2 = np.zeros(self.nx)
        for k in range(0, self.nx):
            qxH_total2[k] = 0.5*(self.mu*CONST.H_MASS)*(self.vth**3)*np.sum(self.dvr_vol*((self.Internal.vr2vx2[:,:,k]*results.fH[:,:,k]) @ (self.mesh.vx*self.dvx)))
        self.Errors.qxH_total_error = np.abs(results.qxH_total - qxH_total2) / np.max(np.abs(np.array([results.qxH_total, qxH_total2])))

        #	Compute error in QH_total
        Q1 = np.zeros(self.nx)
        Q2 = np.zeros(self.nx)
        self.Errors.QH_total_error = np.zeros(self.nx)
        for k in range(0, self.nx-1):
            Q1[k] = (results.qxH_total[k+1] - results.qxH_total[k]) / (self.mesh.x[k+1] - self.mesh.x[k])
            Q2[k] = 0.5*(results.QH_total[k+1] + results.QH_total[k])
        self.Errors.QH_total_error = np.abs(Q1 - Q2) / np.max(np.abs(np.array([Q1, Q2])))

        if self.debrief > 0:
            print(self.prompt+'Maximum particle convervation error of total collision operator: '+sval(max(self.Errors.C_error)))
            print(self.prompt+'Maximum H_P_CX  particle convervation error: '+sval(max(self.Errors.CX_error)))
            print(self.prompt+'Maximum H_H_EL  particle conservation error: '+sval(max_H_H_error[0]))
            print(self.prompt+'Maximum H_H_EL  x-momentum conservation error: '+sval(max_H_H_error[1]))
            print(self.prompt+'Maximum H_H_EL  total energy conservation error: '+sval(max_H_H_error[2]))
            print(self.prompt+'Maximum H_H2_EL particle conservation error: '+sval(max_H_H2_error[0]))
            print(self.prompt+'Maximum H_P_EL  particle conservation error: '+sval(max_H_P_error[0]))
            print(self.prompt+'Average mesh_error = '+str(ave_mesh_error))
            print(self.prompt+'Maximum mesh_error = '+str(max_mesh_error))
            for m in range(5):
                print(self.prompt+'Maximum fH vx^'+sval(m)+' moment error: '+sval(max_moment_error[m]))
            print(self.prompt+'Maximum qxH_total error = '+str(max(self.Errors.qxH_total_error)))
            print(self.prompt+'Maximum QH_total error = '+str(max(self.Errors.QH_total_error)))
            if self.debug > 0:
                input()




    # ------ Testing Functions ------

    def _test_init_parameters(self):
        '''
        Performs compatibility tests for passed parameters when initializing class
        '''

        dx = self.mesh.x - np.roll(self.mesh.x, 1)
        dx = dx[1:]
        notpos = np.argwhere(dx <= 0)
        if notpos.size > 0:
            raise Exception(self.prompt + 'x must be increasing with index!')
        if (self.nvx % 2) != 0:
            raise Exception(self.prompt + 'Number of elements in vx must be even!') 
        if self.mesh.Ti.size != self.nx:
            raise Exception(self.prompt + 'Number of elements in Ti and x do not agree!')
        if self.vxi is None:
            self.vxi = np.zeros(self.nx)
        if self.vxi.size != self.nx:
            raise Exception(self.prompt + 'Number of elements in vxi and x do not agree!')
        if self.mesh.Te.size != self.nx:
            raise Exception(self.prompt + 'Number of elements in Te and x do not agree!')
        if self.mesh.ne.size != self.nx:
            raise Exception(self.prompt + 'Number of elements in n and x do not agree!')
        if self.GammaxHBC is None:
            raise Exception(self.prompt + 'GammaxHBC is not defined!')
        if self.mesh.PipeDia is None:
            self.mesh.PipeDia = np.zeros(self.nx)
        if self.mesh.PipeDia.size != self.nx:
            raise Exception(self.prompt + 'Number of elements in PipeDia and x do not agree!')
        if len(self.fHBC[:,0]) != self.nvr:
            raise Exception(self.prompt + 'Number of elements in fHBC[:,0] and vr do not agree!')
        if len(self.fHBC[0,:]) != self.nvx:
            raise Exception(self.prompt + 'Number of elements in fHBC[0,:] and vx do not agree!')
        if np.sum(abs(self.mesh.vr)) <= 0:
            raise Exception(self.prompt+'vr is all 0!')
        ii = np.argwhere(self.mesh.vr <= 0)
        if ii.size > 0:
            raise Exception(self.prompt+'vr contains zero or negative element(s)!')
        if np.sum(abs(self.mesh.vx)) <= 0:
            raise Exception(self.prompt+'vx is all 0!')
        if np.sum(self.mesh.x) <= 0:
            raise Exception(self.prompt+'Total(x) is less than or equal to 0!')
        if self.mu is None:
            raise Exception(self.prompt+'mu is not defined!')
        if self.mu not in [1,2]:
            raise Exception(self.prompt+'mu must be 1 or 2!')
        
        if np.size(self.vx_neg) < 1:
            raise Exception(self.prompt+'vx contains no negative elements!')
        if np.size(self.vx_pos) < 1:
            raise Exception(self.prompt+'vx contains no positive elements!')
        if np.size(self.vx_zero) > 0:
            raise Exception(self.prompt+'vx contains one or more zero elements!')
        diff = np.nonzero(self.mesh.vx[self.vx_pos] != -np.flipud(self.mesh.vx[self.vx_neg]))[0]
        if diff.size > 0:
            raise Exception(self.prompt + " vx array elements are not symmetric about zero!")
        
        if np.sum(self.fHBC_input) <= 0.0 and abs(self.GammaxHBC) > 0:
            raise Exception(self.prompt+'Values for fHBC[:,:] with vx > 0 are all zero!')
        
        return
        
    
    def _test_input_parameters(self, fH2, fSH, nHP, THP, fH):
        '''
        Performs compatibility tests for passed parameters when calling run_generation
        '''

        if fH2[:,0,0].size != self.nvr:
            raise Exception(self.prompt+'Number of elements in fH2[:,0,0] and vr do not agree!')
        if fH2[0,:,0].size != self.nvx:
            raise Exception(self.prompt+'Number of elements in fH2[0,:,0] and vx do not agree!')
        if fH2[0,0,:].size != self.nx:
            raise Exception(self.prompt+'Number of elements in fH2[0,0,:] and x do not agree!')
        if fSH[:,0,0].size != self.nvr:
            raise Exception(self.prompt+'Number of elements in fSH[:,0,0] and vr do not agree!')
        if fSH[0,:,0].size != self.nvx:
            raise Exception(self.prompt+'Number of elements in fSH[0,:,0] and vx do not agree!')
        if fSH[0,0,:].size != self.nx:
            raise Exception(self.prompt+'Number of elements in fSH[0,0,:] and x do not agree!')
        if nHP.size != self.nx:
            raise Exception(self.prompt+'Number of elements in nHP and x do not agree!')
        if THP.size != self.nx:
            raise Exception(self.prompt+'Number of elements in nHP and x do not agree!')
        if fH[:,0,0].size != self.nvr:
            raise Exception(self.prompt+'Number of elements in fH[:,0,0] and vr do not agree!')
        if fH[0,:,0].size != self.nvx:
            raise Exception(self.prompt+'Number of elements in fH[0,:,0] and vx do not agree!')
        if fH[0,0,:].size != self.nx:
            raise Exception(self.prompt+'Number of elements in fH[0,0,:] and x do not agree!')
        
        return
    

    def _test_grid_spacing(self, alpha_c):
        #	Test x grid spacing based on Eq.(27) in notes
        if self.debrief > 1:
            print(self.prompt+'Testing x grid spacing')
        self.Errors.Max_dx = np.full(self.nx, 1e32)
        for k in range(self.nx):
            for j in self.vx_pos:
                with np.errstate(divide='ignore', invalid='ignore'):
                    local_dx = 2.0 * self.mesh.vx[j] / alpha_c[:, j, k]
                    local_dx = local_dx[np.isfinite(local_dx) & (local_dx > 0)]
                    if local_dx.size > 0:
                        self.Errors.Max_dx[k] = np.minimum(self.Errors.Max_dx[k], np.min(local_dx))

        dx = np.roll(self.mesh.x, -1) - self.mesh.x
        Max_dxL = self.Errors.Max_dx[0:self.nx-1]
        Max_dxR = self.Errors.Max_dx[1:self.nx]
        self.Errors.Max_dx = np.minimum(Max_dxL, Max_dxR)
        ilarge = np.argwhere(self.Errors.Max_dx < dx[0:self.nx-1])

        if ilarge.size > 0:
            print(self.prompt+'x mesh spacing is too large!') #NOTE Check Formatting
            out = ""
            jj = 0

            #	Not sure the output is formatted correctly

            print(' \t    x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)')
            for ii in range(ilarge.size):
                jj += 1
                out += ((str(ilarge[ii])+' \t')[:8]+(str(self.mesh.x[ilarge[ii]+1]-self.mesh.x[ilarge[ii]])+'        ')[:6]+'        '+str(self.Errors.Max_dx[ilarge[ii]])[:4]+'\t')
                if jj>4:
                    print(out)
                    jj = 0
                    out = "\t"
            if jj>0:
                print(out)
            raise Exception("x mesh spacing is too large")
        
        return
        
    

    # ------ Utility Methods ------

    def _debrief_msg(self, message, threshold):
        '''
        Prints debrief message if debrief is over threshold
        '''
        if self.debrief > threshold:
            print(self.prompt+message)