import numpy as np
import copy

from .reverse import reverse
from .make_dvr_dvx import make_dvr_dvx
from .create_shifted_maxwellian_include import create_shifted_maxwellian_include
from .kinetic_mesh import kinetic_mesh

from .collrad_sigmav_ion_h0 import collrad_sigmav_ion_h0
from .jh_related.jhs_coef import jhs_coef
from .sigma.sigmav_ion_h0 import sigmav_ion_h0
from .jh_related.jhalpha_coef import jhalpha_coef
from .sigma.sigmav_rec_h1s import sigmav_rec_h1s
from .sigma.sigma_cx_h0 import sigma_cx_h0
from .sigma.sigma_el_h_h import sigma_el_h_h
from .sigma.sigma_el_h_hh import sigma_el_h_hh
from .sigma.sigma_el_p_h import sigma_el_p_h
from .sigma.sigmav_cx_h0 import sigmav_cx_h0

from .sign import sign #NOTE replace sign with np.sign
from .sval import sval

from .common import constants as CONST
from .common.JH_Coef import JH_Coef
from .common.Kinetic_H import Kinetic_H_Common

# This subroutine is part of the "KN1D" atomic and molecular neutral transport code.

#   This subroutine solves a 1-D spatial, 2-D velocity kinetic neutral transport 
# problem for atomic hydrogen (H) or deuterium by computing successive generations of 
# charge exchange and elastic scattered neutrals. The routine handles electron-impact 
# ionization, proton-atom charge exchange, radiative recombination, and elastic
# collisions with hydrogenic ions, neutral atoms, and molecules.

#   The positive vx half of the atomic neutral distribution function is inputted at x(0) 
# (with arbitrary normalization) and the desired flux of hydrogen atoms entering the slab,
# at x(0) is specified. Background profiles of plasma ions, (e.g., Ti(x), Te(x), n(x), vxi(x),...)
# molecular ions, (nHP(x), THP(x)), and molecular distribution function (fH) are inputted.

# Optionally, the hydrogen source velocity distribution function is also inputted.
# (The H source and fH2 distribution functions can be computed using procedure 
# "Kinetic_H2.pro".) The code returns the atomic hydrogen distribution function, fH(vr,vx,x) 
# for all vx, vr, and x of the specified vr,vx,x grid.

#   Since the problem involves only the x spatial dimension, all distribution functions
# are assumed to have rotational symmetry about the vx axis. Consequently, the distributions
# only depend on x, vx and vr where vr =sqrt(vy^2+vz^2)

#  History:

#    B. LaBombard   First coding based on Kinetic_Neutrals.pro 		22-Dec-2000

#    For more information, see write-up: "A 1-D Space, 2-D Velocity, Kinetic 
#    Neutral Transport Algorithm for Hydrogen Atoms in an Ionizing Plasma", B. LaBombard

# Note: Variable names contain characters to help designate species -
#	atomic neutral (H), molecular neutral (H2), molecular ion (HP), proton (i) or (P) 

def kinetic_h(mesh : kinetic_mesh, mu, vxi, fHBC, GammaxHBC, fH2, fSH, nHP, THP, jh_coeffs : JH_Coef, KH_Common : Kinetic_H_Common, fH = None,
              truncate = 1e-4, Compute_Errors = 0, plot = 0, debug = 0, pause = 0, debrief = 0,
              Simple_CX = 1, Max_Gen = 50, No_Johnson_Hinnov = 0, Use_Collrad_Ionization = 0,
              No_Recomb = 0, H_H_EL = 0, H_P_EL = 0, _H_H2_EL = 0, H_P_CX = 0, ni_correct = 0): # changed fH default to None and Use_Collrad_Ionization capitalization

    vx = mesh.vx
    vr = mesh.vr
    x = mesh.x
    Tnorm = mesh.Tnorm
    Ti = mesh.Ti
    Te = mesh.Te
    n = mesh.ne
    PipeDia = mesh.PipeDia

    #	Input:
    #		vx(*)	- fltarr(nvx), normalized x velocity coordinate 
    #			[negative values, positive values],
    #			monotonically increasing. Note: a nonuniform mesh can be used.
    #			Dimensional velocity (note: Vth is based on ATOM mass)
    #			is v = Vth * vx where Vth=sqrt(2 k Tnorm/(mH*mu))
    #			Note: nvx must be even and vx(*) symmetric about 
    #			zero but not contain a zero element
    #		vr(*)	- fltarr(nvr), normalized radial velocity coordinate 
    #			[positive values], monotonically increasing. Note: a non-uniform mesh can be used.
    #			Dimensional velocity is v = Vth * vr where Vth=sqrt(2 k Tnorm/(mH*mu)) 
    #			Note: vr must not contain a zero element
    #		x(*)	- fltarr(nx), spatial coordinate (meters), 
    #			positive, monontonically increasing. Note: a non-uniform mesh can be used.
    #		Tnorm	- Float, temperature corresponding to the thermal speed (see vx and vr above) (eV)
    #	    mu	- Float, 1=hydrogen, 2=deuterium
    #	    Ti	- fltarr(nx), Ion temperature profile (eV)
    #	    Te	- fltarr(nx), electron temperature profile (eV)
    #	    n	- fltarr(nx), electron density profile (m^-3)
    #	    vxi	- fltarr(nx), x-directed plasma ion and molecular ion flow profile (m s^-1)
    #		fHBC	- fltarr(nvr,nvx), this is an input boundary condition
    #			specifying the shape of the neutral atom velocity distribution 
    #			function at location x(0). Normalization is arbitrary.
    #	        Only values with positive vx, fHBC(*,nvx/2:*) are used
    #	        by the code.
    #		GammaxHBC	- float, desired neutral atom flux density in the +Vx
    #			direction at location x(0) (m^-2 s^-1)
    #			fHBC is scaled to yield this flux density.
    #       PipeDia	- fltarr(nx), effective pipe diameter (meters)
    #			This variable allows collisions with the 'side-walls' to be simulated.
    #			If this variable is undefined, then PipeDia set set to zero. Zero values
    #			of PipeDia are ignored (i.e., treated as an infinite diameter).
    #       fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
    #           function. fH2 is normalized so that the molecular neutral density, nH2(k), is 
    #			defined as the velocity space integration: nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
    #           If this variable is undefined, then it is set equal to zero and
    #           no molecule-atom collisions are included.
    #			Note: dVx is velocity space differential for Vx axis and Vr2pidVr = Vr*!pi*dVr
    #	        with dVr being velocity space differential for Vr axis.
    #       fSH	- fltarr(nvr,nvx,nx), atomic hydrogen source velocity distribution.
    #           fSH must be normalized so that the total atomic neutral
    #           source, SourceH(k), is defined as the velocity space integration:
    #           SourceH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
    #		fSH can be computed from IDL procedure Kinetic_H2.pro
    #           If this variable is undefined, then it is set equal to zero.
    #       nHP	- fltarr(nx), molecular ion density profile (m^-3)
    #           If this parameter is undefined, then it is set equal to zero.
    #		nHP can be computed from IDL procedure Kinetic_H2.pro
    #       THP	- fltarr(nx), molecular ion temperature profile (m^-3)
    #			If this parameter is undefined, then it is set equal to 3 eV at each grid point.
    #			THP can be computed from IDL procedure Kinetic_H2.pro

    #	Input & Output:
    #		fH	- fltarr(nvr,nvx,nx), neutral atom velocity distribution
    #			function. 'Seed' values for this may be specified on input. 
    #	        If this parameter is undefined on input, then a zero 'seed' value will be used. 
    #			The algorithm outputs a self-consistent fH.
    #			fH is normalized so that the neutral density, nH(k), is defined as 
    #			the velocity space integration: nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))

    #   Output:
    #       nH	- fltarr(nx), neutral atom density profile (m^-3)
    #       GammaxH	- fltarr(nx), neutral atom flux profile (# m^-2 s^-1)
    #           computed from GammaxH(k)=Vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))
    #       VxH	- fltarr(nx), neutral atom velocity profile (m s^-1)
    #           computed from GammaxH/nH
    #           To aid in computing the some of the quantities below, the procedure internally
    #           defines the quantities:
    #           vr2vx2_ran(i,j,k)=vr(i)^2+(vx(j)-VxH(k))^2
    #           which is the magnitude of 'random v^2' at each mesh point
    #           vr2vx2(i,j,k)=vr(i)^2+vx(j)^2
    #           which is the magnitude of 'total v^2' at each mesh point
    #           q=1.602177D-19, mH=1.6726231D-27
    #           C(*,*,*) is the right hand side of the Boltzmann equation, evaluated
    #           using the computed neutral distribution function
    #       pH	- fltarr(nx), neutral atom pressure (eV m^-2) computed from:
    #           pH(k)~vth2*total(Vr2pidVr*(vr2vx2_ran(*,*,k)*fH(*,*,k))#dVx))*(mu*mH)/(3*q)
    #       TH	- fltarr(nx), neutral atom temperature profile (eV) computed from: TH=pH/nH
    #       qxH	- fltarr(nx), neutral atom random heat flux profile (watts m^-2) computed from:
    #           qxH(k)~vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH(*,*,k))#(dVx*(vx-VxH(k)))))*0.5*(mu*mH)
    #       qxH_total	- fltarr(nx), total neutral atom heat flux profile (watts m^-2)
    #           This is the total heat flux transported by the neutrals:
    #           qxH_total=(0.5*nH*(mu*mH)*VxH*VxH + 2.5*pH*q)*VxH + piH_xx*VxH + qxH
    #       NetHSource	- fltarr(nx), net H0 source [H0 source - ionization sink - wall sink] (m^-3 s^-1) computed from
    #			NetHSource(k)=total(Vr2pidVr*(C(*,*,k)#dVx))
    #	    Sion	- fltarr(nx), H ionization rate (m^-3 s^-1) 
    #       QH	- fltarr(nx), rate of net thermal energy transfer into neutral atoms (watts m^-3) computed from
    #           QH(k)~vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C(*,*,k))#dVx))*0.5*(mu*mH)
    #       RxH	- fltarr(nx), rate of x momentum transfer to neutral atoms (=force, N m^-2).
    #           RxH(k)~Vth*total(Vr2pidVr*(C(*,*,k)#(dVx*(vx-VxH(k)))))*(mu*mH)
    #       QH_total	- fltarr(nx), net rate of total energy transfer into neutral atoms
    #           = QH + RxH*VxH - 0.5*(mu*mH)*(Sloss-SourceH)*VxH*VxH (watts m^-3)
    #       AlbedoH	- float, Ratio of atomic neutral particle flux with Vx < 0 divided by particle flux
    #           with Vx > 0  at x=x(0)
    #           (Note: For fSH non-zero, the flux with Vx < 0 will include
    #           contributions from molecular hydrogen sources within the 'slab'.
    #           In this case, this parameter does not return the true 'Albedo'.)
    #       WallH	- fltarr(nx), atomic neutral sink rate arising from hitting the 'side walls' (m^-3 s^-1)
    #		    Unlike the molecules in Kinetic_H2, wall collisions result in the destruction of atoms.
    #		    This parameter can be used to specify a resulting source of molecular
    #		    neutrals in Kinetic_H2. (molecular source = 2 times WallH)

    #   KEYWORDS:
    #       Output:
    #           error	- Returns error status: 
    #               0=no error, solution returned
    #			    1=error, no solution returned

    #   COMMON BLOCK Kinetic_H_OUTPUT
    #       Output:
    #           piH_xx	- fltarr(nx), xx element of stress tensor (eV m^-2) computed from:
    #               piH_xx(k)~vth2*total(Vr2pidVr*(fH(*,*,k)#(dVx*(vx-VxH(k))^2)))*(mu*mH)/q - pH
    #           piH_yy	- fltarr(nx), yy element of stress tensor (eV m^-2) computed from:
    #               piH_yy(k)~vth2*total((Vr2pidVr*Vr^2)*(fH(*,*,k)#dVx))*(mu*mH)/q - pH
    #           piH_zz	- fltarr(nx), zz element of stress tensor (eV m^-2) = piH_yy
    #		        Note: cylindrical system relates r^2 = y^2 + z^2. All other stress tensor elements are zero.

    #       The following momentum and energy transfer rates are computed from charge-exchange collsions between species:
    #           RxHCX	- fltarr(nx), rate of x momentum transfer from hydrogren ions to atoms (=force/vol, N m^-3).
    #           EHCX	- fltarr(nx), rate of energy transfer from hydrogren ions to atoms (watts m^-3).
               
    #       The following momentum and energy transfer rates are computed from elastic collsions between species:
    #           RxH2_H	- fltarr(nx), rate of x momentum transfer from neutral molecules to atoms (=force/vol, N m^-3).
    #           RxP_H	- fltarr(nx), rate of x momentum transfer from hydrogen ions to neutral atoms (=force/vol, N m^-3).
    #           EH2_H	- fltarr(nx), rate of energy transfer from neutral molecules to atoms (watts m^-3).
    #           EP_H	- fltarr(nx), rate of energy transfer from hydrogen ions to neutral atoms (watts m^-3).

    #       The following momentum and energy transfer rates are computed from collisions with the 'side-walls'
    #           RxW_H	- fltarr(nx), rate of x momentum transfer from wall to neutral atoms (=force/vol, N m^-3).
    #           EW_H	- fltarr(nx), rate of energy transfer from wall to neutral atoms (watts m^-3).

    #       The following is the rate of parallel to perpendicular energy transfer computed from elastic collisions
    #           Epara_PerpH_H	- fltarr(nx), rate of parallel to perp energy transfer within atomic hydrogen species (watts m^-3).

    #       Source/Sink info:
    #           SourceH	- fltarr(nx), source rate of neutral atoms from H2 dissociation (from integral of inputted fSH) (m^-3 s^-1).
    #           SRecom	- fltarr(nx), source rate of neutral atoms from recombination (m^-3 s^-1).

    #   KEYWORDS:
    #       Input:
    #           truncate	- float, stop computation when the maximum 
    #		        increment of neutral density normalized to 
    #		        inputed neutral density is less than this 
    #    		    value in a subsequent generation. Default value is 1.0e-4
    #           Simple_CX	- if set, then use CX source option (B): Neutrals are born
    #               in velocity with a distribution proportional to the local
    #               ion distribution function. Simple_CX=1 is default.
    #               if not set, then use CX source option (A): The CX source
    #               neutral atom distribution function is computed by evaluating the
    #               the CX cross section for each combination of (vr,vx,vr',vx')
    #               and convolving it with the neutral atom distribution function.
    #               This option requires more CPU time and memory.
    #  	  	    Max_gen	- integer, maximum number of collision generations to try including before giving up.
    #               Default is 50.
    #           No_Johnson_Hinnov	- if set, then compute ionization and recombination rates
    #		        directly from reaction rates published by Janev* for
    #		        ground state hydrogen
    #		        Ionization:    e + H(1s) -> p + e 
    #		        Recombination: e + p -> H(1s) + hv
    #		        *Janev, R.K., et al, "Elementary processes in hydrogen-helium plasmas",
    #		        (Springer-Verlag, Berlin ; New York, 1987)
    #		        Otherwise, compute ionization and recombination rates using
    #	            results from the collisional-radiative model published by Johnson
    #		        and Hinnov [L.C.Johnson and E. Hinnov, J. Quant. Spectrosc. Radiat.
     #		        Transfer. vol. 13 pp.333-358]. This is the default.
    #		        Note: charge exchange is always computed using the ground state reaction
    #	            rates published by Janev:
    #		        Charge Exchange: p + H(1s) -> H(1s) + p
    #           Use_Collrad_Ionization - FS - if set, this overrides
    #               No_Johnson_Hinnov and uses rates from the
    #               COLLRAD code *for ionization only*
    #           No_Recomb	- if set, then DO NOT include recombination as a source of atomic neutrals
    #	            in the algorithm
    # 	        H_H_EL	- if set, then include H -> H elastic self collisions
    #		        Note: if H_H_EL is set, then algorithm iterates fH until
    #               self consistent fH is achieved.
    # 	        H_P_CX	- if set, then include H -> H(+) charge exchange collisions 
    #           H_P_EL	- if set, then include H -> H(+) elastic collisions 
    #           H_H2_EL	- if set, then include H -> H2 elastic collisions 
    #           ni_correct	- if set, then algorithm corrects hydrogen ion density
    #		        according to quasineutrality: ni=ne-nHP. Otherwise, nHP is assumed to be small.
    #           Compute_Errors	- if set, then return error estimates in common block Kinetic_H_ERRORS below
    #	        plot	- 0= no plots, 1=summary plots, 2=detail plots, 3=very detailed plots
    #	        debug	- 0= do not execute debug code, 1=summary debug, 2=detail debug, 3=very detailed debug
    #           debrief	- 0= do not print, 1=print summary information, 2=print detailed information
    #           pause	- if set, then pause between plots

    #	COMMON BLOCK Kinetic_H_ERRORS

    #		if COMPUTE_ERRORS keyword is set then the following is returned in common block Kinetic_H_ERRORS

    #			Max_dx	- float(nx), Max_dx(k) for k=0:nx-2 returns maximum 
    #				allowed x(k+1)-x(k) that avoids unphysical negative 
    #				contributions to fH
    #			Vbar_error	- float(nx), returns numerical error in computing
    #				the speed of ions averged over maxwellian distribution.
    #				The average speed should be:
    #				vbar_exact=2*Vth*sqrt(Ti(*)/Tnorm)/sqrt(!pi)
    #				Vbar_error returns: abs(vbar-vbar_exact)/vbar_exact
    #				where vbar is the numerically computed value.
    #			mesh_error	- fltarr(nvr,nvx,nx), normalized error of solution
    #				based on substitution into Boltzmann equation.
    #			moment_error	- fltarr(nx,m), normalized error of solution
    #				based on substitution into velocity space
    #				moments (v^m) of Boltzmann equation, m=[0,1,2,3,4]
    #  	        C_error	- fltarr(nx), normalized error in charge exchange and elastic scattering collision 
    #			    operator. This is a measure of how well the charge exchange and
    #			    elastic scattering portions of the collision operator
    #			    conserve particles.
    #  			CX_error	- fltarr(nx), normalized particle conservation error in charge exchange collision operator.
    #			H_H_error	-  fltarr(nx,[0,1,2]) return normalized errors associated with 
    #				particle [0], x-momentum [1], and total energy [2] convervation of the elastic self-collision operator

    #			qxH_total_error	- fltarr(nx), normalized error estimate in computation of qxH_total
    #			QH_total_error	- fltarr(nx), normalized error estimate in computation of QH_total

    #	History:
    #		22-Dec-2000 - B. LaBombard - first coding.
    #		11-Feb-2001 - B. LaBombard - added elastic collisions 

    prompt = 'Kinetic_H => '

    #	Set Kinetic_H_input common block variables 
    vx_s = KH_Common.Input.vx_s
    vr_s = KH_Common.Input.vr_s
    x_s = KH_Common.Input.x_s
    Tnorm_s = KH_Common.Input.Tnorm_s
    mu_s = KH_Common.Input.mu_s
    Ti_s = KH_Common.Input.Ti_s
    Te_s = KH_Common.Input.Te_s
    n_s = KH_Common.Input.n_s
    vxi_s = KH_Common.Input.vxi_s
    fHBC_s = KH_Common.Input.fHBC_s
    GammaxHBC_s = KH_Common.Input.GammaxHBC_s # Fixed typo - changed to GammaxHBC_s
    PipeDia_s = KH_Common.Input.PipeDia_s
    fH2_s = KH_Common.Input.fH2_s
    fSH_s = KH_Common.Input.fSH_s
    nHP_s = KH_Common.Input.nHP_s
    THP_s = KH_Common.Input.THP_s # Fixed typo - changed to THP_s
    fH_s = KH_Common.Input.fH_s # Fixed typo - changed to fH_s
    Simple_CX_s = KH_Common.Input.Simple_CX_s
    JH_s = KH_Common.Input.JH_s
    Collrad_s = KH_Common.Input.Collrad_s
    Recomb_s = KH_Common.Input.Recomb_s
    H_H_EL_s = KH_Common.Input.H_H_EL_s
    H_P_EL_s = KH_Common.Input.H_P_EL_s
    H_H2_EL_s = KH_Common.Input.H_H2_EL_s
    H_P_CX_s = KH_Common.Input.H_P_CX_s
    #	FS: added collrad_s

    #	Kinetic_H_internal common block

    vr2vx2 = KH_Common.Internal.vr2vx2
    vr2vx_vxi2 = KH_Common.Internal.vr2vx_vxi2
    fi_hat = KH_Common.Internal.fi_hat
    ErelH_P = KH_Common.Internal.ErelH_P
    Ti_mu = KH_Common.Internal.Ti_mu
    ni = KH_Common.Internal.ni
    sigv = KH_Common.Internal.sigv
    alpha_ion = KH_Common.Internal.alpha_ion
    v_v2 = KH_Common.Internal.v_v2
    v_v = KH_Common.Internal.v_v
    vr2_vx2 = KH_Common.Internal.vr2_vx2
    vx_vx = KH_Common.Internal.vx_vx
    Vr2pidVrdVx = KH_Common.Internal.Vr2pidVrdVx
    SIG_CX = KH_Common.Internal.SIG_CX # Fixed typo - changed to SIG_CX
    SIG_H_H = KH_Common.Internal.SIG_H_H
    SIG_H_H2 = KH_Common.Internal.SIG_H_H2
    SIG_H_P = KH_Common.Internal.SIG_H_P
    Alpha_CX = KH_Common.Internal.Alpha_CX
    Alpha_H_H2 = KH_Common.Internal.Alpha_H_H2
    Alpha_H_P = KH_Common.Internal.Alpha_H_P
    MH_H_sum = KH_Common.Internal.MH_H_sum
    Delta_nHs = KH_Common.Internal.Delta_nHs
    Sn = KH_Common.Internal.Sn
    Rec = KH_Common.Internal.Rec

    #	Kinetic_H_Moments
    
    nH2 = KH_Common.Moments.nH2
    vxH2 = KH_Common.Moments.VxH2
    TH2 = KH_Common.Moments.TH2 # changed to fit current global_vars structure

    #	Internal Debug switches

    Shifted_Maxwellian_Debug = 0
    CI_Test = 1
    Do_Alpha_CX_Test = 0

    #	Internal Tolerances

    DeltaVx_tol = .01
    Wpp_tol = .001

    #	Test input parameters

    if debug>0:
        #plot = np.max(plot, 1)
        debrief = np.maximum(debrief, 1)
        pause = 1
    JH = 1
    if No_Johnson_Hinnov:
        JH = 0
    Recomb = 1
    if No_Recomb:
        Recomb = 0
    error = 0

    nvr = vr.size
    nvx = vx.size
    nx = x.size

    dx = x-np.roll(x,1)
    dx = dx[1:]
    notpos = np.argwhere(dx <= 0)
    if notpos.size > 0:
        raise Exception(prompt+'x[*] must be increasing with index!')
        # error = 1
    if nvx % 2 != 0:
        raise Exception(prompt+'Number of elements in vx must be even!') 
        # error = 1
    if Ti.size != nx:
        raise Exception(prompt+'Number of elements in Ti and x do not agree!')
        # error = 1
    if vxi is None:
        vxi = np.zeros(nx)
    if vxi.size != nx:
        raise Exception(prompt+'Number of elements in vxi and x do not agree!')
        # error = 1
    if Te.size != nx:
        raise Exception(prompt+'Number of elements in Te and x do not agree!')
        # error = 1
    if n.size != nx:
        raise Exception('Number of elements in n and x do not agree!')
        # error = 1
    if GammaxHBC is None:
        raise Exception(prompt+'GammaxHBC is not defined!')
    if PipeDia is None:
        PipeDia = np.zeros(nx)
    if PipeDia.size != nx:
        raise Exception('Number of elements in PipeDia and x do not agree!') # Fixed error message- previously copied from n.size!=nx
        # error = 1
    if len(fHBC[:,0]) != nvr:
        raise Exception(prompt+'Number of elements in fHBC[:,0] and vr do not agree!')
        # error = 1
    if len(fHBC[0,:]) != nvx:
        raise Exception(prompt+'Number of elements in fHBC[0,:] and vx do not agree!')
        # error = 1
    if fH2 is None:
        fH2 = np.zeros((nvr, nvx, nx))
    if fH2[:,0,0].size != nvr:
        raise Exception(prompt+'Number of elements in fH2[:,0,0] and vr do not agree!')
        # error = 1
    if fH2[0,:,0].size != nvx:
        raise Exception(prompt+'Number of elements in fH2[0,:,0] and vx do not agree!')
        # error = 1
    if fH2[0,0,:].size != nx:
        raise Exception(prompt+'Number of elements in fH2[0,0,:] and x do not agree!')
        # error = 1
    if fSH is None:
        fSH = np.zeros((nvr, nvx, nx))
    if fSH[:,0,0].size != nvr:
        raise Exception(prompt+'Number of elements in fSH[:,0,0] and vr do not agree!')
        # error = 1
    if fSH[0,:,0].size != nvx:
        raise Exception(prompt+'Number of elements in fSH[0,:,0] and vx do not agree!')
        # error = 1
    if fSH[0,0,:].size != nx:
        raise Exception(prompt+'Number of elements in fSH[0,0,:] and x do not agree!')
        # error = 1
    if nHP is None:
        nHP = np.zeros(nx)
    if nHP.size != nx:
        raise Exception(prompt+'Number of elements in nHP and x do not agree!')
        # error = 1
    if THP is None:
        THP = np.full(nx, 1.0)
    if THP.size != nx:
        raise Exception(prompt+'Number of elements in nHP and x do not agree!')
        # error = 1
    if fH is None:
        fH = np.zeros((nvr,nvx,nx))
    if fH[:,0,0].size != nvr:
        raise Exception(prompt+'Number of elements in fH[:,0,0] and vr do not agree!')
        # error = 1
    if fH[0,:,0].size != nvx:
        raise Exception(prompt+'Number of elements in fH[0,:,0] and vx do not agree!')
        # error = 1
    if fH[0,0,:].size != nx:
        raise Exception(prompt+'Number of elements in fH[0,0,:] and x do not agree!')
        # error = 1
    if np.sum(abs(vr)) <= 0:
        raise Exception(prompt+'vr is all 0!')
        # error = 1
    ii = np.argwhere(vr <= 0)
    if ii.size > 0:
        raise Exception(prompt+'vr contains zero or negative element(s)!')
        # error = 1
    if np.sum(abs(vx)) <= 0:
        raise Exception(prompt+'vx is all 0!')
        # error = 1
    if np.sum(x) <= 0:
        raise Exception(prompt+'Total(x) is less than or equal to 0!')
        # error = 1
    if mu is None:
        raise Exception(prompt+'mu is not defined!')
    if mu not in [1,2]:
        raise Exception(prompt+'mu must be 1 or 2!')
        # error = 1

    #NOTE Removed Plotting formatting, bring back once the program actually works

    i_n = np.argwhere(vx < 0)
    if np.size(i_n) < 1:
        print(prompt+'vx contains no negative elements!')
        # error = 1
    i_p = np.argwhere(vx > 0)
    if np.size(i_p) < 1:
        print(prompt+'vx contains no positive elements!')
        # error = 1
    i_z = np.argwhere(vx == 0)
    if np.size(i_z) > 0:
        print(prompt+'vx contains one or more zero elements!')
        # error = 1
    diff = np.argwhere(vx[i_p] != -np.flipud(vx[i_n]))
    if diff.size > 0:
        raise Exception(prompt + " vx array elements are not symmetric about zero!")
        # error = 1
    # print("i_p", i_p)
    # print("i_n", i_n)
    # print("i_z", i_z)
    # input()


    fHBC_input = np.zeros(fHBC.shape)
    fHBC_input[:,i_p] = fHBC[:,i_p]
    test = np.sum(fHBC_input)
    if test <= 0.0 and abs(GammaxHBC) > 0:
        raise Exception(prompt+'Values for fHBC[:,:] with vx > 0 are all zero!')
        # error = 1

    #	Output variables

    nH = np.zeros(nx)
    GammaxH = np.zeros(nx)
    VxH = np.zeros(nx)
    pH = np.zeros(nx)
    TH = np.zeros(nx)
    qxH = np.zeros(nx)
    qxH_total = np.zeros(nx)
    NetHSource = np.zeros(nx)
    WallH = np.zeros(nx)
    Sion = np.zeros(nx)
    QH = np.zeros(nx)
    RxH = np.zeros(nx)
    QH_total = np.zeros(nx)
    piH_xx = np.zeros(nx)
    piH_yy = np.zeros(nx)
    piH_zz = np.zeros(nx)
    RxHCX = np.zeros(nx)
    RxH2_H = np.zeros(nx)
    RxP_H = np.zeros(nx)
    RxW_H = np.zeros(nx)
    EHCX = np.zeros(nx)
    EH2_H = np.zeros(nx)
    EP_H = np.zeros(nx)
    EW_H = np.zeros(nx)
    Epara_PerpH_H = np.zeros(nx)
    AlbedoH = 0.0e0
    SourceH = np.zeros(nx)
    SRecomb = np.zeros(nx)

    #	Internal variables

    Work = np.zeros((nvr*nvx))
    fHG = np.zeros((nvr,nvx,nx))
    NHG = np.zeros((nx,Max_Gen+1))
    Vth = np.sqrt((2*CONST.Q*Tnorm) / (mu*CONST.H_MASS))
    Vth2 = Vth*Vth
    Vth3 = Vth2*Vth
    fHs = np.zeros(nx)
    nHs = np.zeros(nx)
    Alpha_H_H = np.zeros((nvr,nvx))
    Omega_H_P = np.zeros(nx)
    Omega_H_H2 = np.zeros(nx)
    Omega_H_H = np.zeros(nx)
    VxHG = np.zeros(nx)
    THG = np.zeros(nx)
    Wperp_paraH = np.zeros(nx)
    vr2vx2_ran2 = np.zeros((nvr,nvx))
    vr2_2vx_ran2 = np.zeros((nvr,nvx))
    vr2_2vx2_2D = np.zeros((nvr,nvx))
    RxCI_CX = np.zeros(nx)
    RxCI_H2_H = np.zeros(nx)
    RxCI_P_H = np.zeros(nx)
    Epara_Perp_CI = np.zeros(nx)
    CI_CX_error = np.zeros(nx)
    CI_H2_H_error = np.zeros(nx)
    CI_P_H_error = np.zeros(nx)
    CI_H_H_error = np.zeros(nx)
    Maxwell = np.zeros((nvr,nvx,nx))

    Vr2pidVr, VrVr4pidVr, dVx, vrL, vrR, vxL, vxR, vol, Vth_DeltaVx, Vx_DeltaVx, Vr_DeltaVr, Vr2Vx2_2D, jpa, jpb, jna, jnb = make_dvr_dvx(vr,vx)
    # print("Vr2pidVr", Vr2pidVr)
    # print("VrVr4pidVr", VrVr4pidVr)
    # print("dVx", dVx)
    # print("vol", vol.T)
    # print("Vth_DeltaVx", Vth_DeltaVx.T)
    # print("Vx_DeltaVx", Vx_DeltaVx.T)
    # print("Vr_DeltaVr", Vr_DeltaVr.T)
    # print("jpa", jpa)
    # print("jpb", jpb)
    # print("jna", jna)
    # print("jnb", jnb)
    # input()

    #	Vr^2-2*Vx^2

    for i in range(nvr):
        vr2_2vx2_2D[i,:] = (vr[i]**2) - 2*(vx**2)
    # print("vr2_2vx2_2D", vr2_2vx2_2D.T)
    # input()

    #	Theta-prime coordinate

    # Theta-prime Coordinate
    ntheta = 5 # use 5 theta mesh points for theta integration
    dTheta = np.ones(ntheta)/ntheta
    theta = np.pi*(np.arange(ntheta)/ntheta + 0.5/ntheta)
    # print("theta", theta)
    cos_theta = np.cos(theta)
    # print("cos", cos_theta)
    # input()

    #	Scale input molecular distribution function to agree with desired flux
    gamma_input = 1.0
    if abs(GammaxHBC) > 0:
        gamma_input = Vth*np.sum(Vr2pidVr*(fHBC_input @ (vx*dVx)))
    ratio = abs(GammaxHBC)/gamma_input
    fHBC_input = fHBC_input*ratio
    if abs(ratio - 1) > 0.01*truncate:
        fHBC = fHBC_input
    fH[:,i_p,0] = fHBC_input[:,i_p]
    # print("fH", fH[:,i_p,0].T)
    # print("shape", fH.shape)
    # input()

    #	if fH2 is zero, then turn off elastic H2 <-> H collisions
    H_H2_EL = _H_H2_EL
    if np.sum(fH2) <= 0:
        H_H2_EL = 0

    #	Set iteration scheme
    fH_iterate = 0
    if (H_H_EL != 0) or (H_P_EL != 0) or (H_H2_EL != 0): 
        fH_iterate = 1

    fH_generations = 0
    if (fH_iterate != 0) or (H_P_CX != 0): 
        fH_generations = 1

    #	Set flags to make use of previously computed local parameters 
    New_Grid = 1
    if vx_s is not None:
        test = 0 
        ii = np.argwhere(vx_s != vx).T ; test = test + np.size(ii)
        ii = np.argwhere(vr_s != vr).T ; test = test + np.size(ii)
        ii = np.argwhere(x_s != x).T ; test = test + np.size(ii)
        ii = np.argwhere(Tnorm_s != Tnorm).T ; test = test + np.size(ii)
        ii = np.argwhere(mu_s != mu).T ; test = test + np.size(ii)
        if test <= 0:
            New_Grid = 0
    New_Protons = 1
    if Ti_s is not None:
        test = 0 
        ii = np.argwhere(Ti_s != Ti).T ; test = test + np.size(ii)
        ii = np.argwhere(n_s != n).T ; test = test + np.size(ii)
        ii = np.argwhere(vxi_s != vxi).T ; test = test + np.size(ii)
        if test <= 0:
            New_protons = 0
    New_Molecular_Ions = 1
    if nHP_s is not None:
        test = 0
        test = np.argwhere(nHP_s != nHP).T ; test = test + np.size(ii)
        test = np.argwhere(THP_s != THP).T ; test = test + np.size(ii)
        if test <= 0:
            New_Molecular_Ions = 0
    New_Electrons = 1
    if Te_s is not None:
        test = 0 
        ii = np.argwhere(Te_s != Te).T ; test = test + np.size(ii)
        ii = np.argwhere(n_s != n).T ; test = test + np.size(ii)
        if test <= 0:
            New_Electrons = 0
    New_fH2 = 1
    if fH2_s is not None:
        ii = np.argwhere(fH2_s != fH2)
        if np.size(ii) <= 0:
            New_fH2 = 0
    New_fSH = 1
    if fSH_s is not None:
        ii = np.argwhere(fSH_s != fSH)
        if np.size(ii) <= 0:
            New_fSH = 0
    New_Simple_CX = 1
    if Simple_CX_s is not None:
        ii = np.argwhere(Simple_CX_s != Simple_CX)
        if np.size(ii) <= 0:
            New_Simple_CX = 0
    New_H_Seed = 1
    if fH_s is not None:
        ii = np.argwhere(fH_s != fH)
        if np.size(ii) <= 0:
            New_H_Seed = 0

    Do_sigv = 			New_Grid | New_Electrons
    Do_ni = 			New_Grid | New_Electrons | New_Protons | New_Molecular_Ions
    Do_fH2_moments = 	(New_Grid | New_fH2) & (np.sum(fH2) > 0.0)
    Do_Alpha_CX = 		(New_Grid | (Alpha_CX is None) | Do_ni | New_Simple_CX) and H_P_CX
    Do_SIG_CX = 		(New_Grid | (SIG_CX is None) | New_Simple_CX) & (Simple_CX == 0) & Do_Alpha_CX
    Do_Alpha_H_H2 = 	(New_Grid | (Alpha_H_H2 is None) | New_fH2) & H_H2_EL
    Do_SIG_H_H2 = 		(New_Grid | (SIG_H_H2 is None)) & Do_Alpha_H_H2
    Do_SIG_H_H = 		(New_Grid | (SIG_H_H is None)) & H_H_EL
    Do_Alpha_H_P = 		(New_Grid | (Alpha_H_P is None) | Do_ni) & H_P_EL
    Do_SIG_H_P = 		(New_Grid | (SIG_H_P is None)) & Do_Alpha_H_P
    Do_v_v2 = 			(New_Grid | (v_v2 is None)) & (CI_Test | Do_SIG_CX | Do_SIG_H_H2 | Do_SIG_H_H | Do_SIG_H_P)
    
    print("Kinetic H Settings")
    print("H_H_EL", H_H_EL)
    print("H_P_EL", H_P_EL)
    print("_H_H2_EL", _H_H2_EL)
    print("H_P_CX", H_P_CX)
    print("New_Grid", New_Grid)
    print("New_Protons", New_Protons)
    print("New_Molecular_Ions", New_Molecular_Ions)
    print("New_Electrons", New_Electrons)
    print("New_fH2", New_fH2)
    print("New_fSH", New_fSH)
    print("New_Simple_CX", New_Simple_CX)
    print("New_H_Seed", New_H_Seed)
    print("Do_sigv", Do_sigv)
    print("Do_ni", Do_ni)
    print("Do_fH2_moments", Do_fH2_moments)
    print("Do_Alpha_CX", Do_Alpha_CX)
    print("Do_SIG_CX", Do_SIG_CX)
    print("Simple_CX", Simple_CX)
    print("Do_Alpha_H_H2", Do_Alpha_H_H2)
    print("Do_SIG_H_H2", Do_SIG_H_H2)
    print("Do_SIG_H_H", Do_SIG_H_H)
    print("Do_Alpha_H_P", Do_Alpha_H_P)
    print("Do_SIG_H_P", Do_SIG_H_P)
    print("Do_v_v2", Do_v_v2)
    input()

    nH2 = np.zeros(nx)
    vxH2 = np.zeros(nx)
    TH2 = np.full(nx, 1.0)

    if Do_fH2_moments:
        if debrief > 1:
            print(prompt+'Computing vx and T moments of fH2')

        #	Compute x flow velocity and temperature of molecular species
        for k in range(nx):
            nH2[k] = np.sum(Vr2pidVr*(fH2[:,:,k] @ dVx))
            if nH2[k] > 0:
                vxH2[k] = Vth*np.sum(Vr2pidVr*(fH2[:,:,k] @ (vx*dVx)))/nH2[k]
                for i in range(nvr):
                    vr2vx2_ran2[i,:] = vr[i]**2 + (vx - vxH2[k]/Vth)**2
                TH2[k] = (2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran2*fH2[:,:,k]) @ dVx)) / (3*CONST.Q*nH2[k])
        # print("TH2", TH2)
        # print("vxH2", vxH2)
        # input()

    if New_Grid:
        if debrief > 1:
            print(prompt+'Computing vr2vx2, vr2vx_vxi2, ErelH_P')

        #	Magnitude of total normalized v^2 at each mesh point
        vr2vx2 = np.zeros((nvr,nvx,nx))
        for i in range(nvr):
            for k in range(nx):
                vr2vx2[i,:,k] = vr[i]**2 + vx**2

        #	Magnitude of total normalized (v-vxi)^2 at each mesh point
        vr2vx_vxi2 = np.zeros((nvr,nvx,nx))
        for i in range(nvr):
            for k in range(nx):
                vr2vx_vxi2[i,:,k] = vr[i]**2 + (vx - vxi[k]/Vth)**2

        #	Atomic hydrogen ion energy in local rest frame of plasma at each mesh point
        ErelH_P = (0.5*CONST.H_MASS*vr2vx_vxi2*Vth2) / CONST.Q
        ErelH_P = np.maximum(ErelH_P, 0.1) # sigmav_cx does not handle neutral energies below 0.1 eV
        ErelH_P = np.minimum(ErelH_P, 2e4) # sigmav_cx does not handle neutral energies above 20 keV

        # print("ErelH_P", ErelH_P[:,10,10].T)
        # print("vr2vx2", vr2vx2[:,10,10].T)
        # print("vr2vx_vxi2", vr2vx_vxi2[:,10,10].T)
        # input()

    if New_Protons:
        if debrief>1:
            print(prompt+'Computing Ti/mu at each mesh point')

        #	Ti/mu at each mesh point
        Ti_mu = np.zeros((nvr,nvx,nx))
        for k in range(nx):
            Ti_mu[:,:,k] = Ti[k]/mu

        #	Compute Fi_hat
        if debrief>1:
            print(prompt+'Computing fi_Hat')
        vx_shift = vxi
        Tmaxwell = Ti
        mol = 1
        Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                      nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2_2D,jpa,jpb,jna,jnb)
        fi_hat = np.copy(Maxwell)
        
        # print("ti_mu", Ti_mu[:,10,10].T)
        # print("fi_hat", fi_hat[:,10,10].T)
        # input()

    if Compute_Errors:
        if debrief>1:
            print(prompt+'Computing Vbar_Error')

        #	Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti[x]/Tnorm)/sqrt(pi)

        vx_shift = np.zeros(nx)
        Tmaxwell = Ti
        mol = 1
        Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                      nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2_2D,jpa,jpb,jna,jnb)
        vbar_test = np.zeros((nvr,nvx,ntheta))
        vbar_error = np.zeros(nx)
        for m in range(ntheta):
            vbar_test[:,:,m] = vr2vx2[:,:,0]
        _vbar_test = np.zeros((nvr*nvx,ntheta))
        _vbar_test[:] = (Vth*np.sqrt(vbar_test)).reshape(_vbar_test.shape, order='F')
        vbar_test = np.zeros((nvr,nvx))
        vbar_test[:] = (_vbar_test @ dTheta).reshape(vbar_test.shape, order='F')
        for k in range(nx):
            vbar = np.sum(Vr2pidVr*((vbar_test*Maxwell[:,:,k]) @ dVx))
            vbar_exact = 2*Vth*np.sqrt(Ti[k]/Tnorm) / np.sqrt(np.pi)
            vbar_error[k] = abs(vbar-vbar_exact) / vbar_exact
        if debrief > 0:
            print(prompt+'Maximum Vbar error = ', sval(max(vbar_error)))
            # input()

    if Do_ni:
        if debrief>1:
            print(prompt+'Computing ni profile')
        ni = n
        if ni_correct:
            ni = n-nHP
        ni = np.maximum(ni, 0.01*n)
        # print("ni", ni)
        # input()

    if Do_sigv:
        if debrief>1:
            print(prompt+'Computing sigv')

        #	Compute sigmav rates for each reaction with option to use rates
        #	from CR model of Johnson-Hinnov

        sigv = np.zeros((nx,3))

        #	Reaction R1:  e + H -> e + H(+) + e   (ionization)
        #NOTE Replace Use_Collrad_Ionization with constant for consistency across program, check with someone who knows what they are doing if this is correct
        if CONST.USE_COLLRAD_IONIZATION:
            sigv[:,1] = collrad_sigmav_ion_h0(n,Te) # from COLLRAD code (DEGAS-2)
        else:
            #NOTE Replace JH with constant for consistency across program, check with someone who knows what they are doing if this is correct
            if CONST.USE_JH:
                sigv[:,1] = jhs_coef(n, Te, jh_coeffs, no_null=True) # Johnson-Hinnov, limited Te range; fixed JHS_coef capitalization #NOTE Not tested yet
            else:
                sigv[:,1] = sigmav_ion_h0(Te) # from Janev et al., up to 20keV #NOTE Not Tested Yet

        #	Reaction R2:  e + H(+) -> H(1s) + hv  (radiative recombination)
        #NOTE Replace JH with constant for consistency across program, check with someone who knows what they are doing if this is correct
        if CONST.USE_JH:
            sigv[:,2] = jhalpha_coef(n, Te, jh_coeffs, no_null=True)
        else:
            sigv[:,2] = sigmav_rec_h1s(Te)

        #	H ionization rate (normalized by vth) = reaction 1
        alpha_ion = (n*sigv[:,1]) / Vth

        #	Recombination rate (normalized by vth) = reaction 2
        Rec = (n*sigv[:,2]) / Vth

        # print("sigv", sigv.T)
        # print("alpha_ion", alpha_ion)
        # print("Rec", Rec)
        # input()

    #	Compute Total Atomic Hydrogen Source
    Sn = np.zeros((nvr,nvx,nx))

    #	Add Recombination (optionally) and User-Supplied Hydrogen Source (velocity space distribution)
    for k in range(nx):
        Sn[:,:,k] = fSH[:,:,k]/Vth
        if Recomb:
            Sn[:,:,k] = Sn[:,:,k] + fi_hat[:,:,k]*ni[k]*Rec[k]
    # print("Sn", Sn.T)
    # input()

    #	Set up arrays for charge exchange and elastic collision computations, if needed

    if Do_v_v2 == 1:
        if debrief > 1:
            print(prompt+'Computing v_v2, v_v, vr2_vx2, and vx_vx')

        #	v_v2=(v-v_prime)^2 at each double velocity space mesh point, including theta angle
        v_v2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))

        #	vr2_vx2=0.125* [ vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2 ]
        #		at each double velocity space mesh point, including theta angle
        vr2_vx2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))
        for m in range(ntheta):
            for l in range(nvx):
                for k in range(nvr):
                    for i in range(nvr):
                        v_starter = vr[i]**2 + vr[k]**2 - 2*vr[i]*vr[k]*cos_theta[m]
                        v_v2[i,:,k,l,m] = v_starter + (vx[:] - vx[l])**2
                        vr2_vx2[i,:,k,l,m] = v_starter - 2*(vx[:] - vx[l])**2

        #	v_v=|v-v_prime| at each double velocity space mesh point, including theta angle
        v_v = np.sqrt(v_v2)
        # print("V_V", v_v.T[0,0,0])
        # print("v_v shape", v_v.shape)
        # input()

        #	vx_vx=(vx-vx_prime) at each double velocity space mesh point
        vx_vx = np.zeros((nvr,nvx,nvr,nvx))
        for j in range(nvx):
            for l in range(nvx):
                vx_vx[:,j,:,l] = vx[j]-vx[l]
        # print("vx_vx", vx_vx.T)
        # input()

        #	Set Vr'2pidVr'*dVx' for each double velocity space mesh point

        Vr2pidVrdVx = np.zeros((nvr,nvx,nvr,nvx))
        for k in range(nvr):
            Vr2pidVrdVx[:,:,k,:] = Vr2pidVr[k]
        for l in range(nvx):
            Vr2pidVrdVx[:,:,:,l] = Vr2pidVrdVx[:,:,:,l]*dVx[l]
        # print("Vr2pidVrdVx", Vr2pidVrdVx.T)
        # input()

    if Simple_CX == 0 and Do_SIG_CX == 1: #NOTE Not Tested Yet
        if debrief>1:
            print(prompt+'Computing SIG_CX')

        #	Option (A) was selected: Compute SigmaV_CX from sigma directly.
        #	In preparation, compute SIG_CX for present velocity space grid, if it has not 
        #	already been computed with the present input parameters

        #	Compute sigma_cx * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx, ntheta))
        _Sig[:] = (v_v*sigma_cx_h0(v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')

        #	Set SIG_CX = vr' x Integral{v_v*sigma_cx} 
        #		over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        SIG_CX = np.zeros((nvr*nvx, nvr*nvx))
        SIG_CX[:] = (Vr2pidVrdVx*((_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F'))).reshape(SIG_CX.shape, order='F')

        #	SIG_CX is now vr' * sigma_cx(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

    if Do_SIG_H_H == 1:
        if debrief>1:
            print(prompt+'Computing SIG_H_H')

        #	Compute SIG_H_H for present velocity space grid, if it is needed and has not already been computed with the present input parameters

        #	Compute sigma_H_H * vr2_vx2 * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
        _Sig[:] = (vr2_vx2*v_v*sigma_el_h_h(v_v2*(0.5*CONST.H_MASS*mu*Vth2/CONST.Q), vis=True) / 8).reshape(_Sig.shape, order='F')
        # print("_SIG", _Sig[239,4].T)
        # input()

        #	Note: For viscosity, the cross section for D -> D is the same function of center of mass energy as H -> H.

        #	Set SIG_H_H = vr' x Integral{vr2_vx2*v_v*sigma_H_H} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        SIG_H_H = np.zeros((nvr*nvx,nvr*nvx))
        SIG_H_H[:] = (Vr2pidVrdVx*(_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F')).reshape(SIG_H_H.shape, order='F')
        # print("SIG_H_H", SIG_H_H.T)
        # input()
        #	SIG_H_H is now vr' * sigma_H_H(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

    if Do_SIG_H_H2 == 1:
        if debrief > 1:
            print(prompt+'Computing SIG_H_H2')

        #	Compute SIG_H_H2 for present velocity space grid, if it is needed and has not
        #		already been computed with the present input parameters

        #	Compute sigma_H_H2 * v_v at all possible relative velocities

        _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
        _Sig[:] = (v_v*sigma_el_h_hh(v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')

        #	NOTE: using H energy here for cross-sections tabulated as H->H2

        #	Set SIG_H_H2 = vr' x vx_vx x Integral{v_v*sigma_H_H2} over theta=0,
        #		2pi times differential velocity space element Vr'2pidVr'*dVx'

        SIG_H_H2 = np.zeros((nvr*nvx,nvr*nvx))
        SIG_H_H2[:] = (Vr2pidVrdVx*vx_vx*(_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F')).reshape(SIG_H_H2.shape, order='F')
        # print("SIG_H_H2", SIG_H_H2.T)
        # print(SIG_H_H2.shape)
        # input()

        #	SIG_H_H2 is now vr' *vx_vx * sigma_H_H2(v_v) * v_v 
        #		(intergated over theta) for all possible ([vr,vx],[vr',vx'])

    if Do_SIG_H_P == 1:
        if debrief>1:
            print(prompt+'Computing SIG_H_P')

        #	Compute SIG_H_P for present velocity space grid, if it is needed and has not 
        #		already been computed with the present input parameters

        #	Compute sigma_H_P * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
        _Sig[:] = (v_v*sigma_el_p_h(v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F') #NOTE Program slows significantly here

        #	Set SIG_H_P = vr' x vx_vx x Integral{v_v*sigma_H_P} over theta=0,
        #		2pi times differential velocity space element Vr'2pidVr'*dVx'

        SIG_H_P = np.zeros((nvr*nvx,nvr*nvx))
        SIG_H_P[:] = (Vr2pidVrdVx*vx_vx*(_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F')).reshape(SIG_H_P.shape, order='F')
        # print("SIG_H_P", SIG_H_P.T)
        # print(SIG_H_P.shape)
        # input()

        #	SIG_H_P is now vr' *vx_vx * sigma_H_P(v_v) * v_v (intergated over theta) 
        #		for all possible ([vr,vx],[vr',vx'])


    #	Compute Alpha_CX for present Ti and ni, if it is needed and has not
    #		already been computed with the present parameters
    if Do_Alpha_CX == 1:
        if debrief > 1:
            print(prompt+'Computing Alpha_CX')

        if Simple_CX:
            #	Option (B): Use maxwellian weighted <sigma v>

            #	Charge Exchange sink rate

            Alpha_CX = sigmav_cx_h0(Ti_mu, ErelH_P) / Vth
            for k in range(nx):
                Alpha_CX[:,:,k] = Alpha_CX[:,:,k]*ni[k]

        else: #NOTE Not Tested Yet
            #	Option (A): Compute SigmaV_CX from sigma directly via SIG_CX

            Alpha_CX = np.zeros((nvr,nvx,nx))
            for k in range(nx):
                Work[:] = (fi_hat[:,:,k]*ni[k]).reshape(Work.shape, order='F')
                Alpha_CX[:,:,k] = (SIG_CX @ Work).reshape(Alpha_CX[:,:,k].shape, order='F')
            
            if Do_Alpha_CX_Test:
                Alpha_CX_Test = sigmav_cx_h0(Ti_mu, ErelH_P) / Vth
                for k in range(nx):
                    Alpha_CX_Test[:,:,k] = Alpha_CX_Test[:,:,k]*ni[k]
                print('Compare alpha_cx and alpha_cx_test')
        # print("Alpha_CX", Alpha_CX.T)
        # print(Alpha_CX.shape)
        # input()
            

    #	Compute Alpha_H_H2 for inputted fH, if it is needed and has not
    #		already been computed with the present input parameters

    if Do_Alpha_H_H2 == 1:
        if debrief > 1:
            print(prompt+'Computing Alpha_H_H2')
        Alpha_H_H2 = np.zeros((nvr,nvx,nx))
        for k in range(nx):
            Work[:] = fH2[:,:,k].reshape(Work.shape, order='F')
            Alpha_H_H2[:,:,k] = (SIG_H_H2 @ Work).reshape(Alpha_H_H2[:,:,k].shape, order='F')
        # print("Alpha_H_H2", Alpha_H_H2[:,:,10].T)
        # print(Alpha_H_H2.shape)
        # input()

    #	Compute Alpha_H_P for present Ti and ni 
    #		if it is needed and has not already been computed with the present parameters
    if Do_Alpha_H_P == 1:
        if debrief > 1:
            print(prompt+'Computing Alpha_H_P')
        Alpha_H_P = np.zeros((nvr,nvx,nx))
        for k in range(nx):
            Work[:] = (fi_hat[:,:,k]*ni[k]).reshape(Work.shape, order='F')
            Alpha_H_P[:,:,k] = (SIG_H_P @ Work).reshape(Alpha_H_P[:,:,k].shape, order='F')
        # print("Alpha_H_P", Alpha_H_P.T)
        # print(Alpha_H_P.shape)
        # input()

    #	Compute nH
    for k in range(nx):
        nH[k] = np.sum(Vr2pidVr*(fH[:,:,k] @ dVx))
    
    if New_H_Seed:
        MH_H_sum = np.zeros((nvr,nvx,nx))
        Delta_nHs = 1

    #	Compute Side-Wall collision rate
    gamma_wall = np.zeros((nvr,nvx,nx))
    for k in range(nx):
        if PipeDia[k] > 0:
            for j in range(nvx):
                gamma_wall[:,j,k] = 2*vr / PipeDia[k]
    # print("nH", nH)
    # print("gamma_wall", gamma_wall.T)
    # input()

    do_fH_Iterate = True

    #	This is the entry point for fH iteration.
    #	Save 'seed' values for comparison later

    while do_fH_Iterate: #NOTE Alpha_CX done before here, but done inside iteration in kh2, does it change per iteration? Is this an error?
        do_fH_Iterate = False
        fHs = copy.copy(fH)
        nHs = copy.copy(nH)

        #	Compute Omega values if nH is non-zero
        ii = np.argwhere(nH <= 0)
        if ii.size <= 0: #NOTE Not Tested Yet, return on iteration

            #	Compute VxH
            if H_P_EL or H_H2_EL or H_H_EL:
                for k in range(nx):
                    VxH[k] = Vth*np.sum(Vr2pidVr*(fH[:,:,k] @ (vx*dVx))) / nH[k]
                # print("VxH", VxH)
                # input()

            #	Compute Omega_H_P for present fH and Alpha_H_P if H_P elastic collisions are included
            if H_P_EL:
                if debrief > 1:
                    print(prompt+'Computing Omega_H_P')
                for k in range(nx):
                    DeltaVx = (VxH[k] - vxi[k]) / Vth
                    MagDeltaVx = np.maximum(abs(DeltaVx), DeltaVx_tol)
                    DeltaVx = sign(DeltaVx)*MagDeltaVx
                    Omega_H_P[k] = np.sum(Vr2pidVr*((Alpha_H_P[:,:,k]*fH[:,:,k]) @ dVx)) / (nH[k]*DeltaVx)
                Omega_H_P = np.maximum(Omega_H_P, 0)
                # print("Omega_H_P", Omega_H_P)
                # input()

            #	Compute Omega_H_H2 for present fH and Alpha_H_H2 if H_H2 elastic collisions are included

            if H_H2_EL:
                if debrief > 1:
                    print(prompt+'Computing Omega_H_H2')
                for k in range(nx):
                    DeltaVx = (VxH[k] - vxH2[k]) / Vth
                    MagDeltaVx = np.maximum(abs(DeltaVx), DeltaVx_tol)
                    DeltaVx = sign(DeltaVx)*MagDeltaVx
                    # print("Mag", fH[:,:,k].T)
                    # input()
                    Omega_H_H2[k] = np.sum(Vr2pidVr*((Alpha_H_H2[:,:,k]*fH[:,:,k]) @ dVx)) / (nH[k]*DeltaVx)
                Omega_H_H2 = np.maximum(Omega_H_H2, 0)
                # print("Omega_H_H2", Omega_H_H2)
                # input()

            #	Compute Omega_H_H for present fH if H_H elastic collisions are included

            if H_H_EL:
                if debrief > 1:
                    print(prompt+'Computing Omega_H_H')
                if np.sum(MH_H_sum) <= 0:
                    for k in range(nx):
                        for i in range(nvr):
                            vr2_2vx_ran2[i,:] = vr[i]**2 - 2*((vx - (VxH[k]/Vth))**2)
                        Wperp_paraH[k] = np.sum(Vr2pidVr*((vr2_2vx_ran2*fH[:,:,k]) @ dVx)) / nH[k]
                else:
                    for k in range(nx):
                        M_fH = MH_H_sum[:,:,k] - fH[:,:,k]
                        Wperp_paraH[k] = -np.sum(Vr2pidVr*((vr2_2vx2_2D*M_fH) @ dVx)) / nH[k]
                for k in range(nx):
                    Work[:] = fH[:,:,k].reshape(Work.shape, order='F')
                    Alpha_H_H[:] = (SIG_H_H @ Work).reshape(Alpha_H_H.shape, order='F')
                    Wpp = Wperp_paraH[k]
                    MagWpp = np.maximum(np.abs(Wpp), Wpp_tol)
                    Wpp = sign(Wpp)*MagWpp
                    Omega_H_H[k] = np.sum(Vr2pidVr*((Alpha_H_H*Work.reshape(Alpha_H_H.shape, order='F')) @ dVx)) / (nH[k]*Wpp)
                Omega_H_H = np.maximum(Omega_H_H, 0)
                # print("Omega_H_H", Omega_H_H)
                # input()

        #	Total Elastic scattering frequency

        Omega_EL = Omega_H_P + Omega_H_H2 + Omega_H_H
        # print("Omega_EL", Omega_EL.T)
        # input()

        #	Total collision frequency

        alpha_c = np.zeros((nvr,nvx,nx))
        if H_P_CX:
            for k in range(nx):
                alpha_c[:,:,k] = Alpha_CX[:,:,k] + alpha_ion[k] + Omega_EL[k] + gamma_wall[:,:,k]
        else:
            for k in range(nx):
                alpha_c[:,:,k] = alpha_ion[k] + Omega_EL[k] + gamma_wall[:,:,k]
        # print("alpha_c", alpha_c.T)
        # input()

        #	Test x grid spacing based on Eq.(27) in notes
        if debrief>1:
            print(prompt+'Testing x grid spacing')
        Max_dx = np.full(nx, 1e32)
        for k in range(nx):
            for j in range(i_p[0][0], nvx): # changed ip to i_p
                Max_dx[k] = np.minimum(Max_dx[k], min(2*vx[j] / alpha_c[:,j,k]))

        dx = np.roll(x,-1) - x
        Max_dxL = Max_dx[0:nx-1]
        Max_dxR = Max_dx[1:nx]
        Max_dx = np.minimum(Max_dxL, Max_dxR)
        ilarge = np.argwhere(Max_dx < dx[0:nx-1])

        if ilarge.size > 0:
            print(prompt+'x mesh spacing is too large!') #NOTE Check Formatting
            debug = 1
            out = ""
            jj = 0

            #	Not sure the output is formatted correctly

            print(' \t    x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)\t   x(k+1)-x(k)   Max_dx(k)')
            for ii in range(ilarge.size):
                jj += 1
                out += ((str(ilarge[ii])+' \t')[:8]+(str(x[ilarge[ii]+1]-x[ilarge[ii]])+'        ')[:6]+'        '+str(Max_dx[ilarge[ii]])[:4]+'\t')
                if jj>4:
                    print(out)
                    jj = 0
                    out = "\t"
            if jj>0:
                print(out)
            error = 1
            raise Exception("x mesh spacing is too large")

        #	Define parameters Ak, Bk, Ck, Dk, Fk, Gk

        Ak = np.zeros((nvr,nvx,nx))
        Bk = np.zeros((nvr,nvx,nx))
        Ck = np.zeros((nvr,nvx,nx))
        Dk = np.zeros((nvr,nvx,nx))
        Fk = np.zeros((nvr,nvx,nx))
        Gk = np.zeros((nvr,nvx,nx))

        for k in range(0, nx-1):
            for j in range(i_p[0][0], nvx): # double check some of the ranges in for statements I might have some typos
                denom = 2*vx[j] + (x[k+1] - x[k])*alpha_c[:,j,k+1]
                Ak[:,j,k] = (2*vx[j] - (x[k+1] - x[k])*alpha_c[:,j,k]) / denom
                Bk[:,j,k] = (x[k+1] - x[k]) / denom
                Fk[:,j,k] = (x[k+1] - x[k])*(Sn[:,j,k+1]+Sn[:,j,k]) / denom
        for k in range(1, nx):
            for j in range(0, i_p[0][0]):
                denom = -2*vx[j] + (x[k] - x[k-1])*alpha_c[:,j,k-1]
                Ck[:,j,k] = (-2*vx[j] - (x[k] - x[k -1])*alpha_c[:,j,k]) / denom
                Dk[:,j,k] = (x[k] - x[k-1]) / denom
                Gk[:,j,k] = (x[k] - x[k-1])*(Sn[:,j,k]+Sn[:,j,k-1]) / denom
        # print("Ak", Ak.T)
        # print("Bk", Bk.T)
        # print("Ck", Ck.T)
        # print("Dk", Dk.T)
        # print("Fk", Fk.T)
        # print("Gk", Gk.T)
        # input()
                        
        #	Compute first-flight (0th generation) neutral distribution function
        Beta_CX_sum = np.zeros((nvr,nvx,nx))
        MH_P_sum = np.zeros((nvr,nvx,nx))
        MH_H2_sum = np.zeros((nvr,nvx,nx))
        MH_H_sum = np.zeros((nvr,nvx,nx))
        igen = 0
        if debrief > 0:
            print(prompt+'Computing atomic neutral generation#'+sval(igen))
        fHG[:,i_p,0] = fH[:,i_p,0]
        for k in range(nx-1):
            fHG[:,i_p,k+1] = fHG[:,i_p,k]*Ak[:,i_p,k] + Fk[:,i_p,k]
        for k in range(nx-1,0,-1):
            fHG[:,i_n,k-1] = fHG[:,i_n,k]*Ck[:,i_n,k] + Gk[:,i_n,k]
        # print("fHG", fHG.T)
        # input()
                
        #	Compute first-flight neutral density profile
        for k in range(nx):
            NHG[k,igen] = np.sum(Vr2pidVr*(fHG[:,:,k] @ dVx))
        # print("NHG", NHG.T)
        # input()

        # NOTE Add plotting once program is working
        # if plot >1:
        #     pass

        #	Set total atomic neutral distribution function to first flight generation

        fH = copy.copy(fHG)
        nH = NHG[:,0]

# next_generation #########################################################################################################################################################################
        while True:
            #print('check next_generation')
            if igen+1 > Max_Gen or fH_generations == 0:
                if debrief > 0:
                    print(prompt+'Completed '+sval(Max_Gen)+' generations. Returning present solution...')
                break
            igen += 1
            if debrief > 0:
                print(prompt+'Computing atomic neutral generation#'+sval(igen))

            #	Compute Beta_CX from previous generation

            Beta_CX = np.zeros((nvr,nvx,nx))
            if H_P_CX:
                if debrief>1:
                    print(prompt+'Computing Beta_CX')

                if Simple_CX:
                    #	Option (B): Compute charge exchange source with assumption that CX source 
                    #		neutrals have ion distribution function
                    for k in range(nx):
                        Beta_CX[:,:,k] = fi_hat[:,:,k]*np.sum(Vr2pidVr*((Alpha_CX[:,:,k]*fHG[:,:,k]) @ dVx))
                else:
                    #	Option (A): Compute charge exchange source using fH and vr x sigma x v_v at 
                    #		each velocity mesh point
                    for k in range(nx):
                        Work[:] = fHG[:,:,k]
                        Beta_CX[:,:,k] = ni[k]*fi_hat[:,:,k]*(SIG_CX @ Work)

                #	Sum charge exchange source over all generations
                Beta_CX_sum += Beta_CX
            # print("Beta_CX_sum", Beta_CX_sum.T)
            # input()

            #	Compute MH from previous generation
            MH_H = np.zeros((nvr,nvx,nx))
            MH_P = np.zeros((nvr,nvx,nx))
            MH_H2 = np.zeros((nvr,nvx,nx))
            OmegaM = np.zeros((nvr,nvx,nx))
            if H_H_EL or H_P_EL or H_H2_EL:

                #	Compute VxHG, THG
                for k in range(0, nx):
                    VxHG[k] = Vth*np.sum(Vr2pidVr*(fHG[:,:,k] @ (vx*dVx))) / NHG[k,igen-1]
                    for i in range(0, nvr):
                        vr2vx2_ran2[i,:] = vr[i]**2 + (vx - VxHG[k]/Vth)**2
                    # if igen >= 2:
                    #     print("vr2", vr2vx2_ran2)
                    #     input()
                    #     print("fHG", fHG[:,:,k].T)
                    #     input()
                    #     print("NHG", NHG[k,igen-1])
                    #     input()
                    THG[k] = (mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran2*fHG[:,:,k]) @ dVx)) / (3*CONST.Q*NHG[k,igen-1])
                # print("THG", THG)
                # input()

                if H_H_EL:
                    if debrief > 1:
                        print(prompt+'Computing MH_H')

                    #	Compute MH_H 
                    vx_shift = VxHG
                    Tmaxwell = THG
                    mol = 1
                    Maxwell = create_shifted_maxwellian_include(vr, vx, Tnorm, vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2_2D,jpa,jpb,jna,jnb)
                    for k in range(nx):
                        MH_H[:,:,k] = Maxwell[:,:,k]*NHG[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_H[k]*MH_H[:,:,k]
                    MH_H_sum += MH_H
                    # print("MH_H_sum", MH_H_sum.T)
                    # print("OmegaM", OmegaM.T)
                    # input()

                if H_P_EL:
                    if debrief>1:
                        print(prompt+'Computing MH_P')

                    #	Compute MH_P 
                    vx_shift = (VxHG+vxi)/2
                    Tmaxwell = THG + (2/4)*(Ti - THG + mu*CONST.H_MASS*((vxi - VxHG)**2) / (6*CONST.Q))
                    mol = 1
                    Maxwell = create_shifted_maxwellian_include(vr, vx, Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2_2D,jpa,jpb,jna,jnb)
                    for k in range(nx):
                        MH_P[:,:,k] = Maxwell[:,:,k]*NHG[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_P[k]*MH_P[:,:,k]
                    MH_P_sum += MH_P
                    # print("MH_P_sum", MH_P_sum.T)
                    # print("OmegaM", OmegaM.T)
                    # input()

                if H_H2_EL:
                    if debrief>1:
                        print(prompt+'Computing MH_H2')

                    #	Compute MH_H2
                    vx_shift = (VxHG + 2*vxH2)/3
                    Tmaxwell = THG + (4./9.)*(TH2 - THG + 2*mu*CONST.H_MASS*((vxH2 - VxHG)**2) / (6*CONST.Q))
                    mol = 1
                    Maxwell = create_shifted_maxwellian_include(vr, vx, Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,Vr2Vx2_2D,jpa,jpb,jna,jnb)
                    for k in range(nx):
                        MH_H2[:,:,k] = Maxwell[:,:,k]*NHG[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_H2[k]*MH_H2[:,:,k]
                    MH_H2_sum += MH_H2
                    # print("MH_H2_sum", MH_H2_sum.T)
                    # print("OmegaM", OmegaM.T)
                    # input()

            #	Compute next generation atomic distribution

            fHG[:] = 0
            # print("FHG", fHG)
            for k in range(0, nx-1):
                fHG[:,i_p,k+1] = Ak[:,i_p,k]*fHG[:,i_p,k] + Bk[:,i_p,k]*(Beta_CX[:,i_p,k+1] + OmegaM[:,i_p,k+1] + Beta_CX[:,i_p,k] + OmegaM[:,i_p,k])
            for k in range(nx-1, 0, -1):
                fHG[:,i_n,k-1] = Ck[:,i_n,k]*fHG[:,i_n,k] + Dk[:,i_n,k]*(Beta_CX[:,i_n,k-1] + OmegaM[:,i_n,k-1] + Beta_CX[:,i_n,k] + OmegaM[:,i_n,k])
                # print("FHG1", fHG[:,i_n,k-1].T)
                # input()
            for k in range(0, nx):
                NHG[k,igen] = np.sum(Vr2pidVr*(fHG[:,:,k] @ dVx))
            # print("fHG", fHG.T)
            # print("NHG", NHG.T)
            # input()

            # NOTE Add plotting once program is working
            # if plot > 1:
            #     pass

            #	Add result to total neutral distribution function
            fH += fHG
            nH += NHG[:,igen]

            #	Compute 'generation error': Delta_nHG=max(NHG(*,igen)/max(nH))
            #		and decide if another generation should be computed
            Delta_nHG = max(NHG[:,igen]/max(nH))
            # print("Delta_nHG", Delta_nHG)
            # input()
            if (Delta_nHG < truncate) or (fH_iterate and (Delta_nHG < 0.003*Delta_nHs)):
                #	If fH 'seed' is being iterated, then do another generation until the 'generation error'
                #		is less than 0.003 times the 'seed error' or is less than TRUNCATE
                break

# fH2_done #########################################################################################################################################################################
        
        # NOTE Add plotting once program is working
        # if plot>0:
        #     pass

        #	Compute H density profile
        for k in range(0, nx):
            nH[k] = np.sum(Vr2pidVr*(fH[:,:,k] @ dVx))
        # print("nH", nH)
        # input()

        if fH_iterate:

            #	Compute 'seed error': Delta_nHs=(|nHs-nH|)/max(nH) 
            #		If Delta_nHs is greater than 10*truncate then iterate fH
            Delta_nHs = np.max(np.abs(nHs - nH))/np.max(nH)
            # print("Delta_nHs", Delta_nHs)
            # input()
            if Delta_nHs > 10*truncate:
                do_fH_Iterate = True 

    #	Update Beta_CX_sum using last generation
    if H_P_CX:
        if debrief > 1:
            print(prompt, 'Computing Beta_CX')
        if Simple_CX:
            # Option (B): Compute charge exchange source with assumption that CX source neutrals have
            # ion distribution function
            for k in range(0, nx):
                Beta_CX[:,:,k] = fi_hat[:,:,k]*np.sum(Vr2pidVr*(Alpha_CX[:,:,k]*fHG[:,:,k] @ dVx))
                # print((Alpha_CX[:,:,k]*fH2G[:,:,k] @ dVx).T)
                # input()
        else:
            # Option (A): Compute charge exchange source using fH2 and vr x sigma x v_v at each velocity mesh point
            for k in range(0, nx):
                Work[:] = fHG[:,:,k]
                Beta_CX[:,:,k] = ni[k]*fi_hat[:,:,k]*(SIG_CX @ Work)
        Beta_CX_sum = Beta_CX_sum + Beta_CX
    # print("Beta_CX", Beta_CX.T)
    # print("Beta_CX_sum", Beta_CX_sum.T)
    # input()
            
    #	Update MH_*_sum using last generation
    MH_H2 = np.zeros((nvr,nvx,nx))
    MH_P = np.zeros((nvr,nvx,nx))
    MH_H = np.zeros((nvr,nvx,nx))
    OmegaM = np.zeros((nvr,nvx,nx))
    if H_H_EL or H_P_EL or H_H2_EL: 
        # Compute VxH2G, TH2G
        for k in range(0, nx):
            VxHG[k] = Vth*np.sum(Vr2pidVr*(fHG[:,:,k] @ (vx*dVx))) / NHG[k,igen]
            for i in range(0, nvr):
                vr2vx2_ran2[i,:] = vr[i]**2 + (vx - VxHG[k]/Vth)**2
            THG[k] = (mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran2*fHG[:,:,k]) @ dVx)) / (3*CONST.Q*NHG[k,igen])
        # print("THG", THG)
        # input()

        if H_H_EL:
            if debrief > 1: 
                print(prompt, 'Computing MH_H')
            # Compute MH_H
            vx_shift = VxHG
            Tmaxwell = np.copy(THG)
            mol = 1
            Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                      nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
            for k in range(0, nx):
                MH_H[:,:,k] = Maxwell[:,:,k]*NHG[k,igen]
                OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_H[k]*MH_H[:,:,k]
            MH_H_sum = MH_H_sum + MH_H
            # print("MH_H_sum", MH_H_sum.T)
            # input()

        if H_P_EL:
            if debrief > 1:
                print(prompt, 'Computing MH_P')
            # Compute MH_P
            vx_shift = (VxHG + vxi)/2
            Tmaxwell = THG + (2/4)*(Ti - THG + mu*CONST.H_MASS*((vxi - VxHG)**2) / (6*CONST.Q))
            mol = 1
            Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                    nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                    Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
            for k in range(0, nx):
                MH_P[:,:,k] = Maxwell[:,:,k]*NHG[k,igen]
                OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_P[k]*MH_P[:,:,k]
            MH_P_sum = MH_P_sum + MH_P
            # print("MH_P_sum", MH_P_sum.T)
            # input()

            if H_H2_EL: #NOTE Not Tested Yet
                if debrief > 1:
                    print(prompt, 'Computing MH_H2')
                # Compute MH_H
                vx_shift = (VxHG + 2*vxH2)/3
                Tmaxwell = THG + (4/9)*(TH2 - THG + 2*mu*CONST.H_MASS*((vxH2 - VxHG)**2) / (6*CONST.Q))
                mol = 1
                Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,Shifted_Maxwellian_Debug,mu,mol,
                                        nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                        Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
                for k in range(0, nx):
                    MH_H2[:,:,k] = Maxwell[:,:,k]*NHG[k,igen]
                    OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H_H2[k]*MH_H2[:,:,k]
                MH_H2_sum = MH_H2_sum + MH_H2
                # print("MH_H2_sum", MH_H2_sum.T)
                # input()

    #	Compute remaining moments

    #NOTE In kinetic_h2, these are calculated in the iteration, does that matter?
    #	GammaxH - particle flux in x direction
    for k in range(0, nx):
            GammaxH[k] = Vth*np.sum(Vr2pidVr*(fH[:,:,k] @ (vx*dVx)))
    # print("GammaxH", GammaxH)
    # input()

    #	VxH - x velocity
    VxH = GammaxH / nH
    _VxH = VxH / Vth 
    # print("VxH", VxH)
    # print("_VxH", _VxH)
    # input()

    #	magnitude of random velocity at each mesh point
    vr2vx2_ran = np.zeros((nvr,nvx,nx))
    for i in range(0, nvr):
        for k in range(0, nx):
            vr2vx2_ran[i,:,k] = vr[i]**2 + (vx - _VxH[k])**2
    # print("vr2vx2_ran", vr2vx2_ran.T)
    # input()

    #	pH - pressure 
    for k in range(nx):
        pH[k] = ((mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*fH[:,:,k]) @ dVx))) / (3*CONST.Q)
    # print("pH", pH)
    # input()

    #	TH - temperature
    TH = pH/nH
    # print("TH", TH)
    # input()

    #	piH_xx
    for k in range(nx):
        piH_xx[k] = (((mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*(fH[:,:,k] @ (dVx*(vx - _VxH[k])**2)))) / CONST.Q) - pH[k]
    #	piH_yy
    for k in range(nx):
        piH_yy[k] = (((mu*CONST.H_MASS)*Vth2*0.5*np.sum((Vr2pidVr*(vr**2))*(fH[:,:,k] @ dVx))) / CONST.Q) - pH[k]
    #	piH_zz
    piH_zz = copy.copy(piH_yy)
    #	qxH
    for k in range(nx):
        qxH[k] = 0.5*(mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*fH[:,:,k]) @ (dVx*(vx - _VxH[k]))))
    # print("piH2_xx", piH_xx)
    # print("piH2_yy", piH_yy)
    # print("qxH", qxH)
    # input()

    #	C = RHS of Boltzman equation for total fH

    for k in range(nx):
        C = Vth*(Sn[:,:,k] + Beta_CX_sum[:,:,k] - alpha_c[:,:,k]*fH[:,:,k] + \
                 Omega_H_P[k]*MH_P_sum[:,:,k] + Omega_H_H2[k]*MH_H2_sum[:,:,k] + Omega_H_H[k]*MH_H_sum[:,:,k])
        # print("C", C.T)
        # input()
        QH[k] = 0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*C) @ dVx))
        RxH[k] = (mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(C @ (dVx*(vx - _VxH[k]))))
        NetHSource[k] = np.sum(Vr2pidVr*(C @ dVx))
        Sion[k] = Vth*nH[k]*alpha_ion[k]
        SourceH[k] = np.sum(Vr2pidVr*(fSH[:,:,k] @ dVx))
        WallH[k] = np.sum(Vr2pidVr*((gamma_wall[:,:,k]*fH[:,:,k]) @ dVx))

        if Recomb:
            SRecomb[k] = Vth*ni[k]*Rec[k]
        else:
            SRecomb[k] = 0

        if H_P_CX:
            CCX = Vth*(Beta_CX_sum[:,:,k] - Alpha_CX[:,:,k]*fH[:,:,k])
            RxHCX[k] = (mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CCX @ (dVx*(vx - _VxH[k]))))
            EHCX[k] = 0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CCX) @ dVx))

        if H_H2_EL:
            CH_H2 = Vth*Omega_H_H2[k]*(MH_H2_sum[:,:,k] - fH[:,:,k])
            RxH2_H[k] = (mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CH_H2 @ (dVx*(vx - _VxH[k]))))
            EH2_H[k] = 0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CH_H2) @ dVx))

        if H_P_EL:
            CH_P = Vth*Omega_H_P[k]*(MH_P_sum[:,:,k] - fH[:,:,k])
            RxP_H[k] = (mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CH_P @ (dVx*(vx - _VxH[k]))))
            EP_H[k] = 0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CH_P) @ dVx))

        CW_H = -Vth*(gamma_wall[:,:,k]*fH[:,:,k])
        RxW_H[k] = (mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CW_H @ (dVx*(vx - _VxH[k]))))
        EW_H[k] = 0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CW_H) @ dVx))
        
        if H_H_EL:
            CH_H = Vth*Omega_H_H[k]*(MH_H_sum[:,:,k] - fH[:,:,k])
            for i in range(0, nvr):
                vr2_2vx_ran2[i,:] = vr[i]**2 - 2*((vx - _VxH[k])**2)
            Epara_PerpH_H[k] = -0.5*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2_2vx_ran2*CH_H) @ dVx))
    # print("QH", QH)
    # print("RxH", RxH)
    # print("NetHSource", NetHSource)
    # print("Sion", Sion)
    # print("SourceH", SourceH)
    # print("WallH", WallH)
    # print("SRecomb", SRecomb)
    # input()
    # print("RxHCX", RxHCX)
    # print("EHCX", EHCX)
    # input()
    # print("RxH2_H", RxH2_H)
    # print("EH2_H", EH2_H)
    # input()
    # print("RxP_H", RxP_H)
    # print("EP_H", EP_H)
    # input()
    # print("RxW_H", RxW_H)
    # print("EW_H", EW_H)
    # input()
    # print("Epara_PerpH_H", Epara_PerpH_H)
    # input()

    #	qxH_total
    qxH_total = (0.5*nH*(mu*CONST.H_MASS)*VxH*VxH + 2.5*pH*CONST.Q)*VxH + CONST.Q*piH_xx*VxH + qxH
    # print("qxH_total", qxH_total)
    # input()

    #	QH_total
    QH_total = QH + RxH*VxH + 0.5*(mu*CONST.H_MASS)*NetHSource*VxH*VxH
    # print("QH_total", QH_total)
    # input()

    #	Albedo
    gammax_plus = Vth*np.sum(Vr2pidVr*(fH[:,i_p[:,0],0] @ (vx[i_p[:,0]]*dVx[i_p[:,0]]))) #NOTE Had to reference i_p and i_n in a weird way, fix how they are called in the first place
    gammax_minus = Vth*np.sum(Vr2pidVr*(fH[:,i_n[:,0],0] @ (vx[i_n[:,0]]*dVx[i_n[:,0]]))) #This is awful and should not be allowed to remain
    if np.abs(gammax_plus) > 0:
        AlbedoH = -gammax_minus/gammax_plus
    # print("AlbedoH", AlbedoH)
    # input()

    #	Compute Mesh Errors

    mesh_error = np.zeros((nvr,nvx,nx))
    max_mesh_error = 0.0
    min_mesh_error = 0.0
    mtest = 5
    moment_error = np.zeros((nx,mtest))
    max_moment_error = np.zeros(mtest)
    C_error = np.zeros(nx)
    CX_error = np.zeros(nx)
    H_H_error = np.zeros((nx, 3))
    H_H2_error = np.zeros((nx, 3))
    H_P_error = np.zeros((nx, 3))
    max_H_H_error = np.zeros(3)
    max_H_H2_error = np.zeros(3)
    max_H_P_error = np.zeros(3)

    if Compute_Errors:
        if debrief > 1:
            print(prompt+'Computing Collision Operator, Mesh, and Moment Normalized Errors')
        NetHSource2 = SourceH + SRecomb - Sion - WallH
        for k in range(nx):
            C_error[k] = abs(NetHSource[k] - NetHSource2[k]) / max(abs(np.array([NetHSource[k], NetHSource2[k]])))
        # print("C_error", C_error)
        # input()

        #	Test conservation of particles for charge exchange operator
        if H_P_CX:
            for k in range(nx):
                CX_A = np.sum(Vr2pidVr*((Alpha_CX[:,:,k]*fH[:,:,k]) @ dVx))
                CX_B = np.sum(Vr2pidVr*(Beta_CX_sum[:,:,k] @ dVx))
                CX_error[k] = np.abs(CX_A - CX_B) / np.max(np.abs(np.array([CX_A, CX_B])))
            # print("CX_error", CX_error)
            # input()

        #	Test conservation of particles, x momentum, and total energy of elastic collision operators
        for m in range(0, 3):
            for k in range(0, nx):
                if m < 2:
                    TfH = np.sum(Vr2pidVr*(fH[:,:,k] @ (dVx*(vx**m))))
                else:
                    TfH = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*fH[:,:,k]) @ dVx))

                if H_H_EL:
                    if m < 2:
                        TH_H = np.sum(Vr2pidVr*(MH_H_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH_H = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH_H_sum[:,:,k]) @ dVx))
                    H_H_error[k,m] = np.abs(TfH - TH_H) / np.max(np.abs(np.array([TfH, TH_H])))
                
                if H_H2_EL:
                    if m < 2:
                        TH_H2 = np.sum(Vr2pidVr*(MH_H2_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH_H2 = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH_H2_sum[:,:,k]) @ dVx))
                    H_H2_error[k,m] = np.abs(TfH - TH_H2) / np.max(np.abs(np.array([TfH, TH_H2])))

                if H_P_EL:
                    if m < 2:
                        TH_P = np.sum(Vr2pidVr*(MH_P_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH_P = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH_P_sum[:,:,k]) @ dVx))
                    H_P_error[k,m] = np.abs(TfH - TH_P) / np.max(np.abs(np.array([TfH, TH_P])))

            max_H_H_error[m] = np.max(H_H_error[:,m])
            max_H_H2_error[m] = np.max(H_H2_error[:,m])
            max_H_P_error[m] = np.max(H_P_error[:,m])
            # print("max_H_H_error", max_H_H_error)
            # print("max_H_H2_error", max_H_H2_error)
            # print("max_H_P_error", max_H_P_error)
            # input()

        if CI_Test:
            #	Compute Momentum transfer rate via full collision integrals for charge exchange and 
            #		mixed elastic scattering.
            #		Then compute error between this and actual momentum transfer 
            #		resulting from CX and BKG (elastic) models.

            if H_P_CX: # P -> H charge exchange momentum transfer via full collision integral
                print(prompt, 'Computing P -> H2 Charge Exchange Momentum Transfer')
                _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
                _Sig[:] = (v_v*sigma_cx_h0(v_v2*(0.5*CONST.H_MASS*Vth2 / CONST.Q))).reshape(_Sig.shape, order='F')
                SIG_VX_CX = np.zeros((nvr*nvx,nvr*nvx))
                SIG_VX_CX[:] = (Vr2pidVrdVx*vx_vx*((_Sig @ dTheta).reshape(vx_vx.shape, order='F'))).reshape(SIG_VX_CX.shape, order='F')
                alpha_vx_cx = np.zeros((nvr,nvx,nx))

                for k in range(0, nx):
                    Work[:] = (fi_hat[:,:,k]*ni[k]).reshape(Work.shape, order='F')
                    alpha_vx_cx[:,:,k] = (SIG_VX_CX @ Work).reshape(alpha_vx_cx[:,:,k].shape, order='F')

                for k in range(0, nx):
                    RxCI_CX[k] = -(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((alpha_vx_cx[:,:,k]*fH[:,:,k]) @ dVx))

                norm = np.max(np.abs(np.array([RxHCX, RxCI_CX])))
                for k in range(0, nx):
                    CI_CX_error[k] = np.abs(RxHCX[k] - RxCI_CX[k]) / norm

                print(prompt,'Maximum normalized momentum transfer error in CX collision operator: ', sval(np.max(CI_CX_error)))

            if H_P_EL: # P -> H momentum transfer via full collision integral
                for k in range(0, nx):
                    RxCI_P_H[k] = -(1/2)*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((Alpha_H_P[:,:,k]*fH[:,:,k]) @ dVx))

                norm = np.max(np.abs(np.array([RxP_H, RxCI_P_H])))
                for k in range(0, nx):
                    CI_P_H_error[k] = np.abs(RxP_H[k] - RxCI_P_H[k]) / norm 

                print(prompt, 'Maximum normalized momentum transfer error in P -> H elastic BKG collision operator: ', sval(np.max(CI_P_H_error)))

            if H_H2_EL: # H2 -> H momentum transfer via full collision integral
                for k in range(0, nx):
                    RxCI_H2_H[k] = -(2/3)*(mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((Alpha_H_H2[:,:,k]*fH[:,:,k]) @ dVx))
                
                norm = np.max(np.abs(np.array([RxH2_H, RxCI_H2_H])))
                for k in range(0, nx):
                    CI_H2_H_error[k] = np.abs(RxH2_H[k] - RxCI_H2_H[k])/norm
                
                print(prompt, 'Maximum normalized momentum transfer error in H2 -> H elastic BKG collision operator: ', sval(np.max(CI_H2_H_error)))

            if H_H_EL: # H -> H perp/parallel energy transfer via full collision integral
                for k in range(0, nx):
                    Work[:] = fH[:,:,k].reshape(Work.shape, order='F')
                    Alpha_H_H[:] = (SIG_H_H @ Work).reshape(Alpha_H_H.shape, order='F')
                    Epara_Perp_CI[k] = 0.5*(mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((Alpha_H_H*fH[:,:,k]) @ dVx)) 
                
                norm = np.max(np.abs(np.array([Epara_PerpH_H, Epara_Perp_CI])))
                for k in range(0, nx):
                    CI_H_H_error[k] = np.abs(Epara_PerpH_H[k] - Epara_Perp_CI[k]) / norm 
                
                print(prompt, 'Maximum normalized perp/parallel energy transfer error in H -> H elastic BKG collision operator: ', sval(np.max(CI_H_H_error)))

        #	Mesh Point Error based on fH satisfying Boltzmann equation

        T1 = np.zeros((nvr,nvx,nx))
        T2 = np.zeros((nvr,nvx,nx))
        T3 = np.zeros((nvr,nvx,nx))
        T4 = np.zeros((nvr,nvx,nx))
        T5 = np.zeros((nvr,nvx,nx))
        for k in range(0, nx-1):
            for j in range(0, nvx):
                T1[:,j,k] = 2*vx[j]*(fH[:,j,k+1] - fH[:,j,k]) / (x[k+1] - x[k]) 
            T2[:,:,k] = (Sn[:,:,k+1] + Sn[:,:,k])
            T3[:,:,k] = Beta_CX_sum[:,:,k+1] + Beta_CX_sum[:,:,k]
            T4[:,:,k] = alpha_c[:,:,k+1]*fH[:,:,k+1] + alpha_c[:,:,k]*fH[:,:,k]
            T5[:,:,k] = Omega_H_P[k+1]*MH_P_sum[:,:,k+1] + Omega_H_H2[k+1]*MH_H2_sum[:,:,k+1] + Omega_H_H[k+1]*MH_H_sum[:,:,k+1] + \
                    Omega_H_P[k]*MH_P_sum[:,:,k] + Omega_H_H2[k]*MH_H2_sum[:,:,k] + Omega_H_H[k]*MH_H_sum[:,:,k]
            mesh_error[:,:,k] = np.abs(T1[:,:,k] - T2[:,:,k] - T3[:,:,k] + T4[:,:,k] - T5[:,:,k]) / \
                                np.max(np.abs(np.array([T1[:,:,k], T2[:,:,k], T3[:,:,k], T4[:,:,k], T5[:,:,k]])))
        ave_mesh_error = np.sum(mesh_error) / np.size(mesh_error)
        max_mesh_error = np.max(mesh_error)
        min_mesh_error = np.min(mesh_error[:,:,0:nx-1])

        #	Moment Error
        for m in range(0, mtest):
            for k in range(0, nx-1):
                MT1 = np.sum(Vr2pidVr*(T1[:,:,k] @ (dVx*(vx**m))))
                MT2 = np.sum(Vr2pidVr*(T2[:,:,k] @ (dVx*(vx**m))))
                MT3 = np.sum(Vr2pidVr*(T3[:,:,k] @ (dVx*(vx**m))))
                MT4 = np.sum(Vr2pidVr*(T4[:,:,k] @ (dVx*(vx**m))))
                MT5 = np.sum(Vr2pidVr*(T5[:,:,k] @ (dVx*(vx**m))))
                #NOTE This is correct for the original code, but is it correct mathematically?
                moment_error[k,m] = np.abs(MT1 - MT2 - MT3 + MT4 - MT5) / np.max(np.abs(np.array([MT1, MT2, MT3, MT4, MT5])))
            max_moment_error[m] = np.max(moment_error[:,m])

        #	Compute error in qxH_total

        #		qxH_total2 total neutral heat flux profile (watts m^-2)
        #			This is the total heat flux transported by the neutrals
        #			computed in a different way from:

        #			qxH_total2(k)=vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#(Vx*dVx)))*0.5*(mu*mH)

        #			This should agree with qxH_total if the definitions of nH, pH, piH_xx,
        #			TH, VxH, and qxH are coded correctly.
        qxH_total2 = np.zeros(nx)
        for k in range(0, nx):
            qxH_total2[k] = 0.5*(mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*fH[:,:,k]) @ (vx*dVx)))
        qxH_total_error = np.abs(qxH_total - qxH_total2) / np.max(np.abs(np.array([qxH_total, qxH_total2])))

        #	Compute error in QH_total
        Q1 = np.zeros(nx)
        Q2 = np.zeros(nx)
        QH_total_error = np.zeros(nx)
        for k in range(0, nx-1):
            Q1[k] = (qxH_total[k+1] - qxH_total[k]) / (x[k+1] - x[k])
            Q2[k] = 0.5*(QH_total[k+1] + QH_total[k])
        QH_total_error = np.abs(Q1 - Q2) / np.max(np.abs(np.array([Q1, Q2])))

        if debrief > 0:
            print(prompt+'Maximum particle convervation error of total collision operator: '+sval(max(C_error)))
            print(prompt+'Maximum H_P_CX  particle convervation error: '+sval(max(CX_error)))
            print(prompt+'Maximum H_H_EL  particle conservation error: '+sval(max_H_H_error[0]))
            print(prompt+'Maximum H_H_EL  x-momentum conservation error: '+sval(max_H_H_error[1]))
            print(prompt+'Maximum H_H_EL  total energy conservation error: '+sval(max_H_H_error[2]))
            print(prompt+'Maximum H_H2_EL particle conservation error: '+sval(max_H_H2_error[0]))
            print(prompt+'Maximum H_P_EL  particle conservation error: '+sval(max_H_P_error[0]))
            print(prompt+'Average mesh_error = '+str(ave_mesh_error))
            print(prompt+'Maximum mesh_error = '+str(max_mesh_error))
            for m in range(5):
                print(prompt+'Maximum fH vx^'+sval(m)+' moment error: '+sval(max_moment_error[m]))
            print(prompt+'Maximum qxH_total error = '+str(max(qxH_total_error)))
            print(prompt+'Maximum QH_total error = '+str(max(QH_total_error)))
            if debug > 0:
                input()

    # NOTE Add plotting once program is working
    # if plot>1:	#	May add later
    #     pass

    #	lines 1665 - 1748 in original code are used for plotting
    #		This may be added later, but has been left out for now

    #	Save input parameters in kinetic_H_input common block

    KH_Common.Input.vx_s = vx
    KH_Common.Input.vr_s = vr
    KH_Common.Input.x_s = x
    KH_Common.Input.Tnorm_s = Tnorm
    KH_Common.Input.mu_s = mu
    KH_Common.Input.Ti_s = Ti
    KH_Common.Input.vxi_s = vxi
    KH_Common.Input.Te_s = Te
    KH_Common.Input.n_s = n
    KH_Common.Input.vxi_s = vxi
    KH_Common.Input.fHBC_s = fHBC
    KH_Common.Input.GammaxHBC_s = GammaxHBC
    KH_Common.Input.PipeDia_s = PipeDia
    KH_Common.Input.fH2_s = fH2
    KH_Common.Input.fSH_s = fSH
    KH_Common.Input.nHP_s = nHP
    KH_Common.Input.THP_s = THP
    KH_Common.Input.fH_s = fH
    KH_Common.Input.Simple_CX_s = Simple_CX
    KH_Common.Input.JH_s = JH
    KH_Common.Input.Collrad_s = Use_Collrad_Ionization
    KH_Common.Input.Recomb_s = Recomb
    KH_Common.Input.H_H_EL_s = H_H_EL
    KH_Common.Input.H_P_EL_s = H_P_EL
    KH_Common.Input.H_H2_EL_s = H_H2_EL
    KH_Common.Input.H_P_CX_s = H_P_CX

    #	Save input parameters in kinetic_H_internal common block

    KH_Common.Internal.vr2vx2 = vr2vx2
    KH_Common.Internal.vr2vx_vxi2 = vr2vx_vxi2
    KH_Common.Internal.fi_hat = fi_hat
    KH_Common.Internal.ErelH_P = ErelH_P
    KH_Common.Internal.Ti_mu = Ti_mu
    KH_Common.Internal.ni = ni
    KH_Common.Internal.sigv = sigv
    KH_Common.Internal.alpha_ion = alpha_ion
    KH_Common.Internal.v_v2 = v_v2
    KH_Common.Internal.v_v = v_v
    KH_Common.Internal.vr2_vx2 = vr2_vx2
    KH_Common.Internal.vx_vx = vx_vx
    KH_Common.Internal.Vr2pidVrdVx = Vr2pidVrdVx
    KH_Common.Internal.SIG_CX = SIG_CX
    KH_Common.Internal.SIG_H_H = SIG_H_H
    KH_Common.Internal.SIG_H_H2 = SIG_H_H2
    KH_Common.Internal.SIG_H_P = SIG_H_P
    KH_Common.Internal.Alpha_CX = Alpha_CX
    KH_Common.Internal.Alpha_H_H2 = Alpha_H_H2
    KH_Common.Internal.Alpha_H_P = Alpha_H_P
    KH_Common.Internal.MH_H_sum = MH_H_sum
    KH_Common.Internal.Delta_nHs = Delta_nHs
    KH_Common.Internal.Sn = Sn
    KH_Common.Internal.Rec = Rec

    #	Save input parameters in kinetic_H_Moments common block

    KH_Common.Moments.nH2 = nH2
    KH_Common.Moments.VxH2 = vxH2
    KH_Common.Moments.TH2 = TH2

    return fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,WallH,error
