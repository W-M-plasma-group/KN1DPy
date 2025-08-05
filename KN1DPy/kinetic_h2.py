import numpy as np
from scipy.ndimage import shift
import matplotlib.pyplot as plt
import copy

from .make_dvr_dvx import make_dvr_dvx
from .create_shifted_maxwellian import create_shifted_maxwellian
from .kinetic_mesh import kinetic_mesh

from .sigma.sigmav_ion_hh import sigmav_ion_hh
from .sigma.sigmav_h1s_h1s_hh import sigmav_h1s_h1s_hh
from .path_interp_2d import path_interp_2d
from .sigma.sigmav_h1s_h2s_hh import sigmav_h1s_h2s_hh
from .sigma.sigmav_p_h1s_hh import sigmav_p_h1s_hh
from .sigma.sigmav_h2p_h2s_hh import sigmav_h2p_h2s_hh
from .sigma.sigmav_h1s_hn3_hh import sigmav_h1s_hn3_hh
from .sigma.sigmav_p_h1s_hp import sigmav_p_h1s_hp
from .sigma.sigmav_p_hn2_hp import sigmav_p_hn2_hp
from .sigma.sigmav_p_p_hp import sigmav_p_p_hp
from .sigma.sigmav_h1s_hn_hp import sigmav_h1s_hn_hp
from .sigma.sigma_cx_hh import sigma_cx_hh
from .sigma.sigma_el_h_hh import sigma_el_h_hh
from .sigma.sigma_el_p_hh import sigma_el_p_hh
from .sigma.sigma_el_hh_hh import sigma_el_hh_hh
from .create_shifted_maxwellian_include import create_shifted_maxwellian_include
from .sigma.sigmav_cx_hh import sigmav_cx_hh

from .sign import sign #NOTE Replace sign with np.sign
from .sval import sval
from .locate import locate

from .common.Kinetic_H2 import Kinetic_H2_Common
from .common import constants as CONST

# This subroutine is part of the "KN1D" atomic and molecular neutal transport code.

#   This subroutine solves a 1-D spatial, 2-D velocity kinetic neutral transport
# problem for molecular hydrogen or deuterium (H2) by computing successive generations of
# charge exchange and elastic scattered neutrals. The routine handles electron-impact
# ionization and dissociation, molecular ion charge exchange, and elastic
# collisions with hydrogenic ions, neutral atoms, and molecules.

#   The positive vx half of the atomic neutral distribution function is inputted at x(0)
# (with arbitrary normalization). The desired flux on molecules entering the slab geometry at
# x(0) is specified. Background profiles of plasma ions, (e.g., Ti(x), Te(x), n(x), vxi(x),...)
# and atomic distribution function (fH) is inputted. (fH can be computed by procedure 
# "Kinetic_H.pro".) Optionally, a molecular source profile (SH2(x)) is also inputted. 
# The code returns the molecular hydrogen distribution function, fH2(vr,vx,x) for all 
# vx, vr, and x of the specified vr,vx,x grid. The atomic (H) and ionic (P) hydrogen 
# source profiles and the atomic source velocity distribution functions 
# resulting from Franck-Condon reaction product energies of H are also returned.

#   Since the problem involves only the x spatial dimension, all distribution functions
# are assumed to have rotational symmetry about the vx axis. Consequently, the distributions
# only depend on x, vx and vr where vr =sqrt(vy^2+vz^2)

#  History:

#    B. LaBombard   First coding based on Kinetic_Neutrals.pro 22-Dec-2000

#    For more information, see write-up: "A 1-D Space, 2-D Velocity, Kinetic
#    Neutral Transport Algorithm for Hydrogen Molecules in an Ionizing Plasma", B. LaBombard

# Note: Variable names contain characters to help designate species -
#       atomic neutral (H), molecular neutral (H2), molecular ion (HP), proton (i) or (P)

def kinetic_h2(mesh : kinetic_mesh, mu, vxi, fH2BC, GammaxH2BC, NuLoss, fH, SH2, fH2, nHP, THP, KH2 : Kinetic_H2_Common, 
               truncate = 1.0e-4, Simple_CX = 1, Max_Gen = 50,  Compute_H_Source = 0,
               No_Sawada = 0, H2_H2_EL = 0, H2_P_EL = 0, H2_H_EL = 0, H2_HP_CX = 0,
               ni_correct = 0, ESH = 0, Eaxis = 0, Compute_Errors = 0,  plot = 0, debug = 0,
               debrief = 0, pause = 0):
    
    #Print Inputs
    # print(mesh)
    # print("mu", mu)
    # print("vxi", vxi)
    # print("fH2BC", fH2BC)
    # print("GammaxH2BC", GammaxH2BC)
    # print("NuLoss", NuLoss)
    # print("fH", fH)
    # print("SH2", SH2)
    # print("fH2", fH2)
    # print("nHP", nHP)
    # print("THP", THP)
    # input()
    
    vx = mesh.vx
    vr = mesh.vr
    x = mesh.x
    Tnorm = mesh.Tnorm
    Ti = mesh.Ti
    Te = mesh.Te
    n = mesh.ne
    PipeDia = mesh.PipeDia
    
    # Kinetic_H2_Output common block
    piH2_xx = KH2.Output.piH2_xx
    piH2_yy = KH2.Output.piH2_yy
    piH2_zz = KH2.Output.piH2_zz
    RxH2CX = KH2.Output.RxH2CX
    RxH_H2 = KH2.Output.RxH_H2
    RxP_H2 = KH2.Output.RxP_H2
    RxW_H2 = KH2.Output.RxW_H2
    EH2CX = KH2.Output.EH2CX
    EH_H2 = KH2.Output.EH_H2
    EP_H2 = KH2.Output.EP_H2
    EW_H2 = KH2.Output.EW_H2
    Epara_PerpH2_H2 = KH2.Output.Epara_PerpH2_H2

    # Kinetic_H2_Errors common block
    Max_dx = KH2.Errors.Max_dx
    vbar_error = KH2.Errors.vbar_error
    mesh_error = KH2.Errors.mesh_error
    C_error = KH2.Errors.C_Error
    CX_error = KH2.Errors.CX_Error
    H2_H2_error = KH2.Errors.H2_H2_error
    qxH2_total_error = KH2.Errors.qxH2_total_error
    QH2_total_error = KH2.Errors.QH2_total_error

    #  Input:
    #		  vx(*)	- fltarr(nvx), normalized x velocity coordinate 
    #			  [negative values, positive values],
    #			  monotonically increasing. Note: a nonuniform mesh can be used.
    #			  Dimensional velocity (note: Vth is based on ATOM mass)
    #			  is v = Vth * vx where Vth=sqrt(2 k Tnorm/(mH*mu))
    #			  Note: nvx must be even and vx(*) symmetric about 
    #			  zero but not contain a zero element
    #		  vr(*)	- fltarr(nvr), normalized radial velocity coordinate 
    #			  [positive values], monotonically increasing. Note: a non-uniform mesh can be used.
    #			  Dimensional velocity is v = Vth * vr where Vth=sqrt(2 k Tnorm/(mH*mu)) 
    #			  Note: vr must not contain a zero element
    #		   x(*)	- fltarr(nx), spatial coordinate (meters), 
    #			  positive, monontonically increasing. Note: a non-uniform mesh can be used.
    #		  Tnorm	- Float, temperature corresponding to the thermal speed (see vx and vr above) (eV)
    #		     mu	- Float, 1=hydrogen, 2=deuterium
    #		     Ti	- fltarr(nx), Ion temperature profile (eV)
    #		     Te	- fltarr(nx), electron temperature profile (eV)
    #		      n	- fltarr(nx), electron density profile (m^-3)
    #		    vxi	- fltarr(nx), x-directed plasma ion and molecular ion flow profile (m s^-1)
    #		  fH2BC	- fltarr(nvr,nvx), this is an input boundary condition
    #			  specifying the shape of the neutral molecule velocity distribution 
    #			  function at location x(0). Normalization is arbitrary.
    #		          Only values with positive vx, fH2BC(*,nvx/2:*) are used
    #		          by the code.
    #	     GammaxH2BC	- float, desired neutral molecule flux density in the +Vx
    #			  direction at location x(0) (m^-2 s^-1)
    #			  fH2BC is scaled to yield this flux density.
    #	     NuLoss	- fltarr(nx), characteristic loss frequency for HP molecules (1/s)
    #			  (for open field lines, this is ~Cs/L). If this variable is undefined,
    #			  then NuLoss set set to zero.
    #	     PipeDia - fltarr(nx), effective pipe diameter (meters)
    #			  This variable allows collisions with the 'side-walls' to be simulated.
    #			  If this variable is undefined, then PipeDia set set to zero. Zero values
    #			  of PipeDia are ignored (i.e., treated as an infinite diameter).
    #        fH	- fltarr(nvr,nvx,nx), neutral atomic velocity distribution
    #             function. fH is normalized so that the atomic neutral density, nH(k), is 
    #			  defined as the velocity space integration:
    #                           nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
    #			  If this variable is undefined, then no molecule-atom collisions are included.
    #			  Note: dVx is velocity space differential for Vx axis and Vr2pidVr = Vr*!pi*dVr
    #		            with dVr being velocity space differential for Vr axis.
    #		    SH2	- fltarr(nx), source profile of wall-temperature (room temperature) H2 molecules (m^-3 s^-1)
    #			  If this variable is undefined, it is set equal to zero.
    #
    #  Input & Output:
    #                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
    #                         function.
    #			  'Seed' values for this may be specified on input. If this parameter
    #			  is undefined, then a zero 'seed' value will be used. 
    #			  The algorithm outputs a self-consistent fH2.
    #			  fH2 is normalized so that the neutral density, nH2(k), is defined as 
    #			  the velocity space integration:
    #                             nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
    #		    nHP	- fltarr(nx), molecular ion density profile (m^-3)
    #			     'Seed' values for this may be specified on input. If this parameter
    #			     is undefined, then a zero 'seed' value will be used. 
    #			     The algorithm outputs a self-consistent profile for nHP.
    #		    THP	- fltarr(nx), molecular ion temperature profile (m^-3)
    #			     'Seed' values for this may be specified on input. If this parameter
    #			     is undefined, then a 'seed' value of 3 eV will be used. 
    #			     The algorithm outputs a self-consistent profile for THP.
    #
    #  Output:
    #                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
    #                         function. fH2 is normalized so that the neutral density, nH2(k), is 
    #			  defined as the velocity space integration:
    #                             nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
    #                   nH2	- fltarr(nx), neutral molecule density profile (m^-3)
    #              GammaxH2	- fltarr(nx), neutral flux profile (# m^-2 s^-1)
    #                             computed from GammaxH2(k)=Vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))
    #                  VxH2	- fltarr(nx), neutral velocity profile (m s^-1)
    #                             computed from GammaxH2/nH2
    #
    #                       To aid in computing the some of the quantities below, the procedure internally
    #                       defines the quantities:
    #                       vr2vx2_ran(i,j,k)=vr(i)^2+(vx(j)-VxH2(k))^2
    #                                     which is the magnitude of 'random v^2' at each mesh point
    #                       vr2vx2(i,j,k)=vr(i)^2+vx(j)^2
    #                                     which is the magnitude of 'total v^2' at each mesh point
    #                       q=1.602177D-19, mH=1.6726231D-27
    #                       C(*,*,*) is the right hand side of the Boltzmann equation, evaluated
    #                                using the computed neutral distribution function
    #
    #                   pH2	- fltarr(nx), neutral pressure (eV m^-2) computed from:
    #                         pH2(k)~vth2*total(Vr2pidVr*(vr2vx2_ran(*,*,k)*fH2(*,*,k))#dVx))*(mu*mH)/(3*q)
    #                   TH2	- fltarr(nx), neutral temperature profile (eV) computed from: TH2=pH2/nH2
    #                  qxH2	- fltarr(nx), neutral random heat flux profile (watts m^-2) computed from:
    #                          qxH2(k)~vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH2(*,*,k))#(dVx*(vx-VxH2(k)))))*0.5*(mu*mH)
    #            qxH2_total	- fltarr(nx), total neutral heat flux profile (watts m^-2)
    #                             This is the total heat flux transported by the neutrals:
    #                         qxH2_total=(0.5*nH2*(mu*mH)*VxH2*VxH2 + 2.5*pH2*q)*VxH2 + piH2_xx*VxH2 + qxH2
    #                 Sloss	- fltarr(nx), H2 loss rate from ionization and dissociation (SH2) (m^-3 s^-1)
    #                   QH2	- fltarr(nx), rate of net thermal energy transfer into neutral molecules (watts m^-3) computed from
    #                               QH2(k)~vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C(*,*,k))#dVx))*0.5*(mu*mH)
    #                  RxH2	- fltarr(nx), rate of x momentum transfer to neutral molecules (=force, N m^-2).
    #                               RxH2(k)~Vth*total(Vr2pidVr*(C(*,*,k)#(dVx*(vx-VxH2(k)))))*(mu*mH)
    #             QH2_total	- fltarr(nx), net rate of total energy transfer into neutral molecules
    #                          = QH2 + RxH2*VxH2 - 0.5*(mu*mH)*(Sloss-SH2)*VxH2*VxH2 (watts m^-3)
    #              AlbedoH2	- float, Ratio of molceular particle flux with Vx < 0 divided by particle flux
    #                          with Vx > 0  at x=x(0)
    #                          (Note: For SH2 non-zero, the flux with Vx < 0 will include
    #                          contributions from molecular hydrogen sources within the 'slab'.
    #                          In this case, this parameter does not return the true 'Albedo'.)
    #	         WallH2 - fltarr(nx), molecular sink rate and source rate from interation with 'side walls' (m^-3 s^-1)
    #			  Note: There is no loss of molecules when they strike the side walls in this model. The 
    #		          molecules are just relaunched with a Maxwellian distribution at the wall temperature.
    #
    #	            nHP	- fltarr(nx), molecular ion density profile (m^-3)
    #			     'Seed' values for this may be specified on input.
    #		    THP	- fltarr(nx), molecular ion temperature profile (m^-3)
    #			     'Seed' values for this may be specified on input.
    #		    fSH	- fltarr(nvr,nvx,nx), H source velocity distribution.
    #			  fSH is normalized so that the total atomic neutral 
    #			  source, SH(k), is defined as the velocity space integration:
    #			      SH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
    #		     SH	- fltarr(nx), atomic neutral source profile (m^-3 s^-1)
    #  		             computed from total(Vr2pidVr*(fSH(*,*,k)#dVx))
    #		     SP	- fltarr(nx), proton source profile (m^-3 s^-1)
    #		    SHP	- fltarr(nx), molecular ion source profile (m^-3 s^-1)
    #		    NuE	- fltarr(nx), energy equilibration frequency of molecular ion with bulk plasma (1/s)
    #		  NuDis	- fltarr(nx), molecular ion dissociation frequency (1/s)
    #
    # KEYWORDS:
    #   Output:
    #		    ESH	- fltarr(nvr,nx) - returns normalized H source energy distribution
    #			     derived from velocity-space distribution:
    #  		            ESH(*,*) = 0.5*mu*mH*Vr^2 FSH(*,vx=0,*)*4*pi*Vr/(mu*mH)
    #		  Eaxis	- fltarr(nvr) - returns energy coordinate for ESH (eV) = 0.5*mu*mH*Vr^2
    #	          error	- Returns error status: 0=no error, solution returned
    #					        1=error, no solution returned
    #
    # COMMON BLOCK KINETIC_H2_OUTPUT
    #    Output:
    #               piH2_xx	- fltarr(nx), xx element of stress tensor (eV m^-2) computed from:
    #                         piH2_xx(k)~vth2*total(Vr2pidVr*(fH2(*,*,k)#(dVx*(vx-VxH2(k))^2)))*(mu*mH)/q - pH2
    #               piH2_yy	- fltarr(nx), yy element of stress tensor (eV m^-2) computed from:
    #                         piH2_yy(k)~vth2*total((Vr2pidVr*Vr^2)*(fH2(*,*,k)#dVx))*(mu*mH)/q - pH2
    #               piH2_zz	- fltarr(nx), zz element of stress tensor (eV m^-2) = piH2_yy
    #			   Note: cylindrical system relates r^2 = y^2 + z^2. All other stress tensor elements are zero.
    #
    #           The following momentum and energy transfer rates are computed from charge-exchange collsions between species:
    #                RxH2CX	- fltarr(nx), rate of x momentum transfer from molecular ions to neutral molecules (=force/vol, N m^-3).
    #                 EH2CX	- fltarr(nx), rate of energy transfer from molecular ions to neutral molecules (watts m^-3).
    #
    #           The following momentum and energy transfer rates are computed from elastic collsions between species:
    #                RxH_H2	- fltarr(nx), rate of x momentum transfer from neutral atoms to molecules (=force/vol, N m^-3).
    #                RxP_H2	- fltarr(nx), rate of x momentum transfer from hydrogen ions to neutral molecules (=force/vol, N m^-3).
    #                 EH_H2	- fltarr(nx), rate of energy transfer from neutral atoms to molecules (watts m^-3).
    #                 EP_H2	- fltarr(nx), rate of energy transfer from hydrogen ions to neutral molecules (watts m^-3).
    #
    #           The following momentum and energy transfer rates are computed from collisions with the 'side-walls'
    #                RxW_H2	- fltarr(nx), rate of x momentum transfer from wall to molecules (=force/vol, N m^-3).
    #                 EW_H2	- fltarr(nx), rate of energy transfer from wall to neutral molecules (watts m^-3).
    #
    #           The following is the rate of parallel to perpendicular energy transfer computed from elastic collisions
    #       Epara_PerpH2_H2	- fltarr(nx), rate of parallel to perp energy transfer within molecular hydrogen species (watts m^-3).
    #
    # KEYWORDS:
    #   Input:
    #	       truncate	- float, stop computation when the maximum 
    #			  increment of neutral density normalized to 
    #			  inputed neutral density is less than this 
    #	    		  value in a subsequent generation. Default value is 1.0e-4
    #
    #             Simple_CX	- if set, then use CX source option (B): Neutral molecules are born
    #                         in velocity with a distribution proportional to the local
    #                         molecular ion distribution function. Simple_CX=1 is default.
    #
    #                         if not set, then use CX source option (A): The CX source
    #                         neutral molecule distribution function is computed by evaluating the
    #                         the CX cross section for each combination of (vr,vx,vr',vx')
    #                         and convolving it with the molecule neutral distribution function.
    #                         This option requires more CPU time and memory.
    #
    #      	        Max_gen	- integer, maximum number of collision generations to try including before giving up.
    #                         Default is 50.
    #      Compute_H_Source	-  if set, then compute fSH, SH, SP, and SHP
    #
    #	      No_Sawada	- if set, then DO NOT correct reaction rates according to
    #			     results from collisional-radiative model of Sawada
    #			     [Sawada, K. and Fujimoto, T., Journal of Applied Physics 78 (1995) 2913.]
    #	       H2_H2_EL	- if set, then include H2 -> H2 elastic self collisions
    #			     Note: if H2_H2_EL is set, then algorithm iterates fH2 until
    #	                     self consistent fH2 is achieved.
    #	       H2_HP_CX	- if set, then include H2 -> H2(+) charge exchange collisions
    #			     Note: if H2_HP_CX is set, then algorithm iterates until
    #	                     self consistent nHp is achieved.
    #	        H2_P_EL	- if set, then include H2 -> H(+) elastic collisions (does not lead to an iteration loop)
    #	        H2_H_EL	- if set, then include H2 -> H elastic collisions (does not lead to an iteration loop)
    #            ni_correct	- if set, then algorithm corrects hydrogen ion density
    #			     according to quasineutrality: ni=ne-nHp
    #			     This is done in an iteration loop.
    #
    #	 Compute_Errors	- if set, then return error estimates in common block KINETIC_H2_ERRORS below
    #
    #		   plot	- 0= no plots, 1=summary plots, 2=detail plots, 3=very detailed plots
    #		  debug	- 0= do not execute debug code, 1=summary debug, 2=detail debug, 3=very detailed debug
    #	        debrief	- 0= do not print, 1=print summary information, 2=print detailed information
    #	          pause	- if set, then pause between plots
    #
    # COMMON BLOCK KINETIC_H2_ERRORS
    #
    #	if COMPUTE_ERRORS keyword is set then the following is returned in common block KINETIC_H2_ERRORS
    #
    #	         Max_dx	- float(nx), Max_dx(k) for k=0:nx-2 returns maximum 
    #			  allowed x(k+1)-x(k) that avoids unphysical negative 
    #			  contributions to fH2
    #	     Vbar_error	- float(nx), returns numerical error in computing
    #			  the speed of neutrals averged over maxwellian distribution#
    #			  over a temperature range spanning the Franck Condon energies
    #		          of reactions R2, R3, R4, R5, R6, R7, R8, R10
    #			  The average speed should be:
    #				 vbar_exact=2*Vth*sqrt(Ti(*)/Tnorm)/sqrt(!pi)
    #			  Vbar_error returns: abs(vbar-vbar_exact)/vbar_exact
    #			  where vbar is the numerically computed value.
    #	     mesh_error	- fltarr(nvr,nvx,nx), normalized error of solution
    #			  based on substitution into Boltzmann equation.
    #	   moment_error	- fltarr(nx,m), normalized error of solution
    #			  based on substitution into velocity space
    #			  moments (v^m) of Boltzmann equation, m=[0,1,2,3,4]
    #      	        C_error	- fltarr(nx), normalized error in charge exchange and elastic scattering collision 
    #				      operator. This is a measure of how well the charge exchange and
    #				      elastic scattering portions of the collision operator
    #				      conserve particles.
    #      	       CX_error	- fltarr(nx), normalized particle conservation error in charge exchange collision operator.
    #	    H2_H2_error	- fltarr(nx,[0,1,2]) return normalized errors associated with 
    #		          particle [0], x-momentum [1], and total energy [2] convervation of the elastic self-collision operator
    #
    #      	   Source_error	- fltarr(nx), normalized error in mass continuity equation. This is a measure
    #				      of how well atomic plus molecular mass continuity is satisfied.
    #      qxH2_total_error	- fltarr(nx), normalized error estimate in computation of qxH2_total
    #       QH2_total_error	- fltarr(nx), normalized error estimate in computation of QH2_total
    #
    # History:
    #	22-Dec-2000 - B. LaBombard - first coding.
    #	06-Feb-2001 - B. LaBombard - added elastic collisions and molecular sources


    prompt = 'Kinetic_H2 => '
    # Kinetic_H2_input common - will change later 
    vx_s = KH2.Input.vx_s
    vr_s = KH2.Input.vr_s
    x_s = KH2.Input.x_s
    Tnorm_s = KH2.Input.Tnorm_s
    mu_s = KH2.Input.mu_s
    Ti_s = KH2.Input.Ti_s
    Te_s = KH2.Input.Te_s
    n_s = KH2.Input.n_s
    vxi_s = KH2.Input.vxi_s
    fH2BC_s = KH2.Input.fH2BC_s
    GammaxH2BC_s = KH2.Input.GammaxH2BC_s
    Nuloss  =  KH2.Input.NuLoss_s
    PipeDia_s = KH2.Input.PipeDia_s
    fH_s = KH2.Input.fH_s
    SH2_s = KH2.Input.SH2_s
    fH2_s = KH2.Input.fH2_s

    nHP_s = KH2.Input.nHP_s
    THP_s = KH2.Input.THP_s
    Simple_CX_s = KH2.Input.Simple_CX_s
    Sawada_s = KH2.Input.Sawada_s
    H2_H2_EL_s = KH2.Input.H2_H2_EL_s
    H2_P_EL_s = KH2.Input.H2_P_EL_s
    H2_H_EL_s = KH2.Input.H2_H_EL_s
    H2_HP_CX_s = KH2.Input.H2_HP_CX_s
    ni_correct_s = KH2.Input.ni_correct_s

    # kinetic_h2_internal common block - will change later 
    vr2vx2 = KH2.Internal.vr2vx2
    vr2vx_vxi2 = KH2.Internal.vr2vx_vxi2
    fw_hat = KH2.Internal.fw_hat
    fi_hat = KH2.Internal.fi_hat
    fHp_hat = KH2.Internal.fHp_hat
    EH2_P = KH2.Internal.EH2_P
    sigv = KH2.Internal.sigv
    Alpha_Loss = KH2.Internal.Alpha_Loss
    v_v2 = KH2.Internal.v_v2
    v_v = KH2.Internal.v_v
    vr2_vx2 = KH2.Internal.vr2_vx2
    vx_vx = KH2.Internal.vx_vx

    Vr2pidVrdVx = KH2.Internal.Vr2pidVrdVx
    SIG_CX = KH2.Internal.SIG_CX
    SIG_H2_H2 = KH2.Internal.SIG_H2_H2
    SIG_H2_H = KH2.Internal.SIG_H2_H
    SIG_H2_P = KH2.Internal.SIG_H2_P
    Alpha_CX = KH2.Internal.Alpha_CX
    Alpha_H2_H = KH2.Internal.Alpha_H2_H
    MH2_H2_sum = KH2.Internal.MH2_H2_sum
    Delta_nH2s = KH2.Internal.Delta_nH2s

    # kinetic_h2_moments common block - will change later 
    nH2 = KH2.Moments.nH
    VxH2 = KH2.Moments.VxH
    TH2 = KH2.Moments.TH

    # internal debug switches 
    shifted_Maxwellian_debug = 0
    CI_Test = 1
    Do_Alpha_CX_Test = 0

    # internal tolerances 
    DeltaVx_tol = .01
    Wpp_tol = .001

    # Test input parameters 
    if debug > 0: 
        plot = np.maximum(plot, 1)
        debrief = np.maximum(debrief, 1)
        pause = 1   # not sure if this is the correct translation 
    if No_Sawada == 0:
        Sawada = 1
    else:
        Sawada = 0
    error = 0 
    nvr = np.size(vr)
    nvx = np.size(vx)
    nx = np.size(x)
    dx = x - np.roll(x, 1) ; dx = dx[1:]
    notpos = np.argwhere(dx <= 0.0)
    count = np.size(notpos)
    # print("vr", vr)
    # print("vx", vx)
    # print("x", x)
    # print("dx", dx)
    # input()
    if count > 0:
        # print(prompt, 'x[:] must be increasing with index!')
        # error = 1
        raise Exception(prompt + " x[:] must be increasing with index!")
    if (nvx % 2) != 0:
        # print(prompt, 'Number of elements in vx must be even!')
        # error = 1
        raise Exception(prompt + " Number of elements in vx must be even!")
    if np.size(Ti) != nx:
        # print(prompt, 'Number of elements in Ti and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in Ti and x do not agree!")
    if vxi is None:
        vxi = np.zeros(nx)
    if np.size(vxi) != nx:
        # print(prompt, 'Number of elements in vxi and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in vxi and x do not agree!")
    if np.size(Te) != nx:
        # print(prompt, 'Number of elements in Te and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in Te and x do not agree!")
    if np.size(n) != nx:
        # print(prompt, 'Number of elements in n and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in n and x do not agree!")
    if NuLoss is None:
        NuLoss = np.zeros(nx)
    if np.size(NuLoss) != nx:
        # print(prompt, 'Number of elements in NuLoss and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in NuLoss and x do not agree!")
    if PipeDia is None:
        PipeDia = np.zeros(nx)
    if np.size(PipeDia) != nx:
        # print(prompt, 'Number of elements in PipeDia and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in PipeDia and x do not agree!")
    if GammaxH2BC is None:
        # print(prompt, 'GammaxH2BC is not defined!')
        # error = 1
        raise Exception(prompt + " GammaxH2BC is not defined!")
    if fH is None:
        fH = np.arange((nvr, nvx, nx))
    if len(fH[:,0,0]) != nvr:
        # print(prompt, 'Number of elements in fH[0][0][:] and vr do not agree!') # come back and double check the error messages 
        # error = 1
        raise Exception(prompt + " Number of elements in fH[:,0,0] and vr do not agree!")
    if len(fH[0,:,0]) != nvx:
        # print(prompt, 'Number of elements in fH[0][:] and vx do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH[0,:,0] and vx do not agree!")
    if len(fH[0,0,:]) != nx:
        # print(prompt, 'Number of elements in fH[0] and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH[0,0,:] and x do not agree!")
    if len(fH2BC[:,0]) != nvr:
        # print(prompt, 'Number of elements in fH2BC[0][:] and vr do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH2BC[:,0] and vr do not agree!")
    if len(fH2BC[0,:]) != nvx:
        # print(prompt, 'Number of elements in fH2BC[0] and vx do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH2BC[0,:] and vx do not agree!")
    if fH2 is None:
        fH2 = np.zeros((nvr,nvx,nx))
    if len(fH2[:,0,0]) != nvr:
        # print(prompt, 'Number of elements in fH2[0][0][:] and vr do not agree!') # come back and double check the error messages 
        # error = 1
        raise Exception(prompt + " Number of elements in fH2[:,0,0] and vr do not agree!")
    if len(fH2[0,:,0]) != nvx:
        # print(prompt, 'Number of elements in fH2[0][:][0] and vx do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH2[0,:,0] and vx do not agree!")
    if len(fH2[0,0,:]) != nx: 
        # print(prompt, 'Number of elements in fH2[:][0][0] and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in fH2[0,0,:] and x do not agree!")
    if SH2 is None:
        SH2 = np.zeros(nx)
    if np.size(SH2) != nx:
        # print(prompt, 'Number of elements in SH2 and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in SH2 and x do not agree!")
    if nHP is None:
        nHP = np.zeros(nx)
    if np.size(nHP) != nx:
        # print(prompt, 'Number of elements in nHP and x do not agree!')
        # error = 1
        raise Exception(prompt + " Number of elements in nHP and x do not agree!")
    if THP is None:
        THP = np.zeros(nx) + 3.0
    if np.size(THP) != nx:
        # print(prompt, 'Number of elements in THP and x do not agree!')
        # error = 1 
        raise Exception(prompt + " Number of elements in THP and x do not agree!")
    if np.sum(abs(vr)) <= 0.0:
        # print(prompt, 'vr is all 0!')
        # error = 1
        raise Exception(prompt + " vr is all 0!")
    
    count = np.size(np.argwhere(vr < 0))
    if count > 0:
        # print(prompt, 'vr contains zero or negative element(s)!')
        # error = 1
        raise Exception(prompt + " vr contains zero or negative element(s)!")
    if np.sum(abs(x)) <= 0.0:
        # print(prompt, 'vx is all 0!')
        # error = 1
        raise Exception(prompt + " vx is all 0!")
    if np.sum(x) <= 0.0:
        # print(prompt, 'Total(x) is less than or equal to 0!')
        # error = 1
        raise Exception(prompt + " Total(x) is less than or equal to 0!")
    if Tnorm is None:
        # print(prompt, 'Tnorm is not defined!')
        # error = 1
        raise Exception(prompt + " Tnorm is not defined!")
    if mu is None:
        # print(prompt, 'mu is not defined!')
        # error = 1
        raise Exception(prompt + " mu is not defined!")
    if mu != 1 and mu != 2:
        # print(prompt, 'mu must be 1 or 2!')
        # error = 1
        raise Exception(prompt + " mu must be 1 or 2!")

    #NOTE Removed Plotting formatting, bring back once the program actually works

    #NOTE Fix how these are set up, temporary fix included later, but this is causing problems
    # Creates an array with shape (N X 1) instead of just an array of length N
    i_n = np.argwhere(vx < 0 ) # fixed typo - GG
    count = np.size(i_n)
    if count < 1:
        # print(prompt, 'vx contains no negative elements!')
        # error = 1
        raise Exception(prompt + " vx contains no negative elements!")
    i_p = np.argwhere(vx > 0)
    count = np.size(i_p)
    if count < 1:
        # print(prompt, 'vx contains no positive elements!')
        # error = 1
        raise Exception(prompt + " vx contains no positive elements!")
    i_z = np.argwhere(vx == 0)
    count = np.size(i_z)
    if count > 0:
        # print(prompt, 'vx contains one or more zero elements!')
        # error = 1
        raise Exception(prompt + " vx contains one or more zero elements!")
    diff = np.argwhere(vx[i_p] != -np.flipud(vx[i_n]))
    count = np.size(diff)
    if count > 0:
        # print(prompt, 'vx array elements are not symmetric about zero!')
        # error = 1
        raise Exception(prompt + " vx array elements are not symmetric about zero!")
    fH2BC_input = np.zeros(fH2BC.shape)
    fH2BC_input[:,i_p] = fH2BC[:,i_p]
    test = np.sum(fH2BC_input)
    if test <= 0.0:
        # print(prompt, 'Values for fH2BC(*,*) with vx > 0 are all zero!')
        # error = 1
        raise Exception(prompt + " Values for fH2BC(:,:) with vx > 0 are all zero!")
    # print("i_p", i_p)
    # print("i_n", i_n)
    # print("i_z", i_z)
    # input()

    # Output variables 
    nH2 = np.zeros(nx)
    GammaxH2 = np.zeros(nx)
    VxH2 = np.zeros(nx)
    pH2 = np.zeros(nx)
    TH2 = np.zeros(nx)
    NuDis = np.zeros(nx)
    NuE = np.zeros(nx)

    qxH2 = np.zeros(nx)
    qxH2_total = np.zeros(nx)
    Sloss = np.zeros(nx)
    WallH2 = np.zeros(nx)
    QH2 = np.zeros(nx)
    RxH2 = np.zeros(nx)
    QH2_total = np.zeros(nx)
    piH2_xx = np.zeros(nx)
    piH2_yy = np.zeros(nx)
    piH2_zz = np.zeros(nx)
    RxH2CX = np.zeros(nx)
    RxH_H2 = np.zeros(nx)
    RxP_H2 = np.zeros(nx)
    RxW_H2 = np.zeros(nx)
    EH2CX = np.zeros(nx)
    EH_H2 = np.zeros(nx)
    EP_H2 = np.zeros(nx)
    EW_H2 = np.zeros(nx)
    Epara_PerpH2_H2 = np.zeros(nx)
    AlbedoH2 = 0.0

    fSH = np.zeros((nvr,nvx,nx))
    SH = np.zeros(nx)
    SP = np.zeros(nx)
    SHP = np.zeros(nx)
    ESH = np.zeros((nvr,nx))
    Eaxis = np.zeros(nx)

    # Internal Varibales 

    Work = np.zeros(nvr * nvx)
    fH2G = np.zeros((nvr,nvx,nx))
    NH2G = np.zeros((nx, Max_Gen + 1))
    Vth = np.sqrt(2 * CONST.Q * Tnorm / (mu * CONST.H_MASS))
    Vth2 = Vth * Vth
    Vth3 = Vth2 * Vth
    fH2s = np.zeros(nx)
    nH2s = np.zeros(nx)
    THPs = np.zeros(nx)
    nHPs = np.zeros(nx)
    Alpha_H2_H2 = np.zeros((nvr,nvx))
    Omega_H2_P = np.zeros(nx)
    Omega_H2_H = np.zeros(nx)
    Omega_H2_H2 = np.zeros(nx)
    VxH2G = np.zeros(nx)
    TH2G = np.zeros(nx)
    Wperp_paraH2 = np.zeros(nx)
    vr2vx2_ran2 = np.zeros((nvr,nvx))
    vr2_2vx_ran2 = np.zeros((nvr,nvx))
    vr2_2vx2_2D = np.zeros((nvr,nvx))
    RxCI_CX = np.zeros(nx)
    RxCI_H_H2 = np.zeros(nx)
    RxCI_P_H2 = np.zeros(nx)
    Epara_Perp_CI = np.zeros(nx)
    CI_CX_error = np.zeros(nx)
    CI_H_H2_error = np.zeros(nx)
    CI_P_H2_error = np.zeros(nx)
    CI_H2_H2_error = np.zeros(nx)
    Maxwell = np.zeros((nvr,nvx,nx))

    Vr2pidVr, VrVr4pidVr, dVx, VrL, VrR, VxL, VxR, vol, Vth_DeltaVx, Vx_DeltaVx, Vr_DeltaVr, Vr2Vx2, jpa, jpb, jna, jnb = make_dvr_dvx(vr, vx)

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
    
    # Vr^2-2*Vx^2
    for i in range(0, nvr):
        vr2_2vx2_2D[i,:] = (vr[i]**2) - 2*(vx**2)
    # print("vr2_2vx2_2D", vr2_2vx2_2D.T)
    # input()
    
    # Theta-prime Coordinate
    ntheta = 5      # use 5 theta mesh points for theta integration
    dTheta = np.ones(ntheta)/ntheta
    theta = np.pi*(np.arange(ntheta)/ntheta + 0.5/ntheta)
    # print("theta", theta)
    cos_theta = np.cos(theta)
    # print("cos", cos_theta)
    # input()

    # Determine Energy Space Differentials 
    Eaxis = Vth2*0.5*mu*CONST.H_MASS*(vr**2)/CONST.Q
    _Eaxis = np.append(Eaxis, 2*Eaxis[nvr - 1] - Eaxis[nvr - 2]) # changed to append to stop error - GG
    Eaxis_mid = np.append(0.0, 0.5*( _Eaxis + np.roll(_Eaxis, -1) )) # changed to append to stop error - GG
    dEaxis = np.roll(Eaxis_mid, -1) - Eaxis_mid
    dEaxis = dEaxis[0:nvr]
    # print("dEaxis", dEaxis)
    # input()

    # Scale input molecular distribution function to agree with desired flux
    gamma_input = 1.0
    if abs(GammaxH2BC) > 0:
        gamma_input = Vth*np.sum(Vr2pidVr*(fH2BC_input @ (vx*dVx)))
    ratio = abs(GammaxH2BC)/gamma_input
    fH2BC_input = fH2BC_input*ratio
    if abs(ratio - 1) > 0.01*truncate:
        fH2BC = fH2BC_input
    fH2[:,i_p,0] = fH2BC_input[:,i_p]
    # print("fH2", fH2.T)
    # print("shape", fH2.shape)
    # input()
    
    # if fh is zero, then turn off elastic H2 <-> H collisions
    H2_H_EL = H2_H_EL # fixed typo - GG #NOTE ???, seems to be input? IDL thing not needed here
    # print("H2_H_EL", H2_H_EL)
    if np.sum(fH) <= 0.0:
        H2_H_EL = 0
    # print("H2_H_EL", H2_H_EL)
    # input()

    # Set iteration Scheme 
    fH2_iterate = 0 
    if (H2_H2_EL != 0) or (H2_HP_CX != 0) or (H2_H_EL != 0) or (H2_P_EL != 0) or (ni_correct != 0):
        fH2_iterate=1
    fH2_generations = 0
    if fH2_iterate != 0:
        fH2_generations = 1
    
    # Set flags to make use of previously computed local parameters 
    New_Grid = 1
    if vx_s is not None:
        test = 0 
        ii = np.argwhere(vx_s != vx) ; test = test + np.size(ii)
        ii = np.argwhere(vr_s != vr) ; test = test + np.size(ii)
        ii = np.argwhere(x_s != x) ; test = test + np.size(ii)
        ii = np.argwhere(Tnorm_s != Tnorm) ; test = test + np.size(ii)
        ii = np.argwhere(mu_s != mu) ; test = test + np.size(ii)
        if test <= 0:
            New_Grid = 0
    New_Protons = 1
    if Ti_s is not None:
        test = 0 
        ii = np.argwhere(Ti_s != Ti)
        test = test + np.size(ii)
        ii = np.argwhere(n_s != n)
        test = test + np.size(ii)
        ii = np.argwhere(vxi_s != vxi)
        test = test + np.size(ii)
        if test <= 0:
            New_Protons = 0
    New_Electrons = 1
    if Te_s is not None:
        test = 0 
        ii = np.argwhere(Te_s != Te) ; test = test + np.size(ii)
        ii = np.argwhere(n_s != n) ; test = test + np.size(ii)
        if test <= 0:
            New_Electrons = 0
    New_fH = 1
    if fH_s is not None:
        ii = np.argwhere(fH_s != fH)
        if np.size(ii) <= 0:
            New_fH = 0
    New_Simple_CX = 1
    if Simple_CX_s is not None:
        ii = np.argwhere(Simple_CX_s != Simple_CX)
        if np.size(ii) <= 0:
            New_Simple_CX = 0
    New_H2_Seed=1
    if fH2_s is not None:
        ii = np.argwhere(fH2_s != fH2)
        if np.size(ii) <= 0:
            New_H2_Seed = 0
    New_HP_Seed=1
    if nHP_s is not None:
        test = 0
        ii = np.argwhere(nHP_s != nHP)
        test = test + np.size(ii)
        ii = np.argwhere(THP_s != THP)
        test = test + np.size(ii)
        if test <= 0:
            New_HP_Seed = 0
    New_ni_correct=1
    if ni_correct_s is not None:
        ii = np.argwhere(ni_correct_s != ni_correct)
        if np.size(ii) <= 0:
            New_ni_correct = 0 

    Do_sigv = New_Grid | New_Electrons
    Do_fH_moments = (New_Grid | New_fH) & (np.sum(fH) > 0.0)
    Do_Alpha_CX =   (New_Grid | (Alpha_CX is None) | New_HP_Seed | New_Simple_CX) & H2_HP_CX
    # Do_Alpha_CX is updated in fH2_iteration loop
    Do_SIG_CX =     (New_Grid | (SIG_CX is None) | New_Simple_CX) & (Simple_CX == 0) & Do_Alpha_CX
    Do_Alpha_H2_H = (New_Grid | (Alpha_H2_H is None) | New_fH) & H2_H_EL
    Do_SIG_H2_H =   (New_Grid | (SIG_H2_H is None)) & Do_Alpha_H2_H
    Do_SIG_H2_H2 =  (New_Grid | (SIG_H2_H2 is None)) & H2_H2_EL
    Do_Alpha_H2_P = (New_Grid | (not 'Alpha_H2_P' in locals()) | New_Protons | New_ni_correct) & H2_P_EL #COME BACK TO THIS Alpha_H2_p is not defined!!!!!
    # Do_Alpha_H2_P is updated in fH2_iteration loop
    Do_SIG_H2_P =   (New_Grid | (SIG_H2_P is None)) & Do_Alpha_H2_P
    Do_v_v2 =      (New_Grid or (v_v2 is None)) & (CI_Test | Do_SIG_CX | Do_SIG_H2_H | Do_SIG_H2_H2 | Do_SIG_H2_P)

    if debug > 0:
        print("Kinetic H2 Settings")
        print("H2_H2_EL", H2_H2_EL)
        print("H2_P_EL", H2_P_EL)
        print("H2_H_EL", H2_H_EL)
        print("H2_HP_CX", H2_HP_CX)
        print("New_Grid", New_Grid)
        print("New_Protons", New_Protons)
        print("New_Electrons", New_Electrons)
        print("New_fH", New_fH)
        print("New_Simple_CX", New_Simple_CX)
        print("New_H2_Seed", New_H2_Seed)
        print("New_HP_Seed", New_HP_Seed)
        print("New_ni_correct", New_ni_correct)
        print("Do_sigv", Do_sigv)
        print("Do_fH_moments", Do_fH_moments)
        print("Do_Alpha_CX", Do_Alpha_CX)
        print("Do_SIG_CX", Do_SIG_CX)
        print("Simple_CX", Simple_CX)
        print("Do_Alpha_H2_H", Do_Alpha_H2_H)
        print("Do_SIG_H2_H", Do_SIG_H2_H)
        print("Do_SIG_H2_H2", Do_SIG_H2_H2)
        print("Do_Alpha_H2_P", Do_Alpha_H2_P)
        print("Do_SIG_H2_P", Do_SIG_H2_P)
        print("Do_v_v2", Do_v_v2)
        input()

    nH = np.zeros(nx)
    VxH = np.zeros(nx)
    TH = np.zeros(nx) + 1.0

    if Do_fH_moments:
        if debrief > 1:
            print(prompt, 'Computing vx and T moments of fH')
        # Compute x flow velocity and temperature of atomic species
        for k in range(0, nx):
            nH[k] = np.sum(Vr2pidVr*(fH[:,:,k] @ dVx))
            if nH[k] > 0:
                VxH[k] = Vth*np.sum(Vr2pidVr*(fH[:,:,k] @ (vx*dVx)))/nH[k]
                for i in range(0, nvr):
                    vr2vx2_ran2[i,:] = vr[i]**2 + (vx - VxH[k]/Vth)**2
                TH[k] = (mu*CONST.H_MASS)*Vth2*np.sum((Vr2pidVr*((vr2vx2_ran2*fH[:,:,k]) @ dVx))) / (3*CONST.Q*nH[k])
    
    if New_Grid:
        if debrief > 1:
            print(prompt, 'Computing vr2vx2, vr2vx_vxi2, EH2_P')
    
        # Magnitude of total normalized v^2 at each mesh point
        vr2vx2 = np.zeros((nvr,nvx,nx))
        for i in range(0, nvr):
            for k in range(0, nx):
                vr2vx2[i,:,k] = vr[i]**2 + vx**2

        # Magnitude of total normalized (v-vxi)^2 at each mesh point
        vr2vx_vxi2 = np.zeros((nvr,nvx,nx))
        for i in range(0, nvr):
            for k in range(0, nx):
                vr2vx_vxi2[i,:,k] = vr[i]**2 + (vx - vxi[k]/Vth)**2

        # Molecular hydrogen ion energy in local rest frame of plasma at each mesh point
        EH2_P = CONST.H_MASS*vr2vx_vxi2*Vth2/ CONST.Q
        EH2_P = np.maximum(EH2_P, 0.1)      # sigmav_cx does not handle neutral energies below 0.1 eV
        EH2_P = np.minimum(EH2_P, 2.0e4)    # sigmav_cx does not handle neutral energies above 20 keV

        # Compute Maxwellian H2 distribution at the wall temperature
        fw_hat = np.zeros((nvr,nvx))
        
        # note: Molecular ions have 'normalizing temperature' of 2 Tnorm, i.e., in order to
        # achieve the same thermal velocity^2, a molecular ion distribution has to have twice the temperature 
        # as an atomic ion distribution

        if (np.sum(SH2) > 0) | (np.sum(PipeDia) > 0): # come back and double check if it should be a bitwise or logical opperator 
            if debrief > 1:
                print(prompt, 'Computing fw_hat')
            vx_shift = np.array([0.0])
            Tmaxwell = np.array([CONST.TWALL])
            mol = 2
            _maxwell = create_shifted_maxwellian(vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm)
            fw_hat = _maxwell[:,:,0]
        
        # print("EH2_P", EH2_P.T)
        # print("vr2vx2", vr2vx2.T)
        # print("vr2vx_vxi2", vr2vx_vxi2.T)
        # print("fw_hat", fw_hat.T)
        # input()

    if New_Protons:
        # Compute fi_hat 
        if debrief > 1:
            print(prompt, 'Computing fi_Hat')
        vx_shift = vxi
        Tmaxwell = Ti
        mol = 1
        Maxwell = create_shifted_maxwellian_include(vr, vx, Tnorm, vx_shift, Tmaxwell, shifted_Maxwellian_debug, mu, mol, nx, nvx, nvr, \
                                          Vth, Vth2, Maxwell, vr2_2vx_ran2, Vr2pidVr, dVx, vol, \
                                          Vth_DeltaVx, Vx_DeltaVx, Vr_DeltaVr, vr2_2vx2_2D, jpa, jpb, jna, jnb)
        fi_hat = copy.copy(Maxwell) #NOTE Modify create_shifted_maxwellian to create new array, so copying not necessary
        # print("fi_hat", fi_hat.T)
        # input()

    if Do_sigv:
        if debrief > 1:
            print(prompt, 'Computing sigv')

        # Compute sigmav rates for each reaction and optionally apply
        # CR model corrections of Sawada

        sigv = np.zeros((nx,11))

        # Reaction R1:  e + H2 -> e + H2(+) + e 
        sigv[:,1] = sigmav_ion_hh(Te)
        if Sawada:
            sigv[:,1] = sigv[:,1] * 3.7 / 2.0

        # Reaction R2:  e + H2 -> H(1s) + H(1s)
        sigv[:,2] = sigmav_h1s_h1s_hh(Te)
        if Sawada:
            # Construct Table 
            Te_table = np.log([5,20,100])
            Ne_table = np.log([1e14,1e17,1e18,1e19,1e20,1e21,1e22])
            fctr_table = np.zeros((7, 3))
            fctr_table[:,0] = np.array([2.2, 2.2, 2.1, 1.9, 1.2,  1.1,  1.05]) / 5.3
            fctr_table[:,1] = np.array([5.1, 5.1, 4.3, 3.1, 1.5,  1.25, 1.25]) / 10.05
            fctr_table[:,2] = np.array([1.3, 1.3, 1.1, 0.8, 0.38, 0.24, 0.22]) / 2.1
            _Te = Te
            _Te = np.maximum(_Te, 5)
            _Te = np.minimum(_Te, 100)
            _n = np.maximum(n, 1e14)
            _n = np.minimum(n, 1e22)
            fctr = path_interp_2d(fctr_table, Ne_table, Te_table, np.log(_n), np.log(_Te))
            sigv[:,2] = (1.0 + fctr)*sigv[:,2]
            # print("Te_table", Te_table)
            # print("Ne_table", Ne_table)
            # print("fctr_table", fctr_table.T)
            # print("_Te", _Te)
            # print("_n", _n)
            # print("sigv2", sigv[:,2])
            # input()
        
        # Reaction R3:  e + H2 -> e + H(1s) + H*(2s)
        sigv[:,3] = sigmav_h1s_h2s_hh(Te)

        # Reaction R4:  e + H2 -> e + p + H(1s)
        sigv[:,4] = sigmav_p_h1s_hh(Te)
        if Sawada:
            sigv[:,4] = sigv[:,4]*1.0/0.6

        # Reaction R5:  e + H2 -> e + H*(2p) + H*(2s)
        sigv[:,5] = sigmav_h2p_h2s_hh(Te)

        # Reaction R6:  e + H2 -> e + H(1s) + H*(n=3)
        sigv[:,6] = sigmav_h1s_hn3_hh(Te)

        # Reaction R7:  e + H2(+) -> e + p + H(1s)
        sigv[:,7] = sigmav_p_h1s_hp(Te)

        # Reaction R8:  e + H2(+) -> e + p + H*(n=2)
        sigv[:,8] = sigmav_p_hn2_hp(Te)

        # Reaction R9:  e + H2(+) -> e + p + p + e
        sigv[:,9] = sigmav_p_p_hp(Te)

        # Reaction R10:  e + H2(+) -> e + H(1s) + H*(n>=2)
        sigv[:,10] = sigmav_h1s_hn_hp(Te)
        
        # Total H2 destruction rate (normalized by vth) = sum of reactions 1-6
        Alpha_Loss = np.zeros(nx)
        Alpha_Loss[:] = n*np.sum(sigv[:,1:7], axis = 1) / Vth

        # print("sigv", sigv.T)
        # print("Alpha_Loss", Alpha_Loss)
        # input()

    # Set up arrays for charge exchange and elastic collision computations, if needed 
    if Do_v_v2 == 1:
        if debrief > 1:
            print(prompt, 'Computing v_v2, v_v, vr2_vx2, and vx_vx')
        # v_v2=(v-v_prime)^2 at each double velocity space mesh point, including theta angle
        v_v2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))

        # vr2_vx2=(vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2
        # at each double velocity space mesh point, including theta angle
        vr2_vx2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))
        for m in range(0, ntheta): # double check this nested for loop 
            for l in range(0, nvx):
                for k in range(0, nvr):
                    for i in range(0, nvr):
                        v_starter = vr[i]**2 + vr[k]**2 - 2*vr[i]*vr[k]*cos_theta[m]
                        v_v2[i,:,k,l,m] = v_starter + (vx[:] - vx[l])**2  # not super confident 
                        vr2_vx2[i,:,k,l,m] = v_starter - 2*(vx[:] - vx[l])**2
        # v_v=|v-v_prime| at each double velocity space mesh point, including theta angle
        v_v = np.sqrt(v_v2)
        # print("V_V", v_v.T[0,0,0])
        # print("v_v shape", v_v.shape)
        # input()

        # vx_vx=(vx-vx_prime) at each double velocity space mesh point
        vx_vx = np.zeros((nvr,nvx,nvr,nvx))
        for j in range(0,nvx):
            for l in range(0, nvx):
                vx_vx[:,j,:,l] = vx[j] - vx[l]  
        # print("vx_vx", vx_vx.T)
        # input()

        # Set Vr'2pidVr'*dVx' for each double velocity space mesh point
        Vr2pidVrdVx = np.zeros((nvr,nvx,nvr,nvx))
        for k in range(0, nvr):
            Vr2pidVrdVx[:,:,k,:] = Vr2pidVr[k]
        for l in range(0, nvx):
            Vr2pidVrdVx[:,:,:,l] = Vr2pidVrdVx[:,:,:,l]*dVx[l]
        # print("Vr2pidVrdVx", Vr2pidVrdVx.T)
        # input()
    
    if Simple_CX == 0 and Do_SIG_CX == 1: #NOTE Not Tested Yet
        if debrief > 1:
            print(prompt, 'Computing SIG_CX')
        #  Option (A) was selected: Compute SigmaV_CX from sigma directly.
        # In preparation, compute SIG_CX for present velocity space grid, if it has not 
        # already been computed with the present input parameters

        # compute SIGMA_CX * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx, ntheta))
        _Sig[:] = (v_v*sigma_cx_hh(v_v2*(CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')

        # Set SIG_CX = vr' x Integral{v_v*sigma_cx} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        SIG_CX = np.zeros((nvr*nvx, nvr*nvx))
        SIG_CX[:] = (Vr2pidVrdVx*((_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F'))).reshape(SIG_CX.shape, order='F') 

        # SIG_CX is now vr' * sigma_cx(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

    if Do_SIG_H2_H == 1: #NOTE Not Tested Yet
        if debrief > 1:
            print(prompt, 'Computing SIG_H2_P')
        # Compute SIG_H2_P for present velocity space grid, if it is needed and has not 
        # already been computed with the present input parameters

        # Compute sigma_H2_H * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx, ntheta))
        _Sig[:] = (v_v*sigma_el_h_hh(v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')

        # Note: using H energy here for cross-section tabulated as H -> H2
        # Set SIG_H2_H = vr' x vx_vx x Integral{v_v * sigma_H2_H} over theta = 0, 2pi times differential velocity space element Vr'2pidVr'*dVx
        SIG_H2_H = np.zeros((nvr*nvx, nvr*nvx))
        SIG_H2_H[:] = (Vr2pidVrdVx*vx_vx*((_Sig @ dTheta).reshape(vx_vx.shape, order='F'))).reshape(SIG_H2_H.shape, order='F')

        # SIG_H2_H is now vr' * vx_vx * sigma_H2_H(v_V) ( integrated over theta ) for all possible ([vr, vx], [vr', vx'])

    if Do_SIG_H2_P == 1:
        if debrief > 1:
            print(prompt, 'Computing SIG_H2_P') 
        #   Compute SIG_H2_P for present velocity space grid, if it is needed and has not 
        # already been computed with the present input parameters

        # Compute sigma_H2_P * v_v at all possible relative velocities
        _Sig = np.zeros((nvr*nvx*nvr*nvx, ntheta))
        _Sig[:] = (v_v*sigma_el_p_hh(v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')
        #energy = v_v2*(0.5*CONST.H_MASS*Vth2/CONST.Q)
        # print("energy", (v_v*sigma_el_p_hh(energy)).T)
        # print("shape", (v_v*sigma_el_p_hh(energy)).shape)
        # input()
        # print("_sig", _Sig.T)
        # print("shape", _Sig.shape)
        # input()

        # Note: using H energy here for cross-section tabulated as p -> H2

        # Set SIG_H2_P = vr' x vx_vx x Integral{v_v * sigma_H2_P} over theta = 0, 2pi times differential velocity space element Vr'2pidVr' * dVx
        SIG_H2_P = np.zeros((nvr*nvx, nvr*nvx))
        SIG_H2_P[:] = (Vr2pidVrdVx*vx_vx*(_Sig @ dTheta).reshape(vx_vx.shape, order='F')).reshape(SIG_H2_P.shape, order='F')

        # print("SIG_H2_P", SIG_H2_P.T)
        # print("shape", SIG_H2_P.shape)
        # input()

        # SIG_H2_P is now vr' * vx_vx * sigma_h2_P(v_v) * v_v (integrated over theta) for all possible ([vr, vx], [vr', vx'])

    if Do_SIG_H2_H2 == 1:
        if debrief > 1:
            print(prompt, 'Computing SIG_H2_H2')
        
        #   Compute SIG_H2_H2 for present velocity space grid, if it is needed and has not 
        # already been computed with the present input parameters

        # Compute sigma_H2_H2 * vr2_vx2 * v_v at all possible relative velocities 
        _Sig = np.zeros((nvr*nvx*nvr*nvx, ntheta))
        _Sig[:] = (vr2_vx2*v_v*sigma_el_hh_hh(v_v2*(CONST.H_MASS*mu*Vth2/CONST.Q), vis = 1)/8.0).reshape(_Sig.shape, order='F')
        # print("vr2_vx2", vr2_vx2)
        # print("v_v", v_v)
        # print("v_v2", v_v2)
        # print("_SIG", _Sig.T)
        # input()

        # Note : For viscosity, the cross section for D -> D is the same function of 
        # center of mass energy as H -> H.

        # Set SIG_H2_H2 = vr' x Integral{vr2_vx2*v_v*sigma_H2_H2} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
        SIG_H2_H2 = np.zeros((nvr * nvx, nvr * nvx))
        SIG_H2_H2[:] = (Vr2pidVrdVx*((_Sig @ dTheta).reshape(Vr2pidVrdVx.shape, order='F'))).reshape(SIG_H2_H2.shape, order='F')

        # print("SIG_H2_H2", SIG_H2_H2)
        # input()

        # SIG_H2_H2 is now vr' * sigma_H2_H2(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
    
    if Do_Alpha_H2_H == 1: #NOTE Not Tested Yet
        if debrief > 1:
            print(prompt, 'Computing Alpha_H2_H')
        
        # Compute Alpha_H2_H for inputed fH, if it is needed and has not
        # already been computed with the present input parameters

        Alpha_H2_H = np.zeros((nvr, nvx, nx))
        for k in range(0, nx):
            Work[:] = fH[:,:,k].reshape(Work.shape, order='F')
            Alpha_H2_H[:,:,k] = (SIG_H2_H @ Work).reshape(Alpha_H2_H[:,:,k].shape, order='F')

        # print("Alpha_H2_H", Alpha_H2_H)
        # input()

        
    # Compute nH2
    for k in range(0, nx):
        nH2[k] = np.sum(Vr2pidVr*(fH2[:,:,k] @ dVx))

    if New_H2_Seed:
        MH2_H2_sum = np.zeros((nvr,nvx,nx))
        Delta_nH2s = 1.0

    gamma_wall = np.zeros((nvr,nvx,nx))
    for k in range(0, nx):
        if PipeDia[k] > 0.0:
            for j in range(0, nvx):
                gamma_wall[:,j,k] = 2*vr/PipeDia[k]
    # print("nH2", nH2)
    # print("gamma_wall", gamma_wall.T)
    # input()

    #fH2 Iteration - I dont know where the fH2_iterate is coming from 
    # This is the iteration entry point for fH2, THP and nHP iteration.
    # Save 'seed' values for comparison later
    do_fH2_iterate = True
    while do_fH2_iterate:
        do_fH2_iterate = False
        fH2s = copy.copy(fH2)
        nH2s = copy.copy(nH2)
        THPs = copy.copy(THP)
        nHPs = copy.copy(nHP)

        # Compute Alpha_CX for present THP and nHP, if it is needed and has not
        # already been computed with the present parameters
        if Do_Alpha_CX == 1:
            if debrief > 1:
                print(prompt, 'Computing Alpha_CX')
            # Set Maxwellian Molecular Ion Distrobution Function (assumed to be drifting with ion velocity, vxi)
            vx_shift = vxi
            Tmaxwell = THP
            mol = 2
            Maxwell = create_shifted_maxwellian_include(vr, vx, Tnorm, vx_shift, Tmaxwell, shifted_Maxwellian_debug, mu, mol, nx, nvx, nvr, \
                                            Vth, Vth2, Maxwell, vr2_2vx_ran2, Vr2pidVr, dVx, vol, \
                                            Vth_DeltaVx, Vx_DeltaVx, Vr_DeltaVr, vr2_2vx2_2D, jpa, jpb, jna, jnb)
            fHp_hat = copy.copy(Maxwell)
            # print("fHp_hat", fHp_hat)
            # input()

            if Simple_CX:
                # Option (B) : Use Maxwellian weighted <sigma v>
                
                # THP/mu at each mesh point
                THP_mu = np.zeros((nvr, nvx, nx))
                for k in range(0, nx):
                    THP_mu[:,:,k] = THP[k]/mu

                # Molecular Charge Exchange sink rate 
                Alpha_CX = sigmav_cx_hh(THP_mu, EH2_P)/Vth
                for k in range(0, nx):
                    Alpha_CX[:,:,k] = Alpha_CX[:,:,k]*nHP[k]
            else:
                # Option (A): Compute SigmaV_CX from sigma directly via SIG_CX
                Alpha_CX = np.zeros((nvr, nvx, nx))
                for k in range(0, nx):
                    Work[:] = (fHp_hat[:,:,k]*nHP[k]).reshape(Work.shape, order='F')
                    Alpha_CX[:,:,k] = SIG_CX @ Work
                if Do_Alpha_CX_Test:
                    alpha_cx_test = sigmav_cx_hh(THP_mu, EH2_P)/Vth
                    for k in range(0, nx):
                        alpha_cx_test[:,:,k] = alpha_cx_test[:,:,k]*nHP[k]
                        print('Compare alpha_cx and alpha_cx_test')
                        input()
            # print("Alpha_CX", Alpha_CX.T)
            # input()

        # Compute Alpha_H2_P for present Ti and ni (optionally correcting for nHP), 
        # if it is needed and has not already been computed with the present parameter
        if Do_Alpha_H2_P == 1:
            if debrief > 1:
                print(prompt, 'Computing Alpha_H2_P')
            Alpha_H2_P = np.zeros((nvr, nvx, nx))
            ni = n
            if ni_correct:
                ni = np.maximum((n-nHP), 0)
            # print("ni", ni)
            for k in range(0, nx):
                Work[:] = (fi_hat[:,:,k]*ni[k]).reshape(Work.shape, order='F')
                Alpha_H2_P[:,:,k] = (SIG_H2_P @ Work).reshape(Alpha_H2_P[:,:,k].shape, order='F')
            # print("fi_hat", fi_hat)
            # print("Work", Work)
            # print("Alpha_H2_P", Alpha_H2_P.T)
            # input()
        
        # Compute Omega values if nH2 is non-zero 

        ii = np.argwhere(nH2 <= 0) 
        if ii.size <= 0: #NOTE NOT Tested yet, return on iteration
            # compute VxH2
            if H2_P_EL or H2_H_EL or H2_H2_EL:
                for k in range(0, nx):
                    VxH2[k] = Vth*np.sum(Vr2pidVr*(fH2[:,:,k] @ (vx*dVx))) / nH2[k]
            # print("VxH2", VxH2)
            # input()

            # compute Omega_H2_P for present fH2 and Alpha_H2_P if H2_P elastic collisions are included
            if H2_P_EL:
                if debrief > 1:
                    print(prompt, 'Computing Omega_H2_P')
                for k in range(0, nx):
                    DeltaVx = (VxH2[k] - vxi[k])/Vth
                    MagDeltaVx = np.maximum(np.abs(DeltaVx), DeltaVx_tol)
                    DeltaVx = sign(DeltaVx)*MagDeltaVx
                    Omega_H2_P[k] = np.sum(Vr2pidVr*(((Alpha_H2_P[:,:,k]*fH2[:,:,k]) @ dVx)))/(nH2[k]*DeltaVx)
                Omega_H2_P =  np.maximum(Omega_H2_P, 0)
                # print("Omega_H2_P", Omega_H2_P)
                # input()

            # Compute Omega_H2_H for present fH2 and Alpha_H2_H if H2_H elastic collisions are included
            if H2_H_EL: #NOTE Note Tested Yet
                if debrief>1:
                    print(prompt+'Computing Omega_H2_H')
                for k in range(nx):
                    DeltaVx = (VxH2[k] - VxH[k])/Vth
                    MagDeltaVx = np.maximum(np.abs(DeltaVx), DeltaVx_tol)
                    DeltaVx = sign(DeltaVx)*MagDeltaVx
                    Omega_H2_H[k] = np.sum(Vr2pidVr*((Alpha_H2_H[:,:,k]*fH2[:,:,k]) @ dVx)/(nH2[k]*DeltaVx))
                Omega_H2_H = np.maximum(Omega_H2_H, 0)
                # print("Omega_H2_H", Omega_H2_H)
                # input()

            # Compute Omega_H2_H2 for present fH2 if H2_H2 elastic collisions are included
            if H2_H2_EL:
                if debrief > 1:
                    print(prompt, 'Computing Omega_H2_H2')
                if np.sum(MH2_H2_sum) < 0:
                    for k in range(0, nx):
                        for i in range(0, nvr):
                            vr2_2vx_ran2[i,:] = vr[i]**2 - 2*(vx - VxH2[k]/Vth)**2
                        Wperp_paraH2[k] = np.sum(Vr2pidVr*((vr2_2vx_ran2*fH2[:,:,k]) @ dVx))/nH2[k]
                else:
                    for k in range(0, nx):
                        M_fH2 = MH2_H2_sum[:,:,k] - fH2[:,:,k]
                        Wperp_paraH2[k] = -np.sum(Vr2pidVr*((vr2_2vx2_2D*M_fH2) @ dVx))/nH2[k]

                for k in range(0, nx):
                    Work[:] = fH2[:,:,k].reshape(Work.shape, order='F')
                    Alpha_H2_H2[:] = (SIG_H2_H2 @ Work).reshape(Alpha_H2_H2.shape, order='F')
                    Wpp = Wperp_paraH2[k]
                    MagWpp = np.maximum(abs(Wpp), Wpp_tol)
                    Wpp = sign(Wpp)*MagWpp  
                    Omega_H2_H2[k] = np.sum(Vr2pidVr*((Alpha_H2_H2*Work.reshape(Alpha_H2_H2.shape, order='F')) @ dVx))/(nH2[k]*Wpp)
                    
                Omega_H2_H2 = np.maximum(Omega_H2_H2, 0)
                # print("Omega_H2_H2", Omega_H2_H2)
                # input()

        # Total Elastic scattering frequency
        Omega_EL = Omega_H2_P + Omega_H2_H + Omega_H2_H2
        # print("Omega_EL", Omega_EL)
        # input()

        # Total collision frequency
        alpha_c = np.zeros((nvr,nvx,nx))
        if H2_HP_CX:
            for k in range(0, nx):
                alpha_c[:,:,k] = Alpha_CX[:,:,k]+Alpha_Loss[k]+Omega_EL[k]+gamma_wall[:,:,k]
        else: 
            for k in range(0, nx):
                alpha_c[:,:,k] = Alpha_Loss[k]+Omega_EL[k]+gamma_wall[:,:,k]
        # print("alpha_c", alpha_c.T)
        # input()

        # Test x grid spacing based on Eq.(27) in notes
        if debrief > 1: 
            print(prompt, 'Testing x grid spacing')
        Max_dx = np.full(nx, 1.0E32)
        for k in range(0, nx) : 
            for j in range(i_p[0][0], nvx):
                denom = alpha_c[:,j,k]
                Max_dx[k] = np.minimum(Max_dx[k], np.min(2*vx[j]/denom)) 
        # print("Max_dx", Max_dx)
        # input()

        dx = np.roll(x,-1) - x
        Max_dxL = Max_dx[0:nx-1]
        Max_dxR = Max_dx[1:nx]
        Max_dx = np.minimum(Max_dxL, Max_dxR)
        ilarge = np.argwhere(Max_dx < dx[0:nx-1])
        # print("dx", dx)
        # print("Max_dxL", Max_dxL)
        # print("Max_dxR", Max_dxR)
        # print("Max_dx", Max_dx)
        # print("ilarge", ilarge.size)
        # input()

        if ilarge.size > 0:
            print(prompt,'x mesh spacing is too large!') #NOTE Check Formatting
            debug = 1
            out = ''
            jj = 0
            print(' x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)')
            for ii in range(0, np.size(ilarge)-1):
                jj = jj + 1
                out = out + str(ilarge[ii]) +' '+ str(x(ilarge[ii]+1)-x(ilarge[ii])) +' '+ Max_dx(ilarge[ii]) +' ' # I didn't include any of the formatting from the original code I thought this is something we could determine later
                if jj > 4:
                    print(out)
                    jj = 0 
                    out = ''
            if jj > 0: 
                print(out) 
            error = 1
            raise Exception("x mesh spacing is too large") 
        
        # Define parameters Ak, Bk, Ck, Dk, Fk, Gk
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
                Fk[:,j,k] = (x[k+1] - x[k])*fw_hat[:,j]*(SH2[k+1]+SH2[k]) / (Vth*denom)
        for k in range(1, nx):
            for j in range(0, i_p[0][0]):
                denom = -2*vx[j] + (x[k] - x[k-1])*alpha_c[:,j,k-1]
                Ck[:,j,k] = (-2*vx[j] - (x[k] - x[k -1])*alpha_c[:,j,k]) / denom
                Dk[:,j,k] = (x[k] - x[k-1]) / denom
                Gk[:,j,k] = (x[k] - x[k-1])*fw_hat[:,j]*(SH2[k]+SH2[k-1]) / (Vth*denom)
        # print("Ak", Ak.T)
        # print("Bk", Bk.T)
        # print("Ck", Ck.T)
        # print("Dk", Dk.T)
        # print("Fk", Fk.T)
        # print("Gk", Gk.T)
        # input()

        # Compute first-flight (0th generation) neutral distribution function
        Swall_sum = np.zeros((nvr,nvx,nx))
        Beta_CX_sum = np.zeros((nvr,nvx,nx))
        MH2_P_sum = np.zeros((nvr,nvx,nx))
        MH2_H_sum = np.zeros((nvr,nvx,nx))
        MH2_H2_sum = np.zeros((nvr,nvx,nx))
        igen = 0
        if debrief > 0:
            print(prompt, 'Computing molecular neutral generation#', sval(igen))
        fH2G[:,i_p,0] = fH2[:,i_p,0]
        for k in range(0, nx-1):
            fH2G[:,i_p,k+1] = fH2G[:,i_p,k]*Ak[:,i_p,k] + Fk[:,i_p,k]
        for k in range(nx-1,0,-1):
            fH2G[:,i_n,k-1] = fH2G[:,i_n,k]*Ck[:,i_n,k] + Gk[:,i_n,k]
        # print("fH2G", fH2G.T)
        # input()

        # Compute first-flight neutral density profile 
        for k in range(0, nx):
            NH2G[k, igen] = np.sum(Vr2pidVr*(fH2G[:,:,k] @ dVx))
        # print("NH2G", NH2G.T)
        # input()

        # if plot > 1: #NOTE Implement Plotting Later
        #     fH21d = np.zeros((nvx, nx)).T
        #     for k in range(0, nx - 1):
        #         fH21d[k] = np.dot(Vr2pidVr, fH2G[k])
        #     for i in range(0, nx-1):
        #         plt.plot(vx, fH21d[i]) #I cant really tell if this is plotting the data rge same way that the idl version does
        #     if debug > 0:
        #         return 


        # Set total molecular neutral distrobution function to first flight generation 
        fH2 = copy.copy(fH2G)
        nH2 = NH2G[:,0]
        # print("fH2", fH2.T)
        # print("nH2", nH2.T)
        # input()
        
        #fH2_done = 0 # I am not sure if this is correct because fH2_done is a function, but I'm not sure what the intention of the original IDL code was so I don't know how to change it - GG
        
# next_generation #########################################################################################################################################################################
        while True:
            if igen+1 > Max_Gen or fH2_generations == 0: 
                if debrief > 1:
                    print(prompt,'Completed ', sval(Max_Gen), ' generations. Returning present solution...')
                break
            igen = igen + 1
            if debrief > 0: 
                print(prompt, 'Computing molecular neutral generation#', sval(igen))
        
            #Compute Swall from previous generation
            Swall = np.zeros((nvr, nvx, nx))
            if np.sum(gamma_wall) > 0:
                if debrief > 1:
                    print(prompt, 'Computing Swall')
                for k in range(0, nx): 
                    Swall[:,:k] = fw_hat*np.sum(Vr2pidVr*((gamma_wall[:,:,k]*fH2G[:,:,k]) @ dVx))
                #Sum wall collision source over all generations
                Swall_sum = Swall_sum + Swall 
            # print("Swall_sum", Swall_sum.T)
            # input()

            #Compute Beta_CX from previous generation
            Beta_CX = np.zeros((nvr,nvx,nx))
            if H2_HP_CX: 
                if debrief > 1:
                    print(prompt, 'Computing Beta_CX')
                if Simple_CX:
                    # Option (B): Compute charge exchange source with assumption that CX source neutrals have 
                    # molecular ion distribution function
                    for k in range(0, nx): 
                        Beta_CX[:,:,k] = fHp_hat[:,:,k]*np.sum(Vr2pidVr*((Alpha_CX[:,:,k]*fH2G[:,:,k]) @ dVx))
                else: 
                    # Option (A): Compute charge exchange source using fH2 and vr x sigma x v_v at each velocity mesh point
                    for k in range(0, nx):
                        Work[:] = fH2G[:,:,k]
                        Beta_CX[:,:,k] = nHP[k]*fHp_hat[:,:,k]*(SIG_CX @ Work)
                #Sum 
                Beta_CX_sum = Beta_CX_sum + Beta_CX
            # print("Beta_CX_sum", Beta_CX_sum.T)
            # input()

            # Compute MH2 from previous generation
            MH2_H2 = np.zeros((nvr,nvx,nx))
            MH2_P = np.zeros((nvr,nvx,nx))
            MH2_H = np.zeros((nvr,nvx,nx))
            OmegaM = np.zeros((nvr,nvx,nx))
            if H2_H2_EL or H2_P_EL or H2_H_EL:
                # Compute VxH2G, TH2G
                for k in range(0, nx):
                    VxH2G[k] = Vth*np.sum(Vr2pidVr*(fH2G[:,:,k] @ (vx * dVx))) / NH2G[k,igen-1]
                    for i in range(0, nvr):
                        vr2vx2_ran2[i,:] = vr[i]**2 + (vx - VxH2G[k]/Vth)**2
                    TH2G[k] = (2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran2*fH2G[:,:,k]) @ dVx)) / (3*CONST.Q*NH2G[k,igen-1])
                # print("TH2G", TH2G)
                # input()

                if H2_H2_EL:
                    if debrief > 1: 
                        print(prompt, 'Computing MH2_H2')
                    # Compute MH2_H2
                    vx_shift = VxH2G
                    Tmaxwell = TH2G
                    mol = 2
                    Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                            nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                            Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
                    # print("Maxwell", Maxwell.T)
                    # input()

                    for k in range(0, nx):
                        MH2_H2[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_H2[k]*MH2_H2[:,:,k]
                    MH2_H2_sum = MH2_H2_sum + MH2_H2
                    # print("MH2_H2_sum", MH2_H2_sum.T)
                    # print("OmegaM", OmegaM.T)
                    # input()

                if H2_P_EL:
                    if debrief > 1:
                        print(prompt, 'Computing MH2_P')
                    # Compute MH2_P
                    vx_shift = (2 * VxH2G + vxi)/3
                    Tmaxwell = TH2G + (4/9)*(Ti - TH2G + ((mu*CONST.H_MASS*(vxi - VxH2G)**2)/(6*CONST.Q)))
                    mol = 2
                    Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                            nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                            Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
                    # print("Maxwell", Maxwell.T)
                    # input()

                    for k in range(0, nx):
                        MH2_P[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_P[k]*MH2_P[:,:,k]
                    MH2_P_sum = MH2_P_sum + MH2_P
                    # print("MH2_P_sum", MH2_P_sum.T)
                    # input()

                if H2_H_EL: #NOTE Not Tested Yet
                    if debrief > 1:
                        print(prompt, 'Computing MH2_H')
                    #Compute MH2_H
                    vx_shift = (2*VxH2G + VxH)/3
                    Tmaxwell = TH2G + (4/9)*(TH - TH2G + ((mu*CONST.H_MASS*(VxH - VxH2G)**2)/(6*CONST.Q)))
                    mol = 2
                    Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                            nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                            Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
                    # print("Maxwell", Maxwell.T)
                    # input()

                    for k in range(0, nx):
                        MH2_H[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen-1]
                        OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_H[k]*MH2_H[:,:,k]
                    MH2_H_sum = MH2_H_sum + MH2_H
                    # print("MH2_H_sum", MH2_H_sum.T)
                    # print("OmegaM", OmegaM.T)
                    # input()

            # Compute next generation molecular distribution
            fH2G[:] = 0.0
            for k in range(0, nx-1):
                fH2G[:,i_p,k+1] = Ak[:,i_p,k]*fH2G[:,i_p,k] + Bk[:,i_p,k]*(Swall[:,i_p,k+1] + Beta_CX[:,i_p,k+1] + OmegaM[:,i_p,k+1] + Swall[:,i_p,k] + Beta_CX[:,i_p,k] + OmegaM[:,i_p,k])
            for k in range(nx-1, 0, -1):
                fH2G[:,i_n,k-1] = Ck[:,i_n,k]*fH2G[:,i_n,k] + Dk[:,i_n,k]*(Swall[:,i_n,k-1] + Beta_CX[:,i_n,k-1] + OmegaM[:,i_n,k-1] + Swall[:,i_n,k] + Beta_CX[:,i_n,k] + OmegaM[:,i_n,k])
            for k in range(0, nx):
                NH2G[k,igen] = np.sum(Vr2pidVr*(fH2G[:,:,k] @ dVx))
            # print("fH2G", fH2G.T)
            # print("NH2G", NH2G.T)
            # input()
            

            # if plot > 1: # NOTE Implement Plotting later
            #     fH21d = np.zeros((nvx, nx)).T
            #     for k in range(0, nx - 1):
            #         fH21d[k] = np.dot(Vr2pidVr, fH2G[k])
            #     plt.plot(vx, fH21d[0, :], label="0", color="b", linewidth=2)
            #     for i in range(1, nx):
            #         if np.any(fH21d[i, :] > 0.9):
            #             plt.plot(vx, fH21d[i, :], label=str(i), color=(i % 8) + 2, linewidth=2)

            #     plt.title(str(igen) + ' Generation ' + 'HH')
            #     plt.xlabel('vx')
            #     plt.ylabel('fH21d')
            #     plt.ylim(0, np.max(fH21d))
            #     plt.legend(title="Curve Index")
            #     plt.show()
            #     if debug > 0:
            #         input()

            # Add result to total neutral distribution function
            fH2 = fH2 + fH2G
            nH2 = nH2 + NH2G[:,igen]
            # print("fH2", fH2.T)
            # print("nH2", nH2.T)
            # input()

            # Compute 'generation error': Delta_nH2G=max(NH2G(*,igen)/max(nH2))
            # and decide if another generation should be computed
            Delta_nH2G = np.max(NH2G[:,igen]/np.max(nH2))
            if (Delta_nH2G < truncate) or (fH2_iterate and (Delta_nH2G < 0.003*Delta_nH2s)):
                # If fH2 'seed' is being iterated, then do another generation until the 'generation error'
                # is less than 0.003 times the 'seed error' or is less than TRUNCATE
                break
            # If fH2 'seed' is NOT being iterated, then do another generation unitl the 'generation error'
            # is less than parameter TRUNCATE


# fH2_done #########################################################################################################################################################################

        # if plot > 0: #NOTE Implement Plotting Later
        #     plt.figure()

        #     # Plot the first curve with logarithmic y-axis
        #     plt.plot(x, NH2G[:, 0], label="0", color="b", linewidth=2)

        #     # Plot additional curves up to generation igen
        #     for i in range(igen + 1):
        #         plt.plot(x, NH2G[:, i], label=str(i), color=(i % 8) + 2, linewidth=2)
            
        #     # Set plot attributes
        #     plt.title(_HH + ' Density by Generation')
        #     plt.xlabel('x (m)')
        #     plt.ylabel('Density (m⁻³)')
        #     plt.yscale('log')
            
        #     plt.legend(title="Generation")

        # Compute H2 density profile
        for k in range(0, nx):
            nH2[k] = np.sum(Vr2pidVr*(fH2[:,:,k] @ dVx))
        # print("nH2", nH2)
        # input()

        # GammaxH2 - particle flux in x direction
        for k in range(0, nx):
            GammaxH2[k] = Vth*np.sum(Vr2pidVr*(fH2[:,:,k] @ (vx*dVx)))
        # print("GammaxH2", GammaxH2)
        # input()

        # VxH2 - x velocity
        VxH2 = GammaxH2 / nH2
        _VxH2 = VxH2 / Vth 
        # print("VxH2", VxH2)
        # print("_VxH2", _VxH2)
        # input()

        # magnitude of random velocity at each mesh point 
        vr2vx2_ran = np.zeros((nvr, nvx, nx))
        for i in range(0, nvr):
            for k in range(0, nx):
                vr2vx2_ran[i,:,k] = vr[i]**2 + (vx - _VxH2[k])**2
        # print("vr2vx2_ran", vr2vx2_ran.T)
        # input()

        # pH2 - pressure 
        for k in range(0, nx):
            pH2[k] = ((2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*fH2[:,:,k]) @ dVx)))/(3*CONST.Q)
        # print("pH2", pH2)
        # input()

        #TH2 - temperature 
        TH2 = pH2/nH2
        # print("TH2", TH2)
        # input()

        # Compute NuDis - Dissociation frequency 
        NuDis = n*np.sum(sigv[:,7:11], 1) # the indexing could be wrong here
        # print("NuDis", NuDis)
        # input()
        
        # Compute NuE (assume np=ne) - Energy equilibration frequency H(+) <-> H2(+)
        NuE = (7.7e-7*n*1.0e-6)/(np.sqrt(mu)*(Ti**1.5)) 
        # print("NuE", NuE)
        # input()
        
        # Compute H2(+) density profile
        nHP = (nH2*n*sigv[:,1])/(NuDis + NuLoss)
        # print("nHP", nHP)
        # input()

        # Compute THP - temperature of molecular ions
        THP = (Ti*NuE)/(NuE + NuDis + NuLoss)
        # print("THP", THP)
        # input()
        
        if fH2_iterate:
            # Compute 'seed error': Delta_nH2s=(|nH2s-nH2|)/max(nH2) 
            # If Delta_nH2s is greater than 10*truncate then iterate fH2

            Delta_nH2s = np.max(np.abs(nH2s - nH2))/np.max(nH2)
            if Delta_nH2s > 10*truncate:
                do_fH2_iterate = True
    # print("Delta_nH2s", Delta_nH2s)
    # input() 
    
    # Update Swall_sum using last generation
    Swall = np.zeros((nvr, nvx, nx))
    if np.sum(gamma_wall) > 0:
        for k in range(0, nx):
            Swall[:,:,k] = fw_hat * np.sum(Vr2pidVr*((gamma_wall[:,:,k]*fH2G[:,:,k]) @ dVx))
            Swall_sum = Swall_sum + Swall 
    # print("Swall", Swall.T)
    # print("Swall_sum", Swall_sum.T)
    # input()
    
    # Update Beta_CX_sum using last generation
    Beta_CX = np.zeros((nvr, nvx, nx))
    if H2_HP_CX:
        if debrief > 1:
            print(prompt, 'Computing Beta_CX')
        if Simple_CX:
            # Option (B): Compute charge exchange source with assumption that CX source neutrals have
            # molecular ion distribution function
            for k in range(0, nx):
                Beta_CX[:,:,k] = fHp_hat[:,:,k]*np.sum(Vr2pidVr*(Alpha_CX[:,:,k]*fH2G[:,:,k] @ dVx))
                # print(fH2G[:,:,k].T)
                # input()
        else: #NOTE Not Tested
            # Option (A): Compute charge exchange source using fH2 and vr x sigma x v_v at each velocity mesh point
            for k in range(0, nx):
                Work[:] = fH2G[:,:,k]
                Beta_CX[:,:,k] = nHP[k]*fHp_hat[:,:,k]*(SIG_CX @ Work)
        Beta_CX_sum = Beta_CX_sum + Beta_CX
    # print("Beta_CX", Beta_CX.T)
    # print("Beta_CX_sum", Beta_CX_sum.T)
    # input()

    # Update MH2_*_sum using last generation
    MH2_H2 = np.zeros((nvr,nvx,nx))
    MH2_P = np.zeros((nvr,nvx,nx))
    MH2_H = np.zeros((nvr,nvx,nx))
    OmegaM = np.zeros((nvr,nvx,nx))
    if H2_H2_EL or H2_P_EL or H2_H_EL: 
        # Compute VxH2G, TH2G
        for k in range(0, nx):
            VxH2G[k] = Vth*np.sum(Vr2pidVr*(fH2G[:,:,k] @ (vx*dVx))) / NH2G[k,igen]
            for i in range(0, nvr):
                vr2vx2_ran2[i,:] = vr[i]**2 + (vx - VxH2G[k]/Vth)**2
            TH2G[k] = (2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran2*fH2G[:,:,k]) @ dVx))/(3*CONST.Q*NH2G[k,igen])
        # print("TH2G", TH2G)
        # input()

        if H2_H2_EL:
            if debrief > 1: 
                print(prompt, 'Computing MH2_H2')
            # Compute MH2_H2
            vx_shift = VxH2G
            Tmaxwell = np.copy(TH2G)
            mol = 2
            Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                      nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
            for k in range(0, nx):
                MH2_H2[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen]
                OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_H2[k]*MH2_H2[:,:,k]
            MH2_H2_sum = MH2_H2_sum + MH2_H2
            # print("MH2_H2_sum", MH2_H2_sum.T)
            # input()

        if H2_P_EL:
            if debrief > 1:
                print(prompt, 'Computing MH2_P')
            # Compute MH2_P
            vx_shift = (2*VxH2G + vxi)/3
            Tmaxwell = TH2G + (4/9)*(Ti - TH2G + mu*CONST.H_MASS*((vxi - VxH2G)**2) / (6*CONST.Q))
            mol = 2
            Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                    nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                    Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
            for k in range(0, nx):
                MH2_P[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen]
                OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_P[k]*MH2_P[:,:,k]
            MH2_P_sum = MH2_P_sum + MH2_P
            # print("MH2_P_sum", MH2_P_sum.T)
            # input()

            if H2_H_EL: #NOTE Not Tested Yet
                if debrief > 1:
                    print(prompt, 'Computing MH2_H')
                # Compute MH2_H
                vx_shift = (2*VxH2G + VxH)/3
                Tmaxwell = TH2G + (4/9)*(TH - TH2G + mu*CONST.H_MASS*((VxH - VxH2G)**2)/(6*CONST.Q))
                mol = 2
                Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                        nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                        Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
                for k in range(0, nx):
                    MH2_H[:,:,k] = Maxwell[:,:,k]*NH2G[k,igen]
                    OmegaM[:,:,k] = OmegaM[:,:,k] + Omega_H2_H[k]*MH2_H[:,:,k]
                MH2_H_sum = MH2_H_sum + MH2_H
                # print("MH2_H_sum", MH2_H_sum.T)
                # input()
        
     # Compute remaining moments
     #NOTE Why are these seperate loops?
    # piH2_xx
    for k in range(0, nx):
        piH2_xx[k] = (((2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*(fH2[:,:,k] @ (dVx*(vx - _VxH2[k])**2))))/CONST.Q) - pH2[k]
    # piH2_yy
    for k in range(0, nx):
        piH2_yy[k] = (((2*mu*CONST.H_MASS)*Vth2*0.5*np.sum((Vr2pidVr*(vr**2))*(fH2[:,:,k] @ dVx)))/CONST.Q) - pH2[k]
    # piH2_zz 
    piH2_zz = copy.copy(piH2_yy) 
    # qxH2
    for k in range(0, nx):
        qxH2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*fH2[:,:,k]) @ (dVx*(vx - _VxH2[k]))))
    # print("piH2_xx", piH2_xx)
    # print("piH2_yy", piH2_yy)
    # print("qxH2", qxH2)
    # input()
    
    # C = RHS of Boltzman equation for total fH2
    for k in range(0, nx):
        C = Vth*((fw_hat[:,:]*SH2[k]/Vth) + Swall_sum[:,:,k] + Beta_CX_sum[:,:,k] - (alpha_c[:,:,k]*fH2[:,:,k]) + \
                   (Omega_H2_P[k]*MH2_P_sum[:,:,k]) + (Omega_H2_H[k]*MH2_H_sum[:,:,k]) + (Omega_H2_H2[k]*MH2_H2_sum[:,:,k]))
        # print("C", C.T)
        # input()
        QH2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2_ran[:,:,k]*C) @ dVx))
        RxH2[k] = (2*mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(C @ (dVx*(vx - _VxH2[k]))))
        Sloss[k] = -np.sum(Vr2pidVr*(C @ dVx)) + SH2[k]
        WallH2[k] = np.sum(Vr2pidVr*((gamma_wall[:,:,k]*fH2[:,:,k]) @ dVx))

        if H2_H_EL:
            CH2_H = Vth*Omega_H2_H[k]*(MH2_H_sum[:,:,k] - fH2[:,:,k])
            RxH_H2[k] = (2*mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CH2_H @ (dVx*(vx - _VxH2[k]))))
            # print("rxH_H2", MH2_H_sum[:,:,k].T)
            # print("b", fH2[:,:,k].T)
            # input()
            EH_H2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CH2_H) @ dVx))

        if H2_P_EL:
            CH2_P = Vth*Omega_H2_P[k]*(MH2_P_sum[:,:,k] - fH2[:,:,k])
            RxP_H2[k] = (2*mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CH2_P @ (dVx*(vx - _VxH2[k]))))
            EP_H2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CH2_P) @ dVx))

        if H2_HP_CX:
            CH2_HP_CX = Vth*(Beta_CX_sum[:,:,k] - Alpha_CX[:,:,k]*fH2[:,:,k])
            RxH2CX[k] = (2*mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CH2_HP_CX @ (dVx*(vx - _VxH2[k]))))
            EH2CX[k] = 0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CH2_HP_CX) @ dVx))

        CW_H2 = Vth*(Swall_sum[:,:,k] - gamma_wall[:,:,k]*fH2[:,:,k])
        RxW_H2[k] = (2*mu*CONST.H_MASS)*Vth*np.sum(Vr2pidVr*(CW_H2 @ (dVx*(vx - _VxH2[k]))))
        EW_H2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*CW_H2) @ dVx))

        if H2_H2_EL:
            CH2_H2 = Vth*Omega_H2_H2[k]*(MH2_H2_sum[:,:,k] - fH2[:,:,k])
            for i in range(0, nvr):
                vr2_2vx_ran2[i,:] = vr[i]**2 - 2*((vx - _VxH2[k])**2)
            Epara_PerpH2_H2[k] = -0.5*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((vr2_2vx_ran2*CH2_H2) @ dVx))
    # print("QH2", QH2)
    # print("RxH2", RxH2)
    # print("Sloss", Sloss)
    # print("WallH2", WallH2)
    # input()
    # print("RxH_H2", RxH_H2)
    # print("EH_H2", EH_H2)
    # input()
    # print("RxP_H2", RxP_H2)
    # print("EP_H2", EP_H2)
    # input()
    # print("RxH2CX", RxH2CX)
    # print("EH2CX", EH2CX)
    # input()
    # print("RxW_H2", RxW_H2)
    # print("EW_H2", EW_H2)
    # input()
    # print("Epara_PerpH2_H2", Epara_PerpH2_H2)
    # input()

    # qxH2_total
    qxH2_total = (0.5*nH2*(2*mu*CONST.H_MASS)*VxH2*VxH2 + 2.5*pH2*CONST.Q)*VxH2 + CONST.Q*piH2_xx*VxH2 + qxH2
    # print("qxH2_total", qxH2_total)
    # input()

    # QH2_total
    QH2_total = QH2 + RxH2*VxH2 - 0.5*(2*mu*CONST.H_MASS)*(Sloss - SH2)*VxH2*VxH2
    # print("QH2_total", QH2_total)
    # input()

    # Albedo
    gammax_plus = Vth*np.sum(Vr2pidVr*(fH2[:,i_p[:,0],0] @ (vx[i_p[:,0]]*dVx[i_p[:,0]]))) #NOTE Had to reference i_p and i_n in a weird way, fix how they are called in the first place
    gammax_minus = Vth*np.sum(Vr2pidVr*(fH2[:,i_n[:,0],0] @ (vx[i_n[:,0]]*dVx[i_n[:,0]]))) #This is awful and should not be allowed to remain
    if np.abs(gammax_plus) > 0:
        AlbedoH2 = -gammax_minus/gammax_plus
    # print("AlbedoH2", AlbedoH2)
    # input()

    # Compute Mesh Errors
    mesh_error = np.zeros((nvr,nvx,nx))
    max_mesh_error = 0.0
    min_mesh_error = 0.0
    mtest = 5
    moment_error = np.zeros((nx,mtest))
    max_moment_error = np.zeros(mtest)
    C_error = np.zeros(nx)
    CX_error = np.zeros(nx)
    Wall_error = np.zeros(nx)
    H2_H2_error = np.zeros((nx, 3))
    H2_H_error = np.zeros((nx, 3))
    H2_P_error = np.zeros((nx, 3))
    max_H2_H2_error = np.zeros(3)
    max_H2_H_error = np.zeros(3)
    max_H2_P_error = np.zeros(3)

    if Compute_Errors: #NOTE Simplify these with lambda functions, not norm and error calculation are repeated and unneccesary

        if debrief > 1:
            print(prompt, 'Computing Collision Operator, Mesh, and Moment Normalized Errors')

        Sloss2 = Vth*Alpha_Loss*nH2 
        for k in range(0, nx):
            C_error[k] = np.abs(Sloss[k] - Sloss2[k])/np.max(np.abs(np.array([Sloss[k], Sloss2[k]])))
        # print("C_error", C_error)
        # input()

        # Test conservation of particles for charge exchange operator
        if H2_HP_CX:
            for k in range(0, nx):
                CX_A = np.sum(Vr2pidVr*((Alpha_CX[:,:,k]*fH2[:,:,k]) @ dVx))
                CX_B = np.sum(Vr2pidVr*(Beta_CX_sum[:,:,k] @ dVx))
                CX_error[k] = np.abs(CX_A - CX_B)/np.max(np.abs(np.array([CX_A, CX_B])))
            # print("CX_error", CX_error)
            # input()

        # Test conservation of particles for wall collision operator
        if np.sum(PipeDia) > 0: #NOTE Not Tested Yet
            for k in range(0, nx):
                Wall_A = WallH2[k]
                Wall_B = np.sum(Vr2pidVr*(Swall_sum[:,:,k] @ dVx))
                if np.max(np.abs(np.array([Wall_A, Wall_B]))) > 0:
                    Wall_error[k] = np.abs(Wall_A - Wall_B)/np.max(np.abs(np.array([Wall_A, Wall_B])))
            # print("Wall_error", Wall_error)
            # input()

        # Test conservation of particles, x momentum, and total energy of elastic collision operators
        for m in range(0, 3):
            for k in range(0, nx):
                if m < 2:
                    TfH2 = np.sum(Vr2pidVr*(fH2[:,:,k] @ (dVx*(vx**m))))
                else:
                    TfH2 = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*fH2[:,:,k]) @ dVx))

                if H2_H2_EL:
                    if m < 2:
                        TH2_H2 = np.sum(Vr2pidVr*(MH2_H2_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH2_H2 = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH2_H2_sum[:,:,k]) @ dVx))
                    H2_H2_error[k,m] = np.abs(TfH2 - TH2_H2)/np.max(np.abs(np.array([TfH2, TH2_H2])))
                
                if H2_H_EL:
                    if m < 2:
                        TH2_H = np.sum(Vr2pidVr*(MH2_H_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH2_H = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH2_H_sum[:,:,k]) @ dVx))
                    H2_H_error[k,m] = np.abs(TfH2 - TH2_H)/np.max(np.abs(np.array([TfH2, TH2_H])))

                if H2_P_EL:
                    if m < 2:
                        TH2_P = np.sum(Vr2pidVr*(MH2_P_sum[:,:,k] @ (dVx*(vx**m))))
                    else:
                        TH2_P = np.sum(Vr2pidVr*((vr2vx2[:,:,k]*MH2_P_sum[:,:,k]) @ dVx))
                    H2_P_error[k,m] = np.abs(TfH2 - TH2_P)/np.max(np.abs(np.array([TfH2, TH2_P])))

            max_H2_H2_error[m] = np.max(H2_H2_error[:,m])
            max_H2_H_error[m] = np.max(H2_H_error[:,m])
            max_H2_P_error[m] = np.max(H2_P_error[:,m])
            # print("max_H2_H2_error", max_H2_H2_error)
            # print("max_H2_H_error", max_H2_H_error)
            # print("max_H2_P_error", max_H2_P_error)
            # input()

        if CI_Test:
            minRx = 1.0e-6
            minEpara_perp = 1.0e-6

            # Compute Momentum transfer rate via full collision integrals for charge exchange and mixed elastic scattering
            # Then compute error between this and actual momentum transfer resulting from CX and BKG (elastic) models

            if H2_HP_CX: # H2(+) -> H2 charge exchange momentum transfer via full collision integral
                print(prompt, 'Computing H2(+) -> H2 Charge Exchange Momentum Transfer')
                _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
                _Sig[:] = (v_v*sigma_cx_hh(v_v2*(CONST.H_MASS*Vth2/CONST.Q))).reshape(_Sig.shape, order='F')
                SIG_VX_CX = np.zeros((nvr*nvx,nvr*nvx))
                SIG_VX_CX[:] = (Vr2pidVrdVx*vx_vx*((_Sig @ dTheta).reshape(vx_vx.shape, order='F'))).reshape(SIG_VX_CX.shape, order='F')
                alpha_vx_cx = np.zeros((nvr,nvx,nx))

                for k in range(0, nx):
                    Work[:] = (nHP[k]*fHp_hat[:,:,k]).reshape(Work.shape, order='F')
                    alpha_vx_cx[:,:,k] = (SIG_VX_CX @ Work).reshape(alpha_vx_cx[:,:,k].shape, order='F')

                for k in range(0, nx):
                    RxCI_CX[k] = -(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((alpha_vx_cx[:,:,k]*fH2[:,:,k]) @ dVx))

                norm = np.max(np.abs(np.array([RxH2CX, RxCI_CX])))
                for k in range(0, nx):
                    CI_CX_error[k] = np.abs(RxH2CX[k] - RxCI_CX[k])/norm

                print(prompt,'Maximum normalized momentum transfer error in CX collision operator: ', sval(np.max(CI_CX_error)))

            if H2_P_EL: # P -> H2 momentum transfer via full collision integral
                for k in range(0, nx):
                    RxCI_P_H2[k] = -(1/3)*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((Alpha_H2_P[:,:,k] * fH2[:,:,k]) @ dVx))

                norm = np.max(np.abs(np.array([RxP_H2, RxCI_P_H2])))
                for k in range(0, nx):
                    CI_P_H2_error[k] = np.abs(RxP_H2[k] - RxCI_P_H2[k])/norm 

                print(prompt, 'Maximum normalized momentum transfer error in P -> H2 elastic BKG collision operator: ', sval(np.max(CI_P_H2_error)))
            
            if H2_H_EL: # H -> H2 momentum transfer via full collision integral
                for k in range(0, nx):
                    RxCI_H_H2[k] = -(1/3)*(2*mu*CONST.H_MASS)*Vth2*np.sum(Vr2pidVr*((Alpha_H2_H[:,:,k]*fH2[:,:,k]) @ dVx))
                
                norm = np.max(np.abs(np.array([RxH_H2, RxCI_H_H2])))
                for k in range(0, nx):
                    CI_H_H2_error[k] = np.abs(RxH_H2[k] - RxCI_H_H2[k])/norm
                
                print(prompt, 'Maximum normalized momentum transfer error in H -> H2 elastic BKG collision operator: ', sval(np.max(CI_H_H2_error)))
            
            if H2_H2_EL: # H2 -> H2 perp/parallel energy transfer via full collision integral
                for k in range(0, nx):
                    Work[:] = fH2[:,:,k].reshape(Work.shape, order='F')
                    Alpha_H2_H2[:] = (SIG_H2_H2 @ Work).reshape(Alpha_H2_H2.shape, order='F')
                    Epara_Perp_CI[k] = 0.5*(2*mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((Alpha_H2_H2*fH2[:,:,k]) @ dVx)) 
                
                norm = np.max(np.abs(np.array([Epara_PerpH2_H2, Epara_Perp_CI])))
                for k in range(0, nx):
                    CI_H2_H2_error[k] = np.abs(Epara_PerpH2_H2[k] - Epara_Perp_CI[k])/norm 
                
                print(prompt, 'Maximum normalized perp/parallel energy transfer error in H2 -> H2 elastic BKG collision operator: ', sval(np.max(CI_H2_H2_error)))
        
        # Mesh Point Error based on fH2 satisfying Boltzmann equation
        T1 = np.zeros((nvr,nvx,nx))
        T2 = np.zeros((nvr,nvx,nx))
        T3 = np.zeros((nvr,nvx,nx))
        T4 = np.zeros((nvr,nvx,nx))
        T5 = np.zeros((nvr,nvx,nx))
        T6 = np.zeros((nvr,nvx,nx))
        for k in range(0, nx-1):
            for j in range(0, nvx):
                T1[:,j,k] = 2*vx[j]*(fH2[:,j,k+1] - fH2[:,j,k])/(x[k+1] - x[k]) 
            T2[:,:,k] = fw_hat[:,:]*(SH2[k+1] + SH2[k])/Vth
            T3[:,:,k] = Beta_CX_sum[:,:,k+1] + Beta_CX_sum[:,:,k]
            T4[:,:,k] = alpha_c[:,:,k+1]*fH2[:,:,k+1] + alpha_c[:,:,k]*fH2[:,:,k]
            T5[:,:,k] = Omega_H2_P[k+1]*MH2_P_sum[:,:,k+1] + Omega_H2_H[k+1]*MH2_H_sum[:,:,k+1] + Omega_H2_H2[k+1]*MH2_H2_sum[:,:,k+1] + \
                    Omega_H2_P[k]*MH2_P_sum[:,:,k] + Omega_H2_H[k]*MH2_H_sum[:,:,k] + Omega_H2_H2[k]*MH2_H2_sum[:,:,k]
            T6[:,:,k] = Swall_sum[:,:,k+1] + Swall_sum[:,:,k]
            mesh_error[:,:,k] = np.abs(T1[:,:,k] - T2[:,:,k] - T3[:,:,k] + T4[:,:,k] - T5[:,:,k] - T6[:,:,k])/ \
                                np.max(np.abs(np.array([T1[:,:,k], T2[:,:,k], T3[:,:,k], T4[:,:,k], T5[:,:,k], T6[:,:,k]])))
        ave_mesh_error = np.sum(mesh_error) / np.size(mesh_error)
        max_mesh_error = np.max(mesh_error)
        min_mesh_error = np.min(mesh_error[:,:,0:nx-1])

        # Moment Error
        for m in range(0, mtest):
            for k in range(0, nx - 1):
                MT1 = np.sum(Vr2pidVr*(T1[:,:,k] @ (dVx*(vx**m))))
                MT2 = np.sum(Vr2pidVr*(T2[:,:,k] @ (dVx*(vx**m))))
                MT3 = np.sum(Vr2pidVr*(T3[:,:,k] @ (dVx*(vx**m))))
                MT4 = np.sum(Vr2pidVr*(T4[:,:,k] @ (dVx*(vx**m))))
                MT5 = np.sum(Vr2pidVr*(T5[:,:,k] @ (dVx*(vx**m))))
                MT6 = np.sum(Vr2pidVr*(T6[:,:,k] @ (dVx*(vx**m))))
                #NOTE This is correct for the original code, but is it correct mathematically?
                moment_error[k,m] = np.abs(MT1 - MT2 - MT3 + MT4 - MT5 - MT6) / np.max(np.abs(np.array([MT1, MT2, MT3, MT4, MT5, MT6])))
            max_moment_error[m] = np.max(moment_error[:,m])

        # Compute error in qxH2_total
        # qxH2_total2 total neutral heat flux profile (watts m^-2)
        #    This is the total heat flux transported by the neutrals
        #    computed in a different way from:
        # 
        #    qxH2_total2(k)=vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH2(*,*,k))#(Vx*dVx)))*0.5*(2*mu*mH)
        # 
        #    This should agree with qxH2_total if the definitions of nH2, pH2, piH2_xx,
        #    TH2, VxH2, and qxH2 are coded correctly.
        qxH2_total2 = np.zeros(nx)
        for k in range(0, nx):
            qxH2_total2[k] = 0.5*(2*mu*CONST.H_MASS)*Vth3*np.sum(Vr2pidVr*((vr2vx2[:,:,k]*fH2[:,:,k]) @ (vx*dVx)))
        qxH2_total_error = np.abs(qxH2_total - qxH2_total2) / np.max(np.abs(np.array([qxH2_total, qxH2_total2])))

        # Compute error in QH2_total
        Q1 = np.zeros(nx)
        Q2 = np.zeros(nx)
        QH2_total_error = np.zeros(nx)
        for k in range(0, nx-1):
            Q1[k] = (qxH2_total[k+1] - qxH2_total[k]) / (x[k+1] - x[k])
            Q2[k] = 0.5*(QH2_total[k+1] + QH2_total[k])
        QH2_total_error = np.abs(Q1 - Q2)/np.max(np.abs(np.array([Q1, Q2])))

        if debrief > 0:
            print(prompt, 'Maximum particle convervation error of total collision operator: ', sval(np.max(C_error)))
            print(prompt, 'Maximum H2_HP_CX particle convervation error: ', sval(np.max(CX_error)))
            print(prompt, 'Maximum H2_Wall particle convervation error: ', sval(np.max(Wall_error)))
            print(prompt, 'Maximum H2_H2_EL particle conservation error: ', sval(max_H2_H2_error[0]))
            print(prompt, 'Maximum H2_H2_EL x-momentum conservation error: ', sval(max_H2_H2_error[1]))
            print(prompt, 'Maximum H2_H2_EL total energy conservation error: ', sval(max_H2_H2_error[2]))
            print(prompt, 'Maximum H2_H_EL  particle conservation error: ', sval(max_H2_H_error[0]))
            print(prompt, 'Maximum H2_P_EL  particle conservation error: ', sval(max_H2_P_error[0]))
            print(prompt, 'Average mesh_error =', ave_mesh_error)
            print(prompt, 'Maximum mesh_error =', max_mesh_error)
            for m in range(0, 5):
                print(prompt, 'Maximum fH2 vx^', sval(m), ' moment error: ', sval(max_moment_error[m]))
            print(prompt, 'Maximum qxH2_total error =', np.max(qxH2_total_error))
            print(prompt, 'Maximum QH2_total error =', np.max(QH2_total_error))
            if debug > 0:
                input()
    
    # NOTE Commented Plotting elements, re-add once program is working
    # mid1 = locate(x, 0.7 * (np.max(x) + np.min(x)) / 2)
    # mid2 = locate(x, 0.85 * (np.max(x) + np.min(x)) / 2)
    # mid3 = locate(x, (np.max(x) + np.min(x)) / 2)
    # mid4 = locate(x, 1.15 * (np.max(x) + np.min(x)) / 2)
    # mid5 = locate(x, 1.3 * (np.max(x) + np.min(x)) / 2)
    # mid6 = locate(x, 1.45 * (np.max(x) + np.min(x)) / 2)

    # if plot > 1:
    #     fH21d = np.zeros((nvx, nx)).T
    #     for k in range(0, nx - 1):
    #         fH21d[k] = np.dot(Vr2pidVr, fH2[k])
    #     plt.figure()
    #     plt.plot(vx, fH21d[:, 0], label="0", linewidth=2)
    #     for i in range(nx):
    #         plt.plot(vx, fH21d[:, i], label=str(i), color=(i % 6) + 2, linewidth=2)

    #     plt.title(_HH + ' Velocity Distribution Function: fH2(Vx)')
    #     plt.xlabel('Vx/Vth')
    #     plt.legend()
    #     plt.show()
    #     if pause:
    #         # press_return 
    #         return
    
    # if plot > 0:
    #     data = np.array([nH, n, nHP, nH2])
    #     jp = np.argwhere(data > 0).shape[0]
    #     yrange = np.array([np.min(data[jp]), np.max(data[jp])])

    #     plt.figure()
    #     plt.plot(x, nH, label='nH', color='b')
    #     plt.plot(x, n, label='n', color='g')
    #     plt.plot(x, nH2, label='nH2', color='r')
    #     plt.plot(x, nHP, label='nHP', color='c')


    #     plt.title('Density Profiles')
    #     plt.xlabel('x (meters)')
    #     plt.ylabel('m⁻³')
    #     plt.yscale('log')
    #     plt.ylim(yrange)
    #     plt.legend()

    #     plt.annotate(_H, (x[mid1], 1.2 * nH[mid1]), color='b')
    #     plt.annotate('e⁻', (x[mid2], 1.2 * n[mid2]), color='g')
    #     plt.annotate(_HH, (x[mid3], 1.2 * nH2[mid3]), color='r')
    #     plt.annotate(_Hp, (x[mid4], 1.2 * nHP[mid4]), color='c')
    #     plt.show()

    #     if pause:
    #         # press_return 
    #         return 
    
    # if plot > 0:
    #     data = np.array([TH, Te, THP, TH2])
    #     jp = np.argwhere(data > 0)
    #     yrange = np.array([np.min(data[jp]), np.max(data[jp])])

    #     plt.figure()
    #     plt.plot(x, TH, label='TH', color='b')
    #     plt.plot(x, Te, label='Te', color='g')
    #     plt.plot(x, TH2, label='TH2', color='r')
    #     plt.plot(x, THP, label='THP', color='c')

    #     plt.title('Temperature Profiles')
    #     plt.xlabel('x (meters)')
    #     plt.ylabel('eV')
    #     plt.yscale('log')
    #     plt.ylim(yrange)
    #     plt.legend()

    #     plt.annotate(_H, (x[mid1], 1.2 * TH[mid1]), color='b')
    #     plt.annotate('e⁻', (x[mid2], 1.2 * Te[mid2]), color='g')
    #     plt.annotate(_HH, (x[mid3], 1.2 * TH2[mid3]), color='r')
    #     plt.annotate(_Hp, (x[mid4], 1.2 * THP[mid4]), color='c')
    #     plt.show()
    #     if pause:
    #         #press_return 
    #         return


    if Compute_H_Source:
        if debrief > 1:
            print(prompt, 'Computing Velocity Distributions of H products...')
        # Set Normalized Franck-Condon Velocity Distributions for reactions R2, R3, R4, R5, R6, R7, R8, R10

        # Make lookup table to select reaction Rn in SFCn
        #   Rn=2 3 4 5 6 7 8   10
        nFC = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 0, 7])
        SFCn = np.zeros((nvr,nvx,nx,8))
        Eave = np.zeros((nx, 8))
        Emax = np.zeros((nx, 8))
        Emin = np.zeros((nx, 8))

        # Reaction R2: e + H2 -> e + H(1s) + H(1s)
        ii = nFC[2]
        Eave[:,ii] = 3.0 
        Emax[:,ii] = 4.25
        Emin[:,ii] = 2

        # Reaction R3: e + H2 -> e + H(1s) + H*(2s)
        ii = nFC[3]
        Eave[:,ii] = 0.3
        Emax[:,ii] = 0.55
        Emin[:,ii] = 0.0 

        # Reaction R4:  e + H2 -> e + H(+) + H(1s) + e
        ii = nFC[4]
        Ee = 3*Te/2     # Note the FC energy depends on electron energy
        kk = np.argwhere(Ee <= 26.0)
        if kk.size > 0:
            Eave[kk,ii] = 0.25
        kk = np.argwhere((Ee > 26.0) & (Ee <= 41.6))
        if kk.size > 0:
            Eave[kk,ii] = 0.5*(Ee[kk] - 26)
            Eave[kk,ii] = np.maximum(Eave[kk, ii], 0.25)
        kk = np.argwhere(Ee > 41.6)
        if kk.size > 0:
            Eave[kk,ii] = 7.8
        Emax[:,ii] = 1.5*Eave[:,ii]   # Note the max/min values here are a guess
        Emin[:,ii] = 0.5*Eave[:,ii]   # Note the max/min values here are a guess

        # Reaction R5: e + H2 -> e + H*(2p) + H*(2s)
        ii = nFC[5]
        Eave[:,ii] = 4.85 
        Emax[:,ii] = 5.85
        Emin[:,ii] = 2.85

        # Reaction R6: e + H2 -> e + H(1s) + H*(n=3)
        ii = nFC[6]
        Eave[:,ii] = 2.5 
        Emax[:,ii] = 3.75 
        Emin[:,ii] = 1.25 

        # Reaction R7: e + H2(+) -> e + H(+) + H(1s)
        ii = nFC[7]
        Eave[:,ii] = 4.3   
        Emax[:,ii] = 4.3 + 2.1     # Note the max/min values here are a guess
        Emin[:,ii] = 4.3 - 2.1     # Note the max/min values here are a guess

        # Reaction R8: e + H2(+) -> e + H(+) + H*(n=2)
        ii = nFC[8]
        Eave[:,ii] = 1.5 
        Emax[:,ii] = 1.5 + 0.75     # Note the max/min values here are a guess
        Emin[:,ii] = 1.5 - 0.75     # Note the max/min values here are a guess 

        # Reaction R10: e + H2(+) -> H(1s) + H*(n>=2)
        ii = nFC[10]

        # Compute relative cross-sections for populating a specific n level for reaction R10
        # (see page 62 in Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987)
        #   n=2   3    4    5    6
        R10rel = np.array([0.1, 0.45, 0.22, 0.12, 0.069])
        for k in range(7, 11): 
            R10rel = np.append(R10rel, 10/(k**3))
        En = 13.58/((2 + np.arange(9))**2) # Energy of Levels

        for k in range(0, nx):
            truncate_point = np.minimum(len(Ee), len(En))
            EHn = 0.5*(Ee[:truncate_point] - En[:truncate_point])*R10rel/np.sum(R10rel) #NOTE For some reason, this is how idl handles subtracting two arrays of different sizes
            EHn = np.maximum(EHn, 0)
            Eave[k,ii] = np.sum(EHn)
            Eave[k,ii] = np.maximum(Eave[k,ii], 0.25)
            Emax[k,ii] = 1.5*Eave[k,ii] # Note the max/min values here are a guess
            Emin[k,ii] = 0.5*Eave[k,ii] # Note the max/min values here are a guess
        # print("Eave", Eave.T)
        # print("Emax", Emax.T)
        # print("Emin", Emin.T)
        # input()
        
        # Set SFCn values for reactions R2, R3, R4, R5, R6, R7, R8, R10
        Vfc = np.zeros((nvr,nvx,nx))
        Tfc = np.zeros((nvr,nvx,nx))
        magV = np.sqrt(vr2vx2)
        _THP = np.zeros((nvr,nvx,nx))
        _TH2 = np.zeros((nvr,nvx,nx)) 
        for k in range(0, nx):
            _THP[:,:,k] = THP[k]/Tnorm
            _TH2[:,:,k] = TH2[k]/Tnorm 
        # print("_THP", _THP.T)
        # print("_TH2", _TH2.T)
        # input()
        
        # The following function is choosen to represent the velocity distribution of the
        # hydrogen products for a given reaction, accounting for the Franck-Condon energy
        # distribution and accounting for additional velocity spreading due to the finite
        # temperature of the molcules (neutral and ionic) prior to breakup:
        # 
        #     f(Vr,Vx) = exp( -0.5*mH*mu*(|v|-Vfc+0.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )

        #       	|v|=sqrt(Vr^2+Vx^2)
        #	        Tfc= Franck-Condon 'temperature' = (Emax-Emin)/4
        #	        Vfc= Franck-Condon  velocity = sqrt(2 Eave/mH/mu)
        #		    Tmol= temperature of H2 molecule (neutral or ionic)

        #    This function is isotropic in velocity space and can be written in terms
        #  of a distribution in particle speed, |v|, 

        #     f(|v|) = exp( -(|v|-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )
        #
        # with velocities normalized by vth and T normalized by Tnorm.

        #  Recognizing the the distribution in energy, f(E), and particle speed, f(|v|),
        #  are related by  f(E) dE = f(|v|) 2 pi v^2 dv, and that dE/dv = mH mu v,
        #  f(E) can be written as

        #     f(E) = f(|v|) 2 pi |v|/(mH mu) = const. |v| exp( -(|v|-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )

        # The function f(Vr,Vx) was chosen because it has has the following characteristics:

        # (1) For Tmol/2 << Tfc,  the peak in the v^2 times the energy distribution, can be found
        #    by finding the |v| where df(E)/d|v| =0

        #    df(E)/d|v|= 0 = 3v^2 exp() - 2(|v|-Vfc+1.5*Tfc/Vfc)/Tfc v^3 exp() 
        #                    2(|v|-Vfc+1.5*Tfc/Vfc)/Tfc |v| = 3
        #    which is satisfied when |v|=Vfc. Thus the energy-weighted energy distribution peaks
        #    at the velocity corresponding to the average Franck-Condon energy.

        # (2) for Tmol/2 >> Tfc ~ Vfc^2, the velocity distribution becomes

        #	f(|v|) = exp( -2(|v|-Vfc+1.5*Tfc/Vfc)^2/Tmol )

        #    which leads to a velocity distribution that approaches the molecular velocity
        #    distribution with the magnitude of the average velocity divided by 2. This
        #    is the appropriate situation for when the Franck-Condon energies are negligible
        #    relative to the thermal speed of the molecules.

        Rn = np.array([2, 3, 4, 5, 6, 7, 8, 10])
        for jRn in range(np.size(Rn)):
            ii = nFC[Rn[jRn]]
            # print("ii", ii)
            Tfc[0,0,:] = 0.25*(Emax[:,ii] - Emin[:,ii])/Tnorm # Franck-Condon 'effective temperature'
            Vfc[0,0,:] = np.sqrt(Eave[:,ii]/Tnorm) # Velocity corresponding to Franck-Condon 'mean evergy'
            for k in range(nx):
                Vfc[:,:,k] = Vfc[0,0,k]
                Tfc[:,:,k] = Tfc[0,0,k]

            if Rn[jRn] <= 6:
                # For R2-R6, the Franck-Condon 'mean energy' is taken equal to Eave
                #	   and the 'temperature' corresponds to the sum of the Franck-Condon 'temperature', Tfc,
                #          and the temperature of the H2 molecules, TH2. (Note: directed neutral molecule velocity
                #	   is not included and assumed to be small)
                arg = -(magV - Vfc + 1.5*Tfc/Vfc)**2 / (Tfc + 0.5*_TH2)
                SFCn[:,:,:,ii] = np.exp(np.maximum(arg, -80))
            else: 
            #   For R7, R8 and R10, the Franck-Condon 'mean energy' is taken equal to Eave
            #	   and the 'temperature' corresponds to the sum of the Franck-Condon 'temperature', Tfc,
            #          and the temperature of the H2(+) molecular ions, THP. (Note: directed molecular ion velocity
            #	   is not included and assumed to be small)    
                arg = -(magV - Vfc + 1.5*Tfc/Vfc)**2 / (Tfc + 0.5*_THP)
                SFCn[:,:,:,ii] = np.exp(np.maximum(arg, -80))
            
            # print("sfcn", SFCn[:,:,0,4].T)
            # input()

            for k in range(nx):
                SFCn[:,:,k,ii] = SFCn[:,:,k,ii] / (np.sum(Vr2pidVr*(SFCn[:,:,k,ii] @ dVx)))
        # print("SFCn", SFCn.T)
        # print("Tfc", Tfc.T)
        # print("Vfc", Vfc.T)
        # input()
        
        #NOTE Removed plotting for now, re-implement once program is working
        # nm = 3
        # if plot > 3:
        #     ii = nFC[Rn[0]]
        #     m = np.arange(nm) * (nx-1) / (nm-1)
        #     kp = 0

        #     for mm in range(nm):
        #         for jRn in range(len(Rn)):
        #             ii = nFC[Rn[jRn]]
        #             color = (kp % 6) + 1
        #             plt.contour(SFCn[ii, int(m[mm])], colors=f'C{color}', linestyles='dashed')

        #             plt.text(0.7, 0.8 - kp * 0.03, f'R{Rn[jRn]} Te:{Te[int(m[mm])]}', transform=plt.gca().transAxes, color=f'C{color}')
        #             kp = kp + 1
        #             plt.pause(0.1)  # Pause for a short while (adjust as needed)

        #     ii = nFC[Rn[0]]
        #     plt.contour(SFCn[:, :, 0, ii], colors=f'C{color}', linestyles='dashed')
        #     plt.title('SFCn - Franck-Condon Velocity Distributions')

        #     plt.show() # I just realized that sum of the indexing in the plotting code might be incorrect 
    
        # if plot > 2:
        #     m = np.arange(nm) * (nx - 1)/(nm - 1)
        #     kp = 0

        #     plt.figure()
        #     plt.semilogx(Eaxis, Eaxis * 0, linestyle='-', color='none')  # Dummy line for plotting
        #     plt.ylim([0, 1])
        #     plt.title('Energy Distribution of ' + _H + ' Products')
        #     plt.xlabel('Energy (eV)')

        #     for mm in range(len(Rn)):
        #         kp += 1
        #         plt.text(0.85, 0.92 - kp * 0.02, f'Te:{Te[int(m[mm])]:.4f}', fontsize=8)
        #         kp += 1
        #         for jRn in range(len(Rn)):
        #             ii = nFC[Rn[jRn]]
        #             color - jRn % len(Rn) + 1

        #             EFC = Eaxis * SFCn[int(m[mm]), i_p[0]] * VrVr4pidVr / dEaxis
        #             EFC /= np.max(EFC)

        #             plt.plot(Eaxis, EFC, color = f'C{color}')
        #             plt.text(0.87, 0.92 - kp * 0.02, f'R{Rn[jRn]}', fontsize=8, color=f'C{color}')
        #             kp += 1
        #     plt.show()
        #     if pause: 
        #         # press_return 
        #         return 
    
        # if plot > 0:
        #     kp = 0
        #     for jRn in range(0, np.size(Rn - 1)):
        #         ii = nFC[Rn[jRn]]
        #         Ebar = np.zeros(nx)
        #         for k in range(0, nx - 1):
        #             Ebar[k] = 0.5 * (mu * CONST.H_MASS) * Vth2 * np.sum(Vr2pidVr * np.dot(vr2vx2[k] * SFCn[ii, k], dVx)) / CONST.Q
        #         color (kp % np.size(Rn)) + 1
        #         if kp == 0:
        #             plt.figure()
        #             plt.plot(x, Ebar, label='R' + str(Rn[jRn]) + ': ' + _Rn[Rn[jRn]])
        #             plt.title('Average Energy of ' + _H + ', ' + _p + ' Products')
        #             plt.xlabel('x (meters)')
        #             plt.ylabel('eV')
        #             plt.yscale('log')
        #             plt.ylim([0.1, 100])
        #             plt.legend()
        #         else:
        #             plt.plot(x, Ebar, label='R' + str(Rn[jRn]) + ': ' + _Rn[Rn[jRn]], color=f'C{color}')
        #         kp += 1
        #     plt.show()
        #     if pause:
        #         # press_return 
        #         return 

        Vbar_Error = np.zeros(nx)
        if Compute_Errors:
            # Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti(x)/Tnorm)/sqrt(!pi)
            TFC = np.min(Eave[0,:]) + ((np.max(Eave[0,:]) - np.min(Eave[0,:]))*np.arange(0, nx) / (nx - 1))
            # vx_shift = TFC
            # vx_shift[:] = 0.0 
            vx_shift = np.zeros_like(TFC)
            Tmaxwell = TFC
            mol = 1
            Maxwell = create_shifted_maxwellian_include(vr,vx,Tnorm,vx_shift,Tmaxwell,shifted_Maxwellian_debug,mu,mol,
                                      nx,nvx,nvr,Vth,Vth2,Maxwell,vr2vx2_ran2,
                                      Vr2pidVr,dVx,vol,Vth_DeltaVx,Vx_DeltaVx,Vr_DeltaVr,vr2_2vx2_2D,jpa,jpb,jna,jnb)
            vbar_test = Vth*np.sqrt(vr2vx2[:,:,0])
            for k in range(0, nx):
                vbar = np.sum(Vr2pidVr*((vbar_test*Maxwell[:,:,k]) @ dVx))
                vbar_exact = 2*Vth*np.sqrt(TFC[k] / Tnorm) / np.sqrt(np.pi)
                Vbar_Error[k] = np.abs(vbar - vbar_exact) / vbar_exact
            if debrief > 0: 
                print(prompt, 'Maximum Vbar error over FC energy range = ', np.max(Vbar_Error))
            # input()

        # Compute atomic hydrogen source distribution function
        # using normalized FC source distributions SFCn
        fSH_calc = lambda k,x : sigv[k,x]*SFCn[:,:,k,nFC[x]]
        for k in range(0, nx):

            fSH[:,:,k] = n[k]*nH2[k]*(2*fSH_calc(k,2) + 2*fSH_calc(k,3) + fSH_calc(k,4) + 2*fSH_calc(k,5) + 2*fSH_calc(k,6))
            
            fSH[:,:,k] = fSH[:,:,k] + n[k]*nHP[k]*(fSH_calc(k,7) + fSH_calc(k,8) + 2*fSH_calc(k,10))
        # print("fSH", fSH.T)
        # input()

        # Compute total H and H(+) sources
        for k in range(nx):
            SH[k] = np.sum(Vr2pidVr*(fSH[:,:,k] @ dVx))
            SP[k] = n[k]*nH2[k]*sigv[k,4] + n[k]*nHP[k]*(sigv[k,7] + sigv[k,8] + 2*sigv[k,9])
        # print("SH", SH)
        # print("SP", SP)
        # input()

        # Compute total HP source
        SHP = n*nH2*sigv[:,1]
        # print("SHP", SHP.T)
        # input()

        # Compute energy distrobution of H source 
        for k in range(0, nx):
            ESH[:,k] = (Eaxis*fSH[:,i_p[0,0],k]*VrVr4pidVr) / dEaxis #NOTE Had to reference i_p and i_n in a weird way, fix how they are called in the first place
            ESH[:,k] = ESH[:,k] / np.max(ESH[:,k])
        # print("ESH", ESH.T)
        # input()
        
        #NOTE Removed plotting for now, re-implement once program is working
        # if plot > 2:
        #     fH21d = np.zeros((nvx, nx)).T
        #     for k in range(0, nx - 1):
        #         fH21d[k] = np.dot(Vr2pidVr, fSH[k])
        #     plt.figure()
        #     plt.plot(vx, fH21d[:, 0], label=_H + ' Source Velocity Distribution Function: fSH(Vx)')
        #     plt.xlabel('Vx/Vth')
        #     plt.ylabel('m⁻³ s⁻¹ dVx⁻¹')
        #     plt.title(_H + ' Source Velocity Distribution Function: fSH(Vx)')
        #     plt.yscale('log')

        #     for i in range(nx):
        #         plt.plot(vx, fH21d[:, i], label=f'fH21d[:, {i}]', color=f'C{(i % 6) + 2}')
        #     plt.legend()
        #     if pause:
        #         input("Press Enter to continue...") 
    
        # if plot > 2:
        #     plt.figure()
        #     plt.semilogx(Eaxis, ESH[:, 0], label=_H + ' Source Energy Distribution: ESH(E) = E fSH(E)')
        #     plt.xlabel('E (eV)')
        #     plt.ylabel('E fSH(E)')
        #     plt.title(_H + ' Source Energy Distribution: ESH(E) = E fSH(E)')
        #     for k in range(nx):
        #         plt.semilogx(Eaxis, ESH[:, k], label=f'ESH[:, {k}]', color=f'C{(k % 6) + 2}')
        #     plt.legend()
        #     if pause:
        #         input("Press Enter to continue...")
    
        # if plot > 1:
        #     plt.figure()
        #     plt.semilogx(Eaxis, ESH[:, 0], label=_H + ' Source Energy Distribution: ESH(E) = E fSH(E)')
        #     plt.xlabel('E (eV)')
        #     plt.ylabel('E fSH(E)')
        #     plt.title(_H + ' Source Energy Distribution: ESH(E) = E fSH(E)')

        #     for k in range(nx):
        #         plt.semilogx(Eaxis, ESH[:, k], label=f'ESH[:, {k}]', color=f'C{(k % 6) + 2}')
        #     plt.legend()
        #     if pause:
        #         input("Press Enter to continue...")
    
        # if plot > 0:
        #     data = [SH, SP, SHP, SH2, NuLoss * nHP, NuDis * nHP]
        #     jp = np.where(np.array(data) > 0)
        #     yrange = [min(np.array(data)[jp]), max(np.array(data)[jp])]

        #     plt.figure()
        #     plt.plot(x, SH, label=_H + ' source', color='C2')
        #     plt.plot(x, SHP, label=_Hp + ' source', color='C3')
        #     plt.plot(x, SP, label=_p + ' source', color='C4')
        #     plt.plot(x, NuLoss * nHP, label=_Hp + ' loss', color='C5')
        #     plt.plot(x, NuDis * nHP, label=_Hp + ' dissoc.', color='C6')
        #     plt.plot(x, SH2, label=_HH + ' source', color='C1')

        #     plt.xlabel('x (meters)')
        #     plt.ylabel('m⁻³ s⁻¹')
        #     plt.yscale('log')
        #     plt.ylim(yrange)
        #     plt.title('Source and Sink Profiles')
        #     plt.legend()

        #     plt.text(x[mid1], 1.2 * SH[mid1], _H + ' source', color='C2')
        #     plt.text(x[mid2], 1.2 * SHP[mid2], _Hp + ' source', color='C3')
        #     plt.text(x[mid3], 0.6 * SP[mid3], _p + ' source', color='C4')
        #     plt.text(x[mid4], 0.6 * NuLoss[mid4] * nHP[mid4], _Hp + ' loss', color='C5')
        #     plt.text(x[mid5], 0.6 * NuDis[mid5] * nHP[mid5], _Hp + ' dissoc.', color='C6')
        #     plt.text(x[mid6], 0.8 * SH2[mid6], _HH + ' source', color='C1')

        #     if pause:
        #         input("Press Enter to continue...")

        #     plt.show()
    
        # if plot > 0:
        #     gammaxH2_plus = np.zeros(nx)
        #     gammaxH2_minus = np.zeros(nx)

        #     for k in range(nx):
        #         gammaxH2_plus[k] = Vth * np.sum(Vr2pidVr * (fH2[:, i_p, k] * (vx[i_p] * dVx[i_p])))
        #         gammaxH2_minus[k] = Vth * np.sum(Vr2pidVr * (fH2[:, i_n, k] * (vx[i_n] * dVx[i_n])))

        #     gammaxH2 = gammaxH2_plus - gammaxH2_minus

        #     data = [gammaxH2_plus, gammaxH2_minus, gammaxH2]
        #     jp = np.where(np.array(data) < 1.0e32)
        #     yrange = [np.min(np.array(data)[jp]), np.max(np.array(data)[jp])]

        #     plt.figure()
        #     plt.plot(x, gammaxH2, label=_HH + ' Fluxes', color='C2')
        #     plt.xlabel('x (meters)')
        #     plt.ylabel('m⁻² s⁻¹')
        #     plt.title(_HH + ' Fluxes')
        #     plt.yscale('log')
        #     plt.ylim(yrange)
    
        Source_Error = np.zeros(nx)

        # Compute Source Error
        if Compute_Errors:
            if debrief > 1:
                print(prompt, 'Computing Source Error')
            # Test Mass Balance
            # The relationship, 2 dGammaxH2/dx - 2 SH2 + SH + SP + 2 nHp x Nuloss = 0, should be satisfied.
            dGammaxH2dx = np.zeros((nx-1))
            SH_p = np.zeros(nx-1)
            for k in range(0, nx-1):
                dGammaxH2dx[k] = (GammaxH2[k+1] - GammaxH2[k]) / (x[k+1] - x[k])
            for k in range(0, nx-1):
                SH_p[k] = 0.5*(SH[k+1] + SP[k+1] + 2*NuLoss[k+1]*nHP[k+1] - 2*SH2[k+1] + SH[k] + SP[k] + 2*NuLoss[k]*nHP[k] - 2*SH2[k])
            max_source = np.max(np.array([SH, 2*SH2]))
            for k in range(0, nx - 1):
                Source_Error[k] = np.abs(2*dGammaxH2dx[k] + SH_p[k]) / np.max(np.abs(np.array([2*dGammaxH2dx[k], SH_p[k], max_source])))
            if debrief > 0:
                print(prompt, 'Maximum Normalized Source_error =', np.max(Source_Error))
        
        # Save input parameters in common block
        vx_s = vx
        vr_s = vr
        x_s = x
        Tnorm_s = Tnorm
        mu_s = mu
        Ti_s = Ti
        vxi_s = vxi
        Te_s = Te
        n_s = n
        vxi_s = vxi
        fH2BC_s = fH2BC
        GammaxH2BC_s = GammaxH2BC
        NuLoss_s = NuLoss
        PipeDia_s = PipeDia
        fH_s = fH
        SH2_s = SH2
        fH2_s = fH2
        nHP_s = nHP
        THP_s = THP
        Simple_CX_s = Simple_CX
        Sawada_s = Sawada
        H2_H2_EL_s = H2_H2_EL
        H2_P_EL_s = H2_P_EL
        H2_H_EL_s = H2_H_EL
        H2_HP_CX_s = H2_HP_CX
        ni_correct_s = ni_correct

    # Set output parameters to single precision

    #   fH2 = float(fH2)
    #   nH2 = float(nH2)
    #   GammaxH2=float(GammaxH2)
    #   VxH2=float(VxH2)
    #   pH2=float(pH2)
    #   TH2=float(TH2)
    #   qxH2=float(qxH2)
    #   qxH2_total=float(qxH2_total)
    #   Sloss=float(Sloss)
    #   QH2=float(QH2)
    #   RxH2=float(RxH2)
    #   QH2_total=float(QH2_total)
    #   AlbedoH2=float(AlbedoH2)    
    #   nHP=float(nHP)
    #   THP=float(THP)
    #   fSH=float(fSH)
    #   SH=float(SH)
    #   SP=float(SP)
    #   SHP=float(SHP)
    #   NuE=float(NuE)
    #   NuDis=float(NuDis)
    #   piH2_xx=float(piH2_xx)
    #   piH2_yy=float(piH2_yy)
    #   piH2_zz=float(piH2_zz)
    #   RxH_H2=float(RxH_H2)
    #   RxHp_H2=float(RxHp_H2)
    #   RxP_H2=float(RxP_H2)
    #   RxP_H2_CX=float(RxP_H2_CX)
    #   RxH_H2=float(RxH_H2)
    #   EHp_H2=float(EHp_H2)
    #   EP_H2=float(EP_H2)
    #   EP_H2_CX=float(EP_H2_CX)
    #   Epara_PerpH2_H2=float(Epara_PerpH2_H2)
    #   ESH=float(ESH)
    #   Eaxis=float(Eaxis)

    #   vr=float(vr)
    #   vx=float(vx)
    #   x=float(x)   
        
    # Set common blocks 
        # Kinetic_H2_Output common block
    KH2.Output.piH2_xx = piH2_xx
    KH2.Output.piH2_yy = piH2_yy
    KH2.Output.piH2_zz = piH2_zz
    KH2.Output.RxH2CX = RxH2CX
    KH2.Output.RxH_H2 = RxH_H2
    KH2.Output.RxP_H2 = RxP_H2
    KH2.Output.RxW_H2 = RxW_H2
    KH2.Output.EH2CX = EH2CX
    KH2.Output.EH_H2 = EH_H2
    KH2.Output.EP_H2 = EP_H2
    KH2.Output.EW_H2 = EW_H2
    KH2.Output.Epara_PerpH2_H2 = Epara_PerpH2_H2

    # Kinetic_H2_Errors common block
    KH2.Errors.Max_dx = Max_dx
    KH2.Errors.vbar_error = vbar_error
    KH2.Errors.mesh_error = mesh_error
    KH2.Errors.C_Error = C_error
    KH2.Errors.CX_Error = CX_error
    KH2.Errors.H2_H2_error = H2_H2_error
    KH2.Errors.qxH2_total_error = qxH2_total_error
    KH2.Errors.QH2_total_error = QH2_total_error

    # Kinetic_H2_input common  
    KH2.Input.vx_s = vx_s
    KH2.Input.vr_s = vr_s
    KH2.Input.x_s = x_s
    KH2.Input.Tnorm_s = Tnorm_s
    KH2.Input.mu_s = mu_s
    KH2.Input.Ti_s = Ti_s
    KH2.Input.Te_s = Te_s
    KH2.Input.n_s = n_s
    KH2.Input.vxi_s = vxi_s
    KH2.Input.fH2BC_s = fH2BC_s
    KH2.Input.GammaxH2BC_s = GammaxH2BC_s
    KH2.Input.NuLoss_s = Nuloss 
    KH2.Input.PipeDia_s = PipeDia_s
    KH2.Input.fH_s = fH_s
    KH2.Input.SH2_s = SH2_s
    KH2.Input.fH2_s = fH2_s

    KH2.Input.nHP_s = nHP
    KH2.Input.THP_s = THP_s
    KH2.Input.Simple_CX_s = Simple_CX_s
    KH2.Input.Sawada_s = Sawada_s
    KH2.Input.H2_H2_EL_s = H2_H2_EL_s
    KH2.Input.H2_P_EL_s = H2_P_EL_s
    KH2.Input.H2_H_EL_s = H2_H_EL_s
    KH2.Input.H2_HP_CX_s = H2_HP_CX_s
    KH2.Input.ni_correct_s = ni_correct_s

    # kinetic_h2_internal common block  
    KH2.Internal.vr2vx2 = vr2vx2
    KH2.Internal.vr2vx_vxi2 = vr2vx_vxi2
    KH2.Internal.fw_hat = fw_hat
    KH2.Internal.fi_hat = fi_hat
    KH2.Internal.fHp_hat = fHp_hat
    KH2.Internal.EH2_P = EH2_P
    KH2.Internal.sigv = sigv
    KH2.Internal.Alpha_Loss = Alpha_Loss
    KH2.Internal.v_v2 = v_v2
    KH2.Internal.v_v = v_v
    KH2.Internal.vr2_vx2 = vr2_vx2
    KH2.Internal.vx_vx = vx_vx

    KH2.Internal.Vr2pidVrdVx = Vr2pidVrdVx
    KH2.Internal.SIG_CX = SIG_CX
    KH2.Internal.SIG_H2_H2 = SIG_H2_H2
    KH2.Internal.SIG_H2_H = SIG_H2_H
    KH2.Internal.SIG_H2_P = SIG_H2_P
    KH2.Internal.Alpha_CX = Alpha_CX
    KH2.Internal.Alpha_H2_H = Alpha_H2_H
    KH2.Internal.MH2_H2_sum = MH2_H2_sum
    KH2.Internal.Delta_nH2s = Delta_nH2s

    # kinetic_h2_moments common block
    KH2.Moments.nH = nH2
    KH2.Moments.VxH = VxH2
    KH2.Moments.TH = TH2
    
    if debug > 0:
        print(prompt, 'Finished')
        # Press_return 
    return (fH2, nHP, THP, nH2, GammaxH2, VxH2, pH2, TH2, qxH2, qxH2_total, Sloss, QH2, RxH2, QH2_total, AlbedoH2, \
        WallH2, fSH, SH, SP, SHP, NuE, NuDis, ESH, Eaxis, error)
