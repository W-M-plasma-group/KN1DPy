from key_default import key_default
from type_of import type_of
from reverse import reverse
from make_dvr_dvx import make_dvr_dvx
from create_shifted_maxwellian import create_shifted_maxwellian
from sigmav_ion_he_goto import sigmav_ion_he_goto
from jhalpha_coef import jhalpha_coef
from sigma_el_h_h import sigma_el_h_h
from sigma_el_h_hh import sigma_el_h_hh
from sigma_el_p_he import sigma_el_p_he
from sigmav_cx_h0 import sigmav_cx_h0
from press_return import press_return
from sign import sign
from sval import sval
from locate import locate
from x import x
from global_vars import mH, q, k_boltz, Twall

import numpy as np
import matplotlib.pyplot as plt


#Kinetic_HE.py
#Original IDL comments are included
# This subroutine is part of the "KN1D" atomic and molecular neutral transport code.

#   This subroutine solves a 1-D spatial, 2-D velocity kinetic neutral transport 
# problem for atomic hydrogen (H) or deuterium by computing successive generations of 
# charge exchange and elastic scattered neutrals. The routine handles electron-impact 
# ionization, proton-atom charge exchange, radiative recombination, and elastic
# collisions with hydrogenic ions, neutral atoms, and molecules.

#   The positive vx half of the atomic neutral distribution function is inputted at x(0) 
#(with arbitrary normalization) and the desired flux of hydrogen atoms entering the slab,
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

def kinetic_he(vx,vr,x,Tnorm,mu,Ti,Te,n,vxi,fHBC,GammaxHBC,PipeDia,fH2,fSH,nHP,THP,fH=None,
               truncate=1e-4,Compute_Errors=0,plot=0,debug=0,pause=0,debrief=0,Simple_CX=1,
	            Max_Gen=50,No_Johnson_Hinnov=0,No_Recomb=0,H_H_EL=0,H_P_EL=0,_H_H2_EL=0,
		        H_P_CX=0,ni_correct=0, g=None):
   
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
#		   fHBC	- fltarr(nvr,nvx), this is an input boundary condition
#			  specifying the shape of the neutral atom velocity distribution 
#			  function at location x(0). Normalization is arbitrary.
#		          Only values with positive vx, fHBC(*,nvx/2:*) are used
#		          by the code.
#	      GammaxHBC	- float, desired neutral atom flux density in the +Vx
#			  direction at location x(0) (m^-2 s^-1)
#			  fHBC is scaled to yield this flux density.
#	        PipeDia	- fltarr(nx), effective pipe diameter (meters)
#			  This variable allows collisions with the 'side-walls' to be simulated.
#			  If this variable is undefined, then PipeDia set set to zero. Zero values
#			  of PipeDia are ignored (i.e., treated as an infinite diameter).
#                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
#                         function. fH2 is normalized so that the molecular neutral density, nH2(k), is 
#			  defined as the velocity space integration: nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
#                         If this variable is undefined, then it is set equal to zero and
#                         no molecule-atom collisions are included.
#			  NOTE: dVx is velocity space differential for Vx axis and Vr2pidVr = Vr*!pi*dVr
#		                with dVr being velocity space differential for Vr axis.
#                   fSH	- fltarr(nvr,nvx,nx), atomic hydrogen source velocity distribution.
#                         fSH must be normalized so that the total atomic neutral
#                         source, SourceH(k), is defined as the velocity space integration:
#                             SourceH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
#			  fSH can be computed from IDL procedure Kinetic_H2.pro
#                         If this variable is undefined, then it is set equal to zero.
#                   nHP	- fltarr(nx), molecular ion density profile (m^-3)
#                         If this parameter is undefined, then it is set equal to zero.
#			  nHP can be computed from IDL procedure Kinetic_H2.pro
#                   THP	- fltarr(nx), molecular ion temperature profile (m^-3)
#                         If this parameter is undefined, then it is set equal to 3 eV at each grid point.
#			  THP can be computed from IDL procedure Kinetic_H2.pro
#
#  Input & Output:
#                    fH	- fltarr(nvr,nvx,nx), neutral atom velocity distribution
#                         function. 'Seed' values for this may be specified on input. 
#		          If this parameter is undefined on input, then a zero 'seed' value will be used. 
#			  The algorithm outputs a self-consistent fH.
#			  fH is normalized so that the neutral density, nH(k), is defined as 
#			  the velocity space integration: nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
#
#  Output:
#                    nH	- fltarr(nx), neutral atom density profile (m^-3)
#               GammaxH	- fltarr(nx), neutral atom flux profile (# m^-2 s^-1)
#                             computed from GammaxH(k)=Vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))
#                   VxH	- fltarr(nx), neutral atom velocity profile (m s^-1)
#                             computed from GammaxH/nH
#
#                       To aid in computing the some of the quantities below, the procedure internally
#                       defines the quantities:
#                       vr2vx2_ran(i,j,k)=vr(i)^2+(vx(j)-VxH(k))^2
#                                     which is the magnitude of 'random v^2' at each mesh point
#                       vr2vx2(i,j,k)=vr(i)^2+vx(j)^2
#                                     which is the magnitude of 'total v^2' at each mesh point
#                       q=1.602177D-19, mH=1.6726231D-27
#                       C(*,*,*) is the right hand side of the Boltzmann equation, evaluated
#                                using the computed neutral distribution function
#
#                    pH	- fltarr(nx), neutral atom pressure (eV m^-2) computed from:
#                         pH(k)~vth2*total(Vr2pidVr*(vr2vx2_ran(*,*,k)*fH(*,*,k))#dVx))*(mu*mH)/(3*q)
#                    TH	- fltarr(nx), neutral atom temperature profile (eV) computed from: TH=pH/nH
#                   qxH	- fltarr(nx), neutral atom random heat flux profile (watts m^-2) computed from:
#                          qxH(k)~vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH(*,*,k))#(dVx*(vx-VxH(k)))))*0.5*(mu*mH)
#             qxH_total	- fltarr(nx), total neutral atom heat flux profile (watts m^-2)
#                             This is the total heat flux transported by the neutrals:
#                         qxH_total=(0.5*nH*(mu*mH)*VxH*VxH + 2.5*pH*q)*VxH + piH_xx*VxH + qxH
#	     NetHSource	- fltarr(nx), net H0 source [H0 source - ionization sink - wall sink] (m^-3 s^-1) computed from
#				NetHSource(k)=total(Vr2pidVr*(C(*,*,k)#dVx))
#		   Sion	- fltarr(nx), H ionization rate (m^-3 s^-1) 
#                    QH	- fltarr(nx), rate of net thermal energy transfer into neutral atoms (watts m^-3) computed from
#                               QH(k)~vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C(*,*,k))#dVx))*0.5*(mu*mH)
#                   RxH	- fltarr(nx), rate of x momentum transfer to neutral atoms (=force, N m^-2).
#                               RxH(k)~Vth*total(Vr2pidVr*(C(*,*,k)#(dVx*(vx-VxH(k)))))*(mu*mH)
#              QH_total	- fltarr(nx), net rate of total energy transfer into neutral atoms
#                          = QH + RxH*VxH - 0.5*(mu*mH)*(Sloss-SourceH)*VxH*VxH (watts m^-3)
#               AlbedoH	- float, Ratio of atomic neutral particle flux with Vx < 0 divided by particle flux
#                          with Vx > 0  at x=x(0)
#                          (Note: For fSH non-zero, the flux with Vx < 0 will include
#                          contributions from molecular hydrogen sources within the 'slab'.
#                          In this case, this parameter does not return the true 'Albedo'.)
#	          WallH	- fltarr(nx), atomic neutral sink rate arising from hitting the 'side walls' (m^-3 s^-1)
#			   Unlike the molecules in Kinetic_H2, wall collisions result in the destruction of atoms.
#			   This parameter can be used to specify a resulting source of molecular
#			   neutrals in Kinetic_H2. (molecular source = 2 times WallH)
#
# KEYWORDS:
#   Output:
#	          error	- Returns error status: 0=no error, solution returned
#					        1=error, no solution returned
#
# COMMON BLOCK Kinetic_H_OUTPUT
#    Output:
#                piH_xx	- fltarr(nx), xx element of stress tensor (eV m^-2) computed from:
#                         piH_xx(k)~vth2*total(Vr2pidVr*(fH(*,*,k)#(dVx*(vx-VxH(k))^2)))*(mu*mH)/q - pH
#                piH_yy	- fltarr(nx), yy element of stress tensor (eV m^-2) computed from:
#                         piH_yy(k)~vth2*total((Vr2pidVr*Vr^2)*(fH(*,*,k)#dVx))*(mu*mH)/q - pH
#                piH_zz	- fltarr(nx), zz element of stress tensor (eV m^-2) = piH_yy
#			   Note: cylindrical system relates r^2 = y^2 + z^2. All other stress tensor elements are zero.
#
#           The following momentum and energy transfer rates are computed from charge-exchange collsions between species:
#                 RxHCX	- fltarr(nx), rate of x momentum transfer from hydrogren ions to atoms (=force/vol, N m^-3).
#                  EHCX	- fltarr(nx), rate of energy transfer from hydrogren ions to atoms (watts m^-3).
#               
#           The following momentum and energy transfer rates are computed from elastic collsions between species:
#                RxH2_H	- fltarr(nx), rate of x momentum transfer from neutral molecules to atoms (=force/vol, N m^-3).
#                 RxP_H	- fltarr(nx), rate of x momentum transfer from hydrogen ions to neutral atoms (=force/vol, N m^-3).
#                 EH2_H	- fltarr(nx), rate of energy transfer from neutral molecules to atoms (watts m^-3).
#                  EP_H	- fltarr(nx), rate of energy transfer from hydrogen ions to neutral atoms (watts m^-3).
#
#           The following momentum and energy transfer rates are computed from collisions with the 'side-walls'
#                 RxW_H	- fltarr(nx), rate of x momentum transfer from wall to neutral atoms (=force/vol, N m^-3).
#                  EW_H	- fltarr(nx), rate of energy transfer from wall to neutral atoms (watts m^-3).
#
#           The following is the rate of parallel to perpendicular energy transfer computed from elastic collisions
#         Epara_PerpH_H	- fltarr(nx), rate of parallel to perp energy transfer within atomic hydrogen species (watts m^-3).
#
#           Source/Sink info:
#               SourceH	- fltarr(nx), source rate of neutral atoms from H2 dissociation (from integral of inputted fSH) (m^-3 s^-1).
#                SRecom	- fltarr(nx), source rate of neutral atoms from recombination (m^-3 s^-1).
#
# KEYWORDS:
#   Input:
#	       truncate	- float, stop computation when the maximum 
#			  increment of neutral density normalized to 
#			  inputed neutral density is less than this 
#	    		  value in a subsequent generation. Default value is 1.0e-4
#
#             Simple_CX	- if set, then use CX source option (B): Neutrals are born
#                         in velocity with a distribution proportional to the local
#                         ion distribution function. Simple_CX=1 is default.
#
#                         if not set, then use CX source option (A): The CX source
#                         neutral atom distribution function is computed by evaluating the
#                         the CX cross section for each combination of (vr,vx,vr',vx')
#                         and convolving it with the neutral atom distribution function.
#                         This option requires more CPU time and memory.
#
#      	  	Max_gen	- integer, maximum number of collision generations to try including before giving up.
#                         Default is 50.
#
#     No_Johnson_Hinnov	- if set, then compute ionization and recombination rates
#			  directly from reaction rates published by Janev* for
#			  ground state hydrogen
#
#			      Ionization:    e + H(1s) -> p + e 
#			      Recombination: e + p -> H(1s) + hv
#
#			  *Janev, R.K., et al, "Elementary processes in hydrogen-helium plasmas",
#			   (Springer-Verlag, Berlin ; New York, 1987)
#
#			  Otherwise, compute ionization and recombination rates using
#		          results from the collisional-radiative model published by Johnson
#			  and Hinnov [L.C.Johnson and E. Hinnov, J. Quant. Spectrosc. Radiat.
# 			  Transfer. vol. 13 pp.333-358]. This is the default.
#			  Note: charge exchange is always computed using the ground state reaction
#		          rates published by Janev:
#
#			      Charge Exchange: p + H(1s) -> H(1s) + p
#			  
#	      No_Recomb	- if set, then DO NOT include recombination as a source of atomic neutrals
#		          in the algorithm
#
#	 	 H_H_EL	- if set, then include H -> H elastic self collisions
#			     Note: if H_H_EL is set, then algorithm iterates fH until
#	                     self consistent fH is achieved.
#	 	 H_P_CX	- if set, then include H -> H(+) charge exchange collisions 
#	         H_P_EL	- if set, then include H -> H(+) elastic collisions 
#	        H_H2_EL	- if set, then include H -> H2 elastic collisions 
#            ni_correct	- if set, then algorithm corrects hydrogen ion density
#			     according to quasineutrality: ni=ne-nHP. Otherwise, nHP is assumed to be small.
#
#	 Compute_Errors	- if set, then return error estimates in common block Kinetic_H_ERRORS below
#
#		   plot	- 0= no plots, 1=summary plots, 2=detail plots, 3=very detailed plots
#		  debug	- 0= do not execute debug code, 1=summary debug, 2=detail debug, 3=very detailed debug
#	        debrief	- 0= do not print, 1=print summary information, 2=print detailed information
#	          pause	- if set, then pause between plots
#
# COMMON BLOCK Kinetic_H_ERRORS
#
#	if COMPUTE_ERRORS keyword is set then the following is returned in common block Kinetic_H_ERRORS
#
#	         Max_dx	- float(nx), Max_dx(k) for k=0:nx-2 returns maximum 
#			  allowed x(k+1)-x(k) that avoids unphysical negative 
#			  contributions to fH
#	     Vbar_error	- float(nx), returns numerical error in computing
#			  the speed of ions averged over maxwellian distribution.
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
#	      H_H_error	-  fltarr(nx,[0,1,2]) return normalized errors associated with 
#		           particle [0], x-momentum [1], and total energy [2] convervation of the elastic self-collision operator
#
#       qxH_total_error	- fltarr(nx), normalized error estimate in computation of qxH_total
#        QH_total_error	- fltarr(nx), normalized error estimate in computation of QH_total

   prompt='Kinetic_H => '

   #Kinetic_H_input common block global variables (may be replaced in the future, but initially we'll use this)
   vx_s = g.Kinetic_H_input_vx_s
   vr_s = g.Kinetic_H_input_vr_s
   x_s = g.Kinetic_H_input_x_s
   Tnorm_s = g.Kinetic_H_input_Tnorm_s
   mu_s = g.Kinetic_H_input_mu_s
   Ti_s = g.Kinetic_H_input_Ti_s
   Te_s = g.Kinetic_H_input_Te_s
   n_s = g.Kinetic_H_input_n_s
   vxi_s = g.Kinetic_H_input_vxi_s
   fHBC_s = g.Kinetic_H_input_fHBC_s
   GammaxHBC_s = g.Kinetic_H_input_GammaxHBC_s
   PipeDia_s = g.Kinetic_H_input_PipeDia_s
   fH2_s = g.Kinetic_H_input_fH2_s
   fSH_s = g.Kinetic_H_input_fSH_s
   nHP_s = g.Kinetic_H_input_nHP_s
   THP_s = g.Kinetic_H_input_THP_s
   fH_s = g.Kinetic_H_input_fH_s
   Simple_CX_s = g.Kinetic_H_input_Simple_CX_s
   JH_s = g.Kinetic_H_input_JH_s
   Recomb_s = g.Kinetic_H_input_Recomb_s
   H_H_EL_s = g.Kinetic_H_input_H_H_EL_s
   H_P_EL_s = g.Kinetic_H_input_H_P_EL_s
   H_H2_EL_s = g.Kinetic_H_input_H_H2_EL_s
   H_P_CX_s = g.Kinetic_H_input_H_P_CX_s

   #Kinetic_H_internal common block
   vr2vx2 = g.Kinetic_H_internal_vr2vx2
   vr2vx_vxi2 = g.Kinetic_H_internal_vr2vx_vxi2
   fi_hat = g.Kinetic_H_internal_fi_hat
   ErelH_P = g.Kinetic_H_internal_ErelH_P
   Ti_mu = g.Kinetic_H_internal_Ti_mu
   ni = g.Kinetic_H_internal_ni
   sigv = g.Kinetic_H_internal_sigv
   alpha_ion = g.Kinetic_H_internal_alpha_ion
   v_v2 = g.Kinetic_H_internal_v_v2
   v_v = g.Kinetic_H_internal_v_v
   vr2_vx2 = g.Kinetic_H_internal_vr2_vx2
   vx_vx = g.Kinetic_H_internal_vx_vx
   Vr2pidVrdVx = g.Kinetic_H_internal_Vr2pidVrdVx
   SIG_CX = g.Kinetic_H_internal_SIG_CX
   SIG_H_H = g.Kinetic_H_internal_SIG_H_H
   SIG_H_H2 = g.Kinetic_H_internal_SIG_H_H2
   SIG_H_P = g.Kinetic_H_internal_SIG_H_P
   Alpha_CX = g.Kinetic_H_internal_Alpha_CX
   Alpha_H_H2 = g.Kinetic_H_internal_Alpha_H_H2
   Alpha_H_P = g.Kinetic_H_internal_Alpha_H_P
   MH_H_sum = g.Kinetic_H_internal_MH_H_sum
   Delta_nHs = g.Kinetic_H_internal_Delta_nHs
   Sn = g.Kinetic_H_internal_Sn
   Rec = g.Kinetic_H_internal_Rec

   #Kinetic_H_H2_Moments common block
   nH2 = g.Kinetic_H_H2_Moments_nH2
   VxH2 = g.Kinetic_H_H2_Moments_VxH2
   TH2 = g.Kinetic_H_H2_Moments_TH2

   #Internal Debug switches

   shifted_Maxwellian_debug = 0
   CI_Test = 1
   Do_Alpha_CX_Test = 0

   #Internal Tolerances

   DeltaVx_tol = .01
   Wpp_tol = .001

   #	Test input parameters (copied from kinetic_h.py because the idl code is the same in this section)

   if debug > 0:
      plot = plot > 1
      debrief = debrief > 1
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
   dx = x-np.roll(x, 1)
   dx = dx[1:] # moved to new line- fixed error

   while error == 0:
      notpos = dx[dx <= 0]
      if notpos.size > 0:
	      print(prompt + 'x[*] must be increasing with index!')
         error = 1
	   if nvx % 2 != 0:
		   print(prompt + 'Number of elements in vx must be even!') 
		   error = 1
	   if Ti.size != nx:
		   print(prompt + 'Number of elements in Ti and x do not agree!')
		   error = 1
		#	if type_of(vxi) eq 0 then vxi=dblarr(nx) - Doesn't really work in Python
	   if vxi.size != nx:
		   print(prompt + 'Number of elements in vxi and x do not agree!')
		   error = 1
	   if Te.size != nx:
		   print(prompt + 'Number of elements in Te and x do not agree!')
		   error = 1
	   if n.size != nx:
		   print('Number of elements in n and x do not agree!')
		   error = 1
	   if PipeDia.size != nx:
		   print('Number of elements in PipeDia and x do not agree!') # Fixed error message- previously copied from n.size!=nx
		   error = 1
		if fHBC[0, :].size != nvr:
			print(prompt + 'Number of elements in fHBC[0,:] and vr do not agree!')
			error = 1
		if fHBC[:, 0].size != nvx:
			print(prompt + 'Number of elements in fHBC[:,0] and vx do not agree!')
			error = 1
		if fH2[0, 0, :].size != nvr:
			print(prompt + 'Number of elements in fH2[0,0,:] and vr do not agree!')
			error = 1
		if fH2[0, :, 0].size != nvx:
			print(prompt + 'Number of elements in fH2[0,:,0] and vx do not agree!')
			error = 1
		if fH2[:, 0, 0].size != nx:
			print(prompt + 'Number of elements in fH2[:,0,0] and x do not agree!')
			error = 1
		if fSH[0, 0, :].size != nvr:
			print(prompt + 'Number of elements in fSH[0,0,:] and vr do not agree!')
			error = 1
		if fSH[0, :, 0].size != nvx:
			print(prompt + 'Number of elements in fSH[0,:,0] and vx do not agree!')
			error = 1
		if fSH[:, 0, 0].size != nx:
			print(prompt + 'Number of elements in fSH[:,0,0] and x do not agree!')
			error = 1
		if nHP.size != nx:
			print(prompt + 'Number of elements in nHP and x do not agree!')
			error = 1
		if THP.size != nx:
			print(prompt + 'Number of elements in nHP and x do not agree!')
			error = 1
		if fH == None:
			fH = np.zeros((nx, nvx, nvr)) # Set fH if not left as None in function call
		if fH[0, 0, :].size != nvr:
			print(prompt + 'Number of elements in fH[0,0,:] and vr do not agree!')
			error = 1
		if fH[0, :, 0].size != nvx:
			print(prompt + 'Number of elements in fH[0,:,0] and vx do not agree!')
			error = 1
		if fH[:, 0, 0].size != nx:
			print(prompt + 'Number of elements in fH[:,0,0] and x do not agree!')
			error = 1
		if np.sum(abs(vr)) == 0:
			print(prompt + 'vr is all 0!')
			error = 1
		ii = vr[vr <= 0]
		if ii.size > 0:
			print(prompt + 'vr contains zero or negative element(s)!')
			error = 1
		if np.sum(abs(vx)) == 0:
			print(prompt + 'vx is all 0!')
			error = 1
		if np.sum(x) <= 0:
			print(prompt + 'Total(x) is less than or equal to 0!')
			error = 1
		if mu not in [1, 2]:
			print(prompt + 'mu must be 1 or 2!')
			error = 1
		break
   if error == 1:
		if debug > 0:
			print(prompt + 'Finished')
		return
	_e = 'e!U-!N'
	_hv = 'hv'
	if mu == 1:
		_p = 'H!U+!N'
		_H = 'H!U0!N'
		_H1s = 'H(1s)'
		_Hs = 'H!U*!N(2s)'
		_Hp = 'H!U*!N(2p)'
		_Hn2 = 'H!U*!N(n=2)'
		_Hn3 = 'H!U*!N(n=3)'
		_Hn =' H!U*!N(n>=2)'
		_HH = 'H!D2!N'
		_Hp = 'H!D2!U+!N'
	else:
		_p = 'D!U+!N'
		_H = 'D!U0!N'
		_H1s = 'D(1s)'
		_Hs = 'D!U*!N(2s)'
		_Hp = 'D!U*!N(2p)'
		_Hn2 = 'D!U*!N(n=2)'
		_Hn3 = 'D!U*!N(n=3)'
		_Hn = 'D!U*!N(n>=2)'
		_HH = 'D!D2!N'
		_Hp = 'D!D2!U+!N'
	plus = ' + '
	arrow = ' -> '
	elastic = ' (elastic)'
	_R1 = _e + plus + _H1s + arrow + _e + plus + _p + plus + _e
	_R2 = _e + plus + _p + arrow + plus + _H1s + plus + _hv
	_R3 = _p + plus + _H1s + arrow + _H1s + plus + _p
	_R4 = _H + plus + _p + arrow + _H + plus + _p + elastic
	_R5 = _H + plus + _HH + arrow + _H + plus + _HH + elastic
	_R6 = _H + plus + _H + arrow + _H + plus + _H + elastic
	_Rn = [' ', _R1, _R2, _R3, _R4, _R5, _R6]
	while error == 0:
		#i_n=vx[vx<0] was original version; now rewritten
		i_n, count=np.nonzero(vx < 0)[0],np.count_nonzero(vx < 0)
		if count < 1:
			print(prompt + 'vx contains no negative elements!')
			error = 1
			break
		#i_p=vx[vx>0] was original version; now rewritten
		i_p, count=np.nonzero(vx > 0)[0],np.count_nonzero(vx > 0)
		if count < 1:
			print(prompt + 'vx contains no positive elements!')
			error = 1
			break
		#i_z=vx[vx==0] was original version; now rewritten
		i_z, count = np.nonzero(vx == 0)[0],np.count_nonzero(vx == 0)
		if count > 0:
			print(prompt + 'vx contains one or more zero elements!')
			error = 1
			break
		#	rewritten
		if i_p.size == i_n.size:
			for i in range(i_n.size):
				if vx[i_n[i]] != -vx[i_p[i_p.size-i-1]]:
					print(prompt + 'vx array elements are not symmetric about zero!')
					error = 1
					break
		fHBC_input = np.zeros(fHBC.shape)
		fHBC_input[i_p, 0] = fHBC[i_p, 0]
		test = np.sum(fHBC_input)
		if test <= 0.0 and abs(GammaxHBC) > 0:
			print(prompt + 'Values for fHBC(*,*) with vx > 0 are all zero!')
			error = 1
			break
		break
	if error == 1:
		if debug > 0:
			print(prompt+'Finished')
		return

	#	Output variables (again, same as kinetic_h.py)

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

	#	Internal variables (first 7 variables are the same as kinetic_h.py)

   Work = np.zeros((nvx*nvr))
   fHG = np.zeros((nx, nvx, nvr))
   NHG = np.zeros((Max_Gen+1, nx))
   if ihe == 1:
      mu_temp = 4
      Vth = sqrt(2*q*Tnorm/(mu_temp*mH))
   else:
      Vth = sqrt(2*q*Tnorm/(mu*mH))
   Vth2 = vth*vth
   Vth3 = Vth2*Vth
   fHs = np.zeros(nx)
   nHs = np.zeros(nx)
   Alpha_H_H = np.zeros((nvr, nvx))
   Omega_H_P = np.zeros(nx)
   Omega_H_H2 = np.zeros(nx)
   Omega_H_H = np.zeros(nx)
   VxHG = np.zeros(nx)
   THG = np.zeros(nx)
	Wperp_paraH = np.zeros(nx)
	vr2vx2_ran2 = np.zeros((nvr. nvx))
	vr2_2vx_ran2 = np.zeros((nvr, nvx))
	vr2_2vx2_2D = np.zeros((nvr, nvx))
	RxCI_CX = np.zeros(nx)
	RxCI_H2_H = np.zeros(nx)
	RxCI_P_H = np.zeros(nx)
	Epara_Perp_CI = np.zeros(nx)
	CI_CX_error = np.zeros(nx)
	CI_H2_H_error = np.zeros(nx)
	CI_P_H_error = np.zeros(nx)
	CI_H_H_error = np.zeros(nx)
	Maxwell = np.zeros((nvr, nvx, nx))

	#make_dvr_dvx is not a translated function, so it is commented out currently
	#Make_dVr_dVx(vr,vx,Vr2pidVr,VrVr4pidVr,dVx,vrL=vrL,vrR=vrR,vxL=vxL,vxR=vxR,
                #Vol=Vol,Vth_DeltaVx=Vth_DVx,Vx_DeltaVx=Vx_DVx,Vr_DeltaVr=Vr_DVr,Vr2Vx2=Vr2Vx2_2D,
                #jpa=jpa,jpb=jpb,jna=jna,jnb=jnb)

	#vr^2-2*vx^2

   for i=0 in range(nvr-1) do vr2_2vx2_2D(i,*)=vr(i)^2-2*vx^2

	#Theta-prime coordinate

   ntheta=5     &; use 5 theta mesh points for theta integration
   dTheta=replicate(1.0d0,ntheta)/ntheta
   theta=!pi*(dindgen(ntheta)/ntheta+0.5/ntheta)
   cos_theta=cos(theta)

	#sgbaek I want to check the temperature of the input neutral
	#distribution function
    
	#vx_shift=dblarr(nx)
	#vx_shift=vxi
	#Tmaxwell=Ti
	#Tmaxwell[*]=0.025
	#mol=1
	#mu=4
	#shifted_Maxwellian_debug=1
	#@create_shifted_maxwellian.include 
	#mu=2




	#Scale input  distribution function to agree with desired flux


   gamma_input = 1.0
   if abs(GammaxHBC) > 0:
		gamma_input = vth * total(Vr2pidVr*(fHBC_input(vx*dVx)))
	ratio = abs(GammaxHBC)/gamma_input        
	fHBC_input = fHBC_input*ratio
	if abs(ratio - 1) > 0.01 * truncate:
    	fHBC = fHBC_input

   fH(None, ip, 0) = fHBC_input(None, ip)

	#if fH2 is zero, then turn off elastic H2 <-> H collisions

   H_H2_EL=_H_H2_EL
   if total(fH2) <= 0.0:
      H_H2_EL=0

	#Set iteration scheme

   fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations)=0
   if (H_H_EL != 0) or (H_P_EL != 0) or (H_H2_EL != 0):
		fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations) = 1

   	fH_generations=0
   if (fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations) != 0) or (H_P_CX != 0):
		fH_generations = 1

	#Set flags to make use of previously computed local parameters 

   New_Grid=1
   if type_of(vx_s) != 0:
		test = 0
		ii = np.where(vx_s != vx)[0]
		test += len(ii)
		ii = np.where(vr_s != vr)[0]
		test += len(ii)
		ii = np.where(x_s != x)[0]
    	test += len(ii)
    	ii = np.where(Tnorm_s != Tnorm)[0]
    	test += len(ii)
    	ii = np.where(mu_s != mu)[0]
    	test += len(ii)
    	if test <= 0:
        	New_Grid = 0

   
	New_Protons = 1
	if type(Ti_s) != 0:
		test = 0
		ii = np.where(Ti_s != Ti)[0]
		test += len(ii)
		ii = np.where(n_s != n)[0]
		test += len(ii)
		ii = np.where(vxi_s != vxi)[0]
		test += len(ii)
		if test <= 0:
			New_Protons = 0

	New_Molecular_Ions = 1
	if type(nHP_s) != 0:
		test = 0
		ii = np.where(nHP_s != nHP)[0]
		test += len(ii)
		ii = np.where(THP_s != THP)[0]
		test += len(ii)
		if test <= 0:
			New_Molecular_Ions = 0

	New_Electrons = 1
	if type(Te_s) != 0:
		test = 0
		ii = np.where(Te_s != Te)[0]
		test += len(ii)
		ii = np.where(n_s != n)[0]
		test += len(ii)
		if test <= 0:
			New_Electrons = 0

	New_fH2 = 1
	if type(fH2_s) != 0:
		ii = np.where(fH2_s != fH2)[0]
		if len(ii) <= 0:
			New_fH2 = 0

	New_fSH = 1
	if type(fSH_s) != 0:
		ii = np.where(fSH_s != fSH)[0]
		if len(ii) <= 0:
			New_fSH = 0

	New_Simple_CX = 1
	if type(Simple_CX_s) != 0:
		ii = np.where(Simple_CX_s != Simple_CX)[0]
		if len(ii) <= 0:
			New_Simple_CX = 0

	New_H_Seed = 1
	if type(fH_s) != 0:
		ii = np.where(fH_s != fH)[0]
		if len(ii) <= 0:
			New_H_Seed = 0

   Do_sigv = New_Grid or New_Electrons
   Do_ni = New_Grid or New_Electrons or New_Protons or New_Molecular_Ions
   Do_fH2_moments = (New_Grid or New_fH2) and total(fH2) > 0.0
   Do_Alpha_CX = (New_Grid or (type(Alpha_CX) == 0) or Do_ni or New_Simple_CX) and H_P_CX
   Do_SIG_CX = (New_Grid or (type(SIG_CX) == 0) or New_Simple_CX) and (Simple_CX == 0) and Do_Alpha_CX
   Do_Alpha_H_H2 = (New_Grid or (type(Alpha_H_H2) == 0) or New_fH2) and H_H2_EL
   Do_SIG_H_H2 = (New_Grid or (type(SIG_H_H2) == 0)) and Do_Alpha_H_H2
   Do_SIG_H_H = (New_Grid or (type(SIG_H_H) == 0)) and H_H_EL
   Do_Alpha_H_P = (New_Grid or (type(Alpha_H_P) == 0) or Do_ni) and H_P_EL
   Do_SIG_H_P = (New_Grid or (type(SIG_H_P) == 0)) and Do_Alpha_H_P
   Do_v_v2 = (New_Grid or (type(v_v2) == 0)) and (CI_Test or Do_SIG_CX or Do_SIG_H_H2 or Do_SIG_H_H or Do_SIG_H_P)

   print('do_sigv         =', Do_sigv)
   print('do_ni           =', Do_ni)
   print('do_fh2_moments  =', Do_fH2_moments)
   print('do_alpha_cx     =', Do_Alpha_CX)
   print('do_sig_cx       =', Do_SIG_CX)
   print('do_alpha_h_h2   =', Do_Alpha_H_H2)
   print('do_sig_h_h2     =', Do_SIG_H_H2)
   print('do_sig_h_h      =', Do_SIG_H_H)
   print('do_alpha_h_p    =', Do_Alpha_H_P)
   print('do_sig_h_p      =', Do_SIG_H_P)
   print('do_v_v2         =', Do_v_v2)

   nH2 = np.zeros(nx)
   vxH2 = np.zeros(nx)
   TH2 = np.zeros(nx) + 1.0
   if Do_fH2_moments:
      if debrief > 1:
         print(prompt + 'Computing vx and T moments of fH2')


   # Compute x flow velocity and temperature of molecular species

   for k in range(nx):
      nH2[k] = np.sum(Vr2pidVr * (fH2[:, :, k] * dVx))
      if nH2[k] > 0:
        vxH2[k] = vth * np.sum(Vr2pidVr * (fH2[:, :, k] * (vx * dVx))) / nH2[k]
        for i in range(nvr):
            vr2vx2_ran2[i, :] = vr[i] ** 2 + (vx - vxH2[k] / vth) ** 2
        TH2[k] = (2 * mu * mH) * vth2 * np.sum(Vr2pidVr * ((vr2vx2_ran2 * fH2[:, :, k]) * dVx)) / (3 * q * nH2[k])

   if New_Grid:
      if debrief > 1:
         print(prompt + 'Computing vr2vx2, vr2vx_vxi2, ErelH_P')


   # Magnitude of total normalized v^2 at each mesh point

   vr2vx2 = np.zeros((nvr,nvx,nx))
   for i in range(nvr):
      for k in range(nx):
         vr2vx2[i, :, k] = vr[i] ** 2 + vx ** 2


   # Magnitude of total normalized (v-vxi)^2 at each mesh point

   vr2vx_vxi2 = np.zeros((nvr,nvx,nx))
   for i in range(nvr):
      for k in range(nx):
         vr2vx_vxi2[i, :, k] = vr[i]  ** 2 + (vx-vxi[k]/vth) ** 2

   # Atomic hydrogen ion energy in local rest frame of plasma at each mesh point

   ErelH_P = 0.5 * mH * vr2vx_vxi2 * vth2 / q
   ErelH_P = np.logical_and(ErelH_P > 0.1, ErelH_P < 2.0E4) # sigmav_cx does not handle neutral energies below 0.1 eV or above 20 keV

   if New_Protons:
      if debrief > 1:
         print(prompt + 'Computing Ti/mu at each mesh point')


   # Ti/mu at each mesh point

   Ti_mu = np.zeros((nvr,nvx,nx))
   for k in range(nx):
      Ti_mu[:, :, k] = Ti[k]/mu

   # Compute Fi_hat

   # sgbaek - this is for plasma ions (so mu=2 for deuterium)
   if debrief > 1:
      print(prompt + 'Computing fi_Hat')
   vx_shift = vxi
   Tmaxwell = Ti
   mol = 1
   # @create_shifted_maxwellian.include (not sure how to deal with this)
   fi_hat = Maxwell

   if Compute_errors:
      if debrief > 1:
         print(prompt + 'Computing Vbar_Error')


   # Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti(x)/Tnorm)/sqrt(!pi)

   vx_shift = np.zeros(nx)
   Tmaxwell = Ti
   mol = 1
   Shifted_Maxwellian_Debug = 0

   if ihe == 1:
      mu = 4
      # @create_shifted_maxwellian.include (figure this out later)
      mu = 2
   else:
      # @create_shifted_maxwellian.include
      pass

   vbar_test = np.zeros((nvr, nvx, ntheta))
   for m in range(ntheta):
      vbar_test[:, :, m] = vr2vx2[:, :, 0]

   _vbar_test = vth * np.sqrt(vbar_test.flatten())  # sgbaek - careful! vth is for helium
   vbar_test = np.zeros((nvr, nvx))
   vbar_test[:] = _vbar_test[:, np.newaxis] * dtheta

   Vbar_Error = np.zeros(nx)
   for k in range(nx):
      vbar = np.sum(Vr2pidVr * (vbar_test * Maxwell[:, :, k]) * dVx)
      vbar_exact = 2 * Vth * np.sqrt(Ti[k] / Tnorm) / np.sqrt(np.pi)  # sgbaek, vth is for helium
      Vbar_Error[k] = np.abs(vbar - vbar_exact) / vbar_exact

   if debrief > 0:
      print(prompt + 'Maximum Vbar error =', np.max(Vbar_Error))



   if Do_ni:
      if debrief > 1:
         print(prompt + 'Computing ni profile')
      ni = n
      if ni_correct:
         ni = n - nHP
      ni = np.where(ni > 0.01 * n, ni, 0)

   if Do_sigv:
      if debrief > 1:
         print(prompt + 'Computing sigv')


      # Compute sigmav rates for each reaction with option to use rates
      # from CR model of Johnson-Hinnov

      sigv = np.zeros((nx,3))
      #________________________________________________________________________________
      # Reaction R1:  e + He -> e + He(+) + e 
      #________________________________________________________________________________
      if JH:
         sigv[:, 0] = sigmav_ion_he_goto(n, Te, no_null=True)
         print('sigmav_ion_he_goto(n,Te)')
      else:
         sigv[:, 0] = sigmav_ion_he(Te)  # sigv[:, 0] = sigmav_ion_h0(Te)
         print('sigmav_ion_he(Te)')

      #________________________________________________________________________________
      # Reaction R2:  e + H(+) -> H(1s) + hv
      #________________________________________________________________________________
      if JH:
         sigv[:, 2] = JHAlpha_Coef(n, Te, no_null=True)
      else:
         sigv[:, 2] = SigmaV_rec_H1s(Te)

      # H ionization rate (normalized by vth) = reaction 1

      alpha_ion = n*sigv[:, 1]/vth

      # Recombination rate (normalized by vth) = reaction 2

      Rec = n*sigv[:, 2]/vth

   # Compute Total Atomic Hydrogen Source

   Sn = np.zeros((nvr,nvx,nx))

   # Add Recombination (optionally) and User-Supplied Hydrogen Source (velocity space distribution)

   for k in range(nx):
      Sn[:, ;, k] = fSH[:, :, k]/vth
      if Recomb:
         Sn[:, :, k] = Sn[:, :, k]+fi_hat[:, :, k]*ni(k)*Rec[k]
   print('min and max (Sn)= ', np.min(Sn), np.max(Sn))

   #________________________________________________________________________________
   # Set up arrays for charge exchange and elastic collision computations, if needed
   #________________________________________________________________________________

   if Do_v_v2 == 1:
      if debrief > 1:
         print(prompt + 'Computing v_v2, v_v, vr2_vx2, and vx_vx')

      # v_v2=(v-v_prime)^2 at each double velocity space mesh point, including theta angle

      v_v2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))

      # vr2_vx2=0.125* [ vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2 ]
      #         at each double velocity space mesh point, including theta angle

      vr2_vx2 = np.zeros((nvr,nvx,nvr,nvx,ntheta))
      for m in range(ntheta):
         for l in range(nvx):
            for k in range(nvr):
               for i in range(nvr):
                  v_v2[i, :, k, l, m] = vr[i]**2+vr[k]**2-2*vr[i]*vr[k]*np.cos_theta(m)+(vx[:]-vx[l])**2
                  vr2_vx2[i, :, k, l, m] = vr[i]**2+vr[k]**2-2*vr[i]*vr[k]*np.cos_theta(m)-2*(vx[:]-vx[l])**2

      # v_v=|v-v_prime| at each double velocity space mesh point, including theta angle

      v_v = np.sqrt(v_v2)

      # vx_vx=(vx-vx_prime) at each double velocity space mesh point

      vx_vx = np.zeros((nvr,nvx,nvr,nvx))
      for j in range(nvx):
         for l in range(nvx):
            vx_vx[:, j, :, l] = vx[j]-vx[l]

      # Set vr'2pidVr'*dVx' for each double velocity space mesh point

      Vr2pidVrdVx = np.zeros((nvr,nvx,nvr,nvx))
      for k in range(nvr):
         Vr2pidVrdVx[:, :, k, :] = Vr2pidVr[k]
      for l in range(nvx):
         Vr2pidVrdVx[:, :, :, l] = Vr2pidVrdVx[:, :, :, l]*dVx[l]

   if Simple_CX == 0 and Do_SIG_CX == 1:
      if debrief > 1:
         print(prompt+'Computing SIG_CX')
      #________________________________________________________________________________
      # Option (A) was selected: Compute SigmaV_CX from sigma directly.
      # In preparation, compute SIG_CX for present velocity space grid, if it has not 
      # already been computed with the present input parameters
      #________________________________________________________________________________

      # Compute sigma_cx * v_v at all possible relative velocities

      _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
      _Sig[:] = v_v*sigma_cx_h0(v_v2*(0.5*mH*vth2/q))

      # Set SIG_CX = vr' x Integral{v_v*sigma_cx} over theta=0,2pi times differential velocity space element vr'2pidVr'*dVx'

      SIG_CX = np.zeros((nvr*nvx,nvr*nvx))
      SIG_CX[:] = Vr2pidVrdVx*np.concatenate((_Sig, dtheta))

      # SIG_CX is now vr' * sigma_cx(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])


   if Do_SIG_H_H == 1:
      if debrief > 1:
         print(prompt+'Computing SIG_H_H')
      #________________________________________________________________________________
      # Compute SIG_H_H for present velocity space grid, if it is needed and has not 
      # already been computed with the present input parameters
      #________________________________________________________________________________

      # Compute sigma_H_H * vr2_vx2 * v_v at all possible relative velocities

      _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
      _Sig[:] = vr2_vx2*v_v*sigma_EL_H_H(v_v2*(0.5*mH*mu*vth2/q),/VIS)/8.0

      # Note: For viscosity, the cross section for D -> D is the same function of
      #       center of mass energy as H -> H.

      # Set SIG_H_H = vr' x Integral{vr2_vx2*v_v*sigma_H_H} over theta=0,2pi times differential velocity space element vr'2pidVr'*dVx'

      SIG_H_H = np.zeros((nvr*nvx,nvr*nvx))
      SIG_H_H[:] = Vr2pidVrdVx*np.concatenate((_Sig, dtheta))

      # SIG_H_H is now vr' * sigma_H_H(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

   if Do_SIG_H_H2 == 1:
      if debrief > 1:
         print(prompt+'Computing SIG_H_H2')
      #________________________________________________________________________________
      # Compute SIG_H_H2 for present velocity space grid, if it is needed and has not 
      # already been computed with the present input parameters
      #________________________________________________________________________________

      # Compute sigma_H_H2 * v_v at all possible relative velocities

      _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
      _Sig[:] = v_v*sigma_EL_H_HH(v_v2*(0.5*mH*vth2/q))

      # NOTE: using H energy here for cross-sections tabulated as H->H2

      # Set SIG_H_H2 = vr' x vx_vx x Integral{v_v*sigma_H_H2} over theta=0,2pi times differential velocity space element vr'2pidVr'*dVx'

      SIG_H_H2 = np.zeros((nvr*nvx,nvr*nvx))
      SIG_H_H2[:] = Vr2pidVrdVx*vx_vx*np.concatenate((_Sig, dtheta))

      # SIG_H_H2 is now vr' *vx_vx * sigma_H_H2(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

   print('sgbaek do sig_h_p' + do_sig_h_p)
   if Do_SIG_H_P == 1:
      if debrief > 1:
         print(prompt+'Computing SIG_H_P')
      #________________________________________________________________________________
      # Compute SIG_H_P for present velocity space grid, if it is needed and has not 
      # already been computed with the present input parameters
      #________________________________________________________________________________

      # Compute sigma_H_P * v_v at all possible relative velocities

      _Sig = np.zeros((nvr*nvx*nvr*nvx,ntheta))
      print('Replace sigma_H_P with sigma_HE_P0000000000000000000000000000000000')
      mHelium=4*mH
      dum=(0.5*mHelium*vth2/q)/2 # to scale the cross section: sigma_D (E)= sigma_H (E/2); kn1d manual, page 58
      _Sig[:] = v_v*sigma_EL_P_HE(v_v2*(dum)) #sgbaek v_v*sigma_EL_P_H(v_v2*(0.5*mH*vth2/q))
      #_Sig[:]=v_v*sigma_EL_P_H(v_v2*(0.5*mH*vth2/q))


      # Set SIG_H_P = vr' x vx_vx x Integral{v_v*sigma_H_P} over theta=0,2pi times differential velocity space element vr'2pidVr'*dVx'

      SIG_H_P = np.zeros((nvr*nvx,nvr*nvx))
      SIG_H_P[:] = Vr2pidVrdVx*vx_vx*np.concatenate((_Sig, dtheta))

      # SIG_H_P is now vr' *vx_vx * sigma_H_P(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])

   #________________________________________________________________________________ 
   # Compute Alpha_CX for present Ti and ni, if it is needed and has not
   # already been computed with the present parameters
   #________________________________________________________________________________ 
   if Do_Alpha_CX == 1:
      if debrief > 1:
         print(prompt+'Computing Alpha_CX')

      if Simple_CX:
         #________________________________________________________________________________
         # Option (B): Use maxwellian weighted <sigma v>
         #________________________________________________________________________________

         # Charge Exchange sink rate

         alpha_cx=sigmav_cx_H0[Ti_mu,ErelH_P]/vth
         for k in range(nx):
            alpha_cx[:, :, k] = alpha_cx[:, :, k]*ni[k]
         #________________________________________________________________________________

      else:
         #________________________________________________________________________________
         # Option (A): Compute SigmaV_CX from sigma directly via SIG_CX
         #________________________________________________________________________________

         alpha_cx = np.zeros((nvr,nvx,nx))
         for k in range(nx):
            Work[:] = fi_hat[:, :, k]*ni[k]
            alpha_cx[:, :, k] = np.concatenate((SIG_CX, Work))
         if do_alpha_cx_test:
            alpha_cx_test=sigmav_cx_H0[Ti_mu,ErelH_P]/vth
            for k in range(nx):
               alpha_cx_test[:, :, k] = alpha_cx_test[:, :, k]*ni[k]
            print('Compare alpha_cx and alpha_cx_test')
            press_return()
   #________________________________________________________________________________ 
   # Compute Alpha_H_H2 for inputted fH, if it is needed and has not
   # already been computed with the present input parameters
   #________________________________________________________________________________ 
   if Do_Alpha_H_H2 == 1:
      if debrief > 1:
         print(prompt+'Computing Alpha_H_H2')
      Alpha_H_H2 = np.zeros((nvr,nvx,nx))
      for k in range(nx):
         Work[:] = fH2[:, :, k]
         Alpha_H_H2[:, :, k] = np.concatenate((SIG_H_H2, Work))
   #________________________________________________________________________________ 
   # Compute Alpha_H_P for present Ti and ni 
   # if it is needed and has not already been computed with the present parameters
   #________________________________________________________________________________ 
   if Do_Alpha_H_P == 1:
      if debrief > 1:
         print(prompt+'Computing Alpha_H_P')
      Alpha_H_P = np.zeros((nvr,nvx,nx))
      for k in range(nx):
         Work[:] = fi_hat[:, :, k]*ni[k]
         Alpha_H_P[:, :, k] = np.concatenate((SIG_H_P, Work))

   #________________________________________________________________________________
   # Compute nH
   #________________________________________________________________________________
   for k in range(nx):
      nH[k] = np.sum(Vr2pidVr*np.concatenate((fH[:, :, k], dVx)))

   if New_H_Seed:
      MH_H_sum = np.zeros((nvr,nvx,nx))
      Delta_nHs = 1.0

   # Compute Side-Wall collision rate

   gamma_wall = np.zeros((nvr,nvx,nx))
   for k in range(nx):
      if PipeDia(k) > 0.0:
         for j in range(nvx):
            gamma_wall[:, j, k] = 2*vr/PipeDia[k]

   print('max and min (gamma_wall) = ', np.max(gamma_wall), np.min(gamma_wall))


def fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, 
               nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,
               gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations):

   #  This is the entry point for fH iteration.
   #  Save 'seed' values for comparison later

   fHs=fH
   nHs=nH
   #________________________________________________________________________________
   # Compute Omega values if nH is non-zero
   #________________________________________________________________________________
   ii = np.where(nH <= 0, count)
   if count <= 0:

      # Compute VxH

      if H_P_EL or H_H2_EL or H_H_EL:
        VxH = np.zeros(nx)
        for k in range(nx):
            VxH[k] = vth * np.sum(Vr2pidVr * (fH[:, :, k] * (vx * dVx))) / nH[k]
      #________________________________________________________________________________
      # Compute Omega_H_P for present fH and Alpha_H_P if H_P elastic collisions are included
      #________________________________________________________________________________
      if H_P_EL:
        if debrief > 1:
            print(prompt + 'Computing Omega_H_P')
        Omega_H_P = np.zeros(nx)
        for k in range(nx):
            DeltaVx = (VxH[k] - vxi[k]) / vth
            MagDeltaVx = np.abs(DeltaVx) > DeltaVx_tol
            DeltaVx = np.sign(DeltaVx) * MagDeltaVx
            Omega_H_P[k] = np.sum(Vr2pidVr * ((Alpha_H_P[:, :, k] * fH[:, :, k]) * dVx)) / (nH[k] * DeltaVx)
        Omega_H_P = Omega_H_P > 0.0
      #________________________________________________________________________________
      # Compute Omega_H_H2 for present fH and Alpha_H_H2 if H_H2 elastic collisions are included
      #________________________________________________________________________________
      if H_H2_EL:
         if debrief > 1:
            print(prompt+'Computing Omega_H_H2')
         for k in range(nx):
            DeltaVx = (VxH[k]-vxH2[k])/vth
            MagDeltaVx = np.abs(DeltaVx) > DeltaVx_tol
            DeltaVx = np.sign(DeltaVx)*MagDeltaVx
            Omega_H_H2[k] = np.sum(Vr2pidVr*np.concatenate(((Alpha_H_H2[:, :, k]*fH[:, :, k]), dVx)))/(nH[k]*DeltaVx)
         Omega_H_H2=Omega_H_H2 > 0.0
      #________________________________________________________________________________
      # Compute Omega_H_H for present fH if H_H elastic collisions are included
      #________________________________________________________________________________
      if H_H_EL:
         if debrief > 1:
            print(prompt+'Computing Omega_H_H')
         if np.sum(MH_H_sum) <= 0:
            for k in range(nx):
               for i in range(nvr):
                  vr2_2vx_ran2[i, :] = vr[i]**2-2*(vx-VxH[k]/vth)**2
                  Wperp_paraH[k] = np.sum(Vr2pidVr*np.concatenate(((vr2_2vx_ran2*fH[:, :, k]), dVx)))/nH[k]
         else:
            for k in range(nx):
               M_fH = MH_H_sum[:, :, k]-fH[:, :, k]
               Wperp_paraH[k] = -np.sum(Vr2pidVr*np.concatenate(((vr2_2vx2_2D*M_fH), dVx)))/nH[k]
         for k in range(nx):
            Work[:] = fH[:,:,k]
            Alpha_H_H[:] = np.concatenate((SIG_H_H, Work))
            Wpp = Wperp_paraH[k]
            MagWpp = np.abs(Wpp) > Wpp_tol
            Wpp = np.sign(Wpp)*MagWpp
            Omega_H_H[k] = np.sum(Vr2pidVr*np.concatenate(((Alpha_H_H*Work), dVx)))/(nH[k]*Wpp)
         Omega_H_H=Omega_H_H > 0.0

      # Total Elastic scattering frequency

      Omega_EL=Omega_H_P+Omega_H_H2+Omega_H_H

      # Total collision frequency

   alpha_c = np.zeros((nvr,nvx,nx))
   if H_P_CX:
      for k in range(nx):
         alpha_c[:, :, k] = alpha_cx[:, :, k] + alpha_ion[k] + Omega_EL[k] + gamma_wall[:, :, k]
   else:
      for k in range(nx):
         alpha_c[:, :, k] = alpha_ion[k] + Omega_EL[k] + gamma_wall[:, :, k]
   
   # Test x grid spacing based on Eq.(27) in notes

   if debrief > 1:
      print(prompt+'Testing x grid spacing')
   Max_dx = np.zeros(nx)
   Max_dx[:] = 1.0e32
   for k in range(nx):
      for j in range(ip[0], nvx): 
         Max_dx[k] = Max_dx[k] < np.min(2*vx[j]/alpha_c[:, j, k])
   dx = np.roll(x, -1) - x
   Max_dxL = Max_dx[:nx-1]
   Max_dxR = Max_dx[1:nx]
   Max_dx = Max_dxL < Max_dxR
   ilarge = np.where(Max_dx < dx[:nx-1], count)
   if count > 0:
      print(prompt+'x mesh spacing is too large!')
      debug=1
      out=''
      jj=0
      print('   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)')
      for ii in range(len(ilarge)+1):
         jj=jj+1
         out = out + ' '.join(['({:3d}) {:9.2E} {:9.2E}'.format(ilarge[ii], x[ilarge[ii]+1] - x[ilarge[ii]], Max_dx[ilarge[ii]]) for ii in range(len(ilarge))])
         if jj > 4:
            print(out)
            jj=0
            out=''
      if jj > 0:
         print(out)
      error=1
      return

   # Define parameters Ak, Bk, Ck, Dk, Fk, Gk

   Ak = np.zeros((nvr,nvx,nx))
   Bk = np.zeros((nvr,nvx,nx))
   Ck = np.zeros((nvr,nvx,nx))
   Dk = np.zeros((nvr,nvx,nx))
   Fk = np.zeros((nvr,nvx,nx))
   Gk = np.zeros((nvr,nvx,nx))

   for k in range(nx-1):
      for j in range(ip[0], (nvx)):
         denom = 2*vx[j]+(x[k+1]-x[k])*alpha_c[:, j, k+1]
         Ak[:, j, k] = (2*vx[j]-(x[k+1]-x[k])*alpha_c[:, j, k])/denom
         Bk[:, j, k] = (x[k+1]-x[k])/denom
         Fk[:, j, k] = (x[k+1]-x[k])*(Sn[:, j, k+1]+Sn[:, j, k])/denom
   for k in range(1, nx):
      for j in range(ip[0]):
         denom = -2*vx[j]+(x[k]-x[k-1])*alpha_c[:, j, k-1]
         Ck[:, j, k] = (-2*vx[j]-(x[k]-x[k-1])*alpha_c[:, j, k])/denom
         Dk[:, j, k] = (x[k]-x[k-1])/denom
         Gk[:, j, k] = (x[k]-x[k-1])*(Sn[:, j, k]+Sn[:, j, k-1])/denom

   # Compute first-flight (0th generation) neutral distribution function

   Beta_CX_sum = np.zeros((nvr,nvx,nx))
   MH_P_sum = np.zeros((nvr,nvx,nx))
   MH_H2_sum = np.zeros((nvr,nvx,nx))
   MH_H_sum = np.zeros((nvr,nvx,nx))
   igen=0
   if debrief > 0:
      print(prompt+'Computing atomic neutral generation#'+sval(igen))
   fHG[:, ip, 0] = fH[:, ip, 0]
   for k in range(nx-1):
      fHG[:, ip, k+1] = fHG[:, ip, k]*Ak[:, ip, k]+Fk[:, ip, k]
   for k in range(nx,1,-1):
      fHG[:,i_n,k-1] = fHG[:,i_n,k] * Ck[:,i_n,k] + Gk[:,i_n,k]

   # Compute first-flight neutral density profile
   #the plotting may not be a perfect translation which may need to be redone -gcheyde

   for k in range(nx):
      NHG[k,igen] = np.sum(Vr2pidVr*np.concatenate((fHG[:, :, k],dVx)))
   if plot > 1:
      fH1d = np.zeros((nvx,nx))
      for k in range(nx):
         fH1d[:, k] = np.concatenate((Vr2pidVr, fHG[:, :, k]))
      plt.plot(vx, fH1d[:, 0], label='First Generation ' + _H)
      for i in range(nx):
        plt.plot(vx, fH1d[:, i], color=(i % 8) + 2)

      plt.ylim([0, np.max(fH1d)])
      plt.title('First Generation ' + _H)
      plt.legend()
      if debug > 0:
        press_return()
      plt.show()

   # Set total atomic neutral distribution function to first flight generation

   fH=fHG
   nH=NHG[:, 0]

   if fH_generations == 0:
      fH_done(plot, NHG, igen, pause, nx, nH, Vr2pidVr, fH, dVx, truncate, nHs, nvr, nvx, H_P_CX, debrief, prompt,Simple_CX, fi)

def next_generation(max_gen, debrief, prompt, nvr, nvx, nx, H_P_CX, Simple_CX, fi_hat, Vr2pidVr, alpha_cx, Work,
                    fHG, ni, SIG_CX, H_H_EL, H_P_EL, H_H2_EL, VxHG, vth, vr2vx2_ran2, vr, vx, ihe, mH, THG, vth2,
                    dVx, NHG, Maxwell, Omega_H_H, vxi, Omega_H_P, vxH2, TH2, Omega_H_H2, ip, Ak, Bk, Ck, i_n,
                    Dk, nHG, plot, debug, Delta_nHs, truncate, Ti):
   if igen+1 > max_gen:
      if debrief > 0:
         print(prompt+'Completed '+sval(max_gen)+' generations. Returning present solution...')
      fH_done(plot, NHG, igen, pause, nx, nH, Vr2pidVr, fH, dVx, truncate, nHs, nvr, nvx, H_P_CX, debrief, prompt,Simple_CX, fi)
   igen=igen+1
   if debrief > 0:
      print(prompt+'Computing atomic neutral generation#'+sval(igen))
   #________________________________________________________________________________
   # Compute Beta_CX from previous generation
   #________________________________________________________________________________
   Beta_CX = np.zeros((nvr,nvx,nx))
   if H_P_CX:
      if debrief > 1:
         print(prompt+'Computing Beta_CX')
      if Simple_CX:

         # Option (B): Compute charge exchange source with assumption that CX source neutrals have
         #             ion distribution function

         for k in range(nx):
            Beta_CX[:, :, k] = fi_hat[:, :, k]*np.sum(Vr2pidVr*np.concatenate(((alpha_cx(*,*,k)*fHG[:,:,k]), dVx)))
      else:

         # Option (A): Compute charge exchange source using fH and vr x sigma x v_v at each velocity mesh point

         for k in range(nx):
            Work[:] = fHG[:, :, k]
            Beta_CX[:, :, k] = ni[k]*fi_hat[:, :, k]*np.concatenate((SIG_CX, Work))

      # Sum charge exchange source over all generations

      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   #________________________________________________________________________________
   # Compute MH from previous generation
   #________________________________________________________________________________

   MH_H = np.zeros((nvr,nvx,nx))
   MH_P = np.zeros((nvr,nvx,nx))
   MH_H2 = np.zeros((nvr,nvx,nx))
   OmegaM = np.zeros((nvr,nvx,nx))
   if H_H_EL or H_P_EL or H_H2_EL:

      # Compute VxHG, THG

      for k in range(nx):
         VxHG[k] = vth*np.sum(Vr2pidVr*np.concatenate((fHG[:,:,k],(vx*dVx))))/NHG[k,igen-1]
         for i in range(nvr):
            vr2vx2_ran2[i, :]=vr[i]**2+(vx-VxHG[k]/vth)**2
         
         if ihe == 1:
           mu_temp=4 # helium mass
           THG[k] = (mu_temp*mH)*vth2*np.sum(Vr2pidVr*(np.concatenate((vr2vx2_ran2*fHG[:, :, k], dVx))))/(3*q*NHG[k,igen-1])
         else:
           THG[k] = (mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2_ran2*fHG[:, :, k]),dVx)))/(3*q*NHG[k,igen-1])
      if H_H_EL:
         if debrief > 1:
            print(prompt+'Computing MH_H')

         # Compute MH_H  

         vx_shift=VxHG
         Tmaxwell=THG
         mol=1
         #@create_shifted_maxwellian.include
         for k in range(nx):
            MH_H[:, :, k] = Maxwell[:, :, k]*NHG[k,igen-1]
            OmegaM[:, :, k] = OmegaM[:, :, k]+Omega_H_H[k]*MH_H[:, :, k]
         MH_H_sum=MH_H_sum+MH_H
      if H_P_EL:
         if debrief > 1:
            print(prompt+'Computing MH_P')

         # Compute MH_P  

         if ihe == 1:
            # kn1d manual page 67, the first line of eq(3.8). I am replacing H -> helium, and H+ -> D+
            m_helium=4
            vx_shift = (m_helium*VxHG+mu*vxi)/(m_helium+mu)
            Tmaxwell = THG+(2.*m_helium*mu/(m_helium+mu)**2)*(Ti-THG +mu*mH*(vxi-VxHG)**2/(6*q))
            mol=1 #1 for atom, 2 for molecule
            mu=4 # for helium
            Shifted_Maxwellian_Debug=0
            #@create_shifted_maxwellian.include
            mu=2 # back to input.  I know it is a deuterium
         else:
            vx_shift=(VxHG+vxi)/2
            Tmaxwell=THG+(2./4.)*(Ti-THG +mu*mH*(vxi-VxHG)**2/(6*q))
            mol=1
            #@create_shifted_maxwellian.include

         for k in range(nx):
            MH_P[:, :, k] = Maxwell[:, :, k]*NHG[k,igen-1]
            OmegaM[:, :, k] = OmegaM[:, :, k]+Omega_H_P[k]*MH_P[:, :, k]
         MH_P_sum=MH_P_sum+MH_P
         print('construct first MH_P_sum')
         
      if H_H2_EL:
         if debrief > 1:
            print(prompt+'Computing MH_H2')

         # Compute MH_H2

         vx_shift = (VxHG+2*vxH2)/3
         Tmaxwell = THG+(4./9.)*(TH2-THG +2*mu*mH*(vxH2-VxHG)**2/(6*q))
         mol=1
         #@create_shifted_maxwellian.include
         for k in range(nx):
            MH_H2[:, :, k] = Maxwell[:, :, k]*NHG[k, igen-1]
            OmegaM[:, :, k] = OmegaM[:, :, k]+Omega_H_H2[k]*MH_H2[:, :, k]
         MH_H2_sum=MH_H2_sum+MH_H2
   #________________________________________________________________________________
   # Compute next generation atomic distribution
   #________________________________________________________________________________

   fHG[:] = 0.0
   for k in range(nx-1):
      fHG[:, ip, k+1] = Ak[:, ip, k]*fHG[:, ip, k] + Bk[:, ip, k]*(Beta_CX[:, ip, k+1]+OmegaM[:, ip, k+1]+Beta_CX[:, ip, k]+OmegaM[:, ip, k])

   for k in range(nx-1,1,-1):
      fHG[:,i_n,k-1]=Ck[:,i_n,k]*fHG[:,i_n,k] + Dk[:,i_n,k]*(Beta_CX[:,i_n,k-1]+OmegaM[:,i_n,k-1]+Beta_CX[:,i_n,k]+OmegaM[:,i_n,k])

   for k in range(nx):
      nHG[k,igen]=np.sum(Vr2pidVr*np.concatenate((fHG[:,:,k], dVx)))

   if plot > 1:
      fH1d=np.zeros((nvx,nx))
      for k in range(nx):
         fH1d[:,k]=np.concatenate((Vr2pidVr, fHG[:,:,k]))
      #probably replacing plot function with pyplot - gcheyde
      plot(vx,fH1d[:,0],/nodata,yrange=[0,max(fH1d)],title=sval(igen)+' Generation '+_H)
      for i in range(nx):
         #likely also need to replace oplot with pyplot - gcheyde
         oplot(vx,(fH1d[:,i] > 0.9),color=(i % 8)+2)
      if debug > 0:
         press_return()

   # Add result to total neutral distribution function

   fH=fH+fHG
   nH=nH+nHG[:,igen]
   #________________________________________________________________________________
   # Compute 'generation error': Delta_nHG=max(NHG(*,igen)/max(nH))
   # and decide if another generation should be computed
   #________________________________________________________________________________
   Delta_nHG=np.max(NHG[:,igen]/np.max(nH))
   if fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations):
      print('Compute generation error, Delta_nHG= ', Delta_nHG,  ',  0.003*Delta_nHs= ', 0.003*Delta_nHs,  ',   truncate=', truncate)

      # If fH 'seed' is being iterated, then do another generation until the 'generation error'
      # is less than 0.003 times the 'seed error' or is less than TRUNCATE

      if (Delta_nHG < 0.003*Delta_nHs) or (Delta_nHG < truncate):
         fH_done(plot, NHG, igen, pause, nx, nH, Vr2pidVr, fH, dVx, truncate, nHs, nvr, nvx, H_P_CX, debrief, prompt,Simple_CX, fi)
      else:

         # If fH 'seed' is NOT being iterated, then do another generation unitl the 'generation error'
         # is less than parameter TRUNCATE

         if Delta_nHG < truncate:
            fH_done(plot, NHG, igen, pause, nx, nH, Vr2pidVr, fH, dVx, truncate, nHs, nvr, nvx, H_P_CX, debrief, prompt,Simple_CX, fi)
      #________________________________________________________________________________

   next_generation(max_gen, debrief, prompt, nvr, nvx, nx, H_P_CX, Simple_CX, fi_hat, Vr2pidVr, alpha_cx, Work,fHG, ni, SIG_CX, H_H_EL, H_P_EL, H_H2_EL, VxHG, vth, vr2vx2_ran2, vr, vx, ihe, mH, THG, vth2,dVx, NHG, Maxwell, Omega_H_H, vxi, Omega_H_P, vxH2, TH2, Omega_H_H2, ip, Ak, Bk, Ck, i_n,Dk, nHG, plot, debug, Delta_nHs, truncate, Ti)

def fH_done(plot, NHG, igen, pause, nx, nH, Vr2pidVr, fH, dVx, truncate, nHs, nvr, nvx, H_P_CX, debrief, prompt,
            Simple_CX, fi):

   if plot > 0:
      plot(x,NHG[:,0],/nodata,yrange=[np.max(NHG)*truncate,np.max(NHG)],title=_H+' Density by Generation',xtitle='x (m)',ytitle='Density (m!U-3!N)',/ylog)
      for i in range(igen+1):
         oplot(x,NHG[:,i],color=i+1)
      print('Total number of genrations = ', igen)
      if pause:
         press_return()

   # Compute H density profile
   for k in range(nx):
      nH[k]=np.sum(Vr2pidVr*np.concatenate((fH[:,:,k],dVx)))

   if fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations):
      #________________________________________________________________________________
      # Compute 'seed error': Delta_nHs=(|nHs-nH|)/max(nH) 
      # If Delta_nHs is greater than 10*truncate then iterate fH
      #________________________________________________________________________________
      print('Compute Seed Error, Delta_nHS = ', Delta_nHs, ',    10*truncate= ',  10*truncate)
      Delta_nHs=np.max(np.abs(nHs-nH))/np.max(nH)
      if Delta_nHs > 10*truncate:
         fH_Iterate(count, H_P_EL, H_H2_EL, H_H_EL, nx, vth, Vr2pidVr, vx, dVx, debrief, prompt, vxi, DeltaVx_tol, Alpha_H_P, vxH2, Alpha_H_H2, nvr, vr2_2vx_ran2, Wperp_paraH, vr, vr2_2vx2_2D, Work, Alpha_H_H, SIG_H_H, Wpp_tol, nvx, H_P_CX, alpha_cx, alpha_ion,gamma_wall, Sn, ip, fHG, NHG, plot, i_n, _H, fH_generations)
   #________________________________________________________________________________
   # Update Beta_CX_sum using last generation
   #________________________________________________________________________________
   Beta_CX=np.zeros((nvr,nvx,nx))
   if H_P_CX:
      if debrief > 1:
         print(prompt+'Computing Beta_CX')
      if Simple_CX:

         # Option (B): Compute charge exchange source with assumption that CX source neutrals have
         #             ion distribution function

         for k in range(nx):
            Beta_CX[:,:,k]=fi_hat[:,:,k]*np.sum(Vr2pidVr*np.concatenate(((alpha_cx[:,:,k]*fHG[:,:,k]),dVx)))

      else:

         # Option (A): Compute charge exchange source using fH and vr x sigma x v_v at each velocity mesh point

         for k in range(nx):
            Work[:]=fHG[:,:,k]
            Beta_CX[:,:,k]=ni[k]*fi_hat[:,:,k]*np.concatenate((SIG_CX,Work))

      # Sum charge exchange source over all generations

      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   #________________________________________________________________________________
   # Update MH_*_sum using last generation
   #________________________________________________________________________________

   MH_H=np.zeros((nvr,nvx,nx))
   MH_P=np.zeros((nvr,nvx,nx))
   MH_H2=np.zeros((nvr,nvx,nx))
   OmegaM=np.zeros((nvr,nvx,nx))
   if H_H_EL or H_P_EL or H_H2_EL:

      # Compute VxHG, THG

      for k in range(nx):
         VxHG[k]=vth*np.sum(Vr2pidVr*np.concatenate((fHG[:,:,k],(vx*dVx))))/NHG[k,igen]
         for i in range(nvr):
            vr2vx2_ran2[i,:]=vr[i]**2+(vx-VxHG[k]/vth)**2
         if ihe == 1: #sgbaek, check for mu
            mu_helium=4
            THG[k]=(mu_helium*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2_ran2*fHG[:,:,k]),dVx)))/(3*q*NHG[k,igen])
         else:
            THG[k]=(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2_ran2*fHG[:,:,k]),dVx)))/(3*q*NHG[k,igen])
      if H_H_EL:
         if debrief > 1:
            print(prompt+'Computing MH_H')

         # Compute MH_H  

         vx_shift=VxHG
         Tmaxwell=THG
         mol=1
         #@create_shifted_maxwellian.include
         for k in range(nx):
            MH_H[:,:,k]=Maxwell[:,:,k]*NHG[k,igen]
            OmegaM[:,:,k]=OmegaM[:,:,k]+Omega_H_H[k]*MH_H[:,:,k]
         MH_H_sum=MH_H_sum+MH_H
      if H_P_EL:
         if debrief > 1:
            print(prompt+'Computing MH_P')

         # Compute MH_P  

         if ihe == 1:
            # kn1d manual page 67, the first line of eq(3.8). I am replacing H -> helium, and H+ -> D+
            m_helium=4
            vx_shift=(m_helium*VxHG+mu*vxi)/(m_helium+mu)
            Tmaxwell=THG+(2.*m_helium*mu/(m_helium+mu)**2)*(Ti-THG +mu*mH*(vxi-VxHG)**2/(6*q))
            mol=m_helium
            mu=4
            Shifted_Maxwellian_Debug=1
            #@create_shifted_maxwellian.include
            print('FINAL UPDATE MH_P')
            mu=2
         else:
            vx_shift=(VxHG+vxi)/2
            Tmaxwell=THG+(2./4.)*(Ti-THG +mu*mH*(vxi-VxHG)**2/(6*q))
            mol=1
            #@create_shifted_maxwellian.include
         for k in range(nx):
            MH_P[:,:,k]=Maxwell[:,:,k]*NHG[k,igen]
            OmegaM[:,:,k]=OmegaM[:,:,k]+Omega_H_P[k]*MH_P[:,:,k]
         MH_P_sum=MH_P_sum+MH_P
         print('second MH_P_SUM')
      if H_H2_EL:
         if debrief > 1:
            print(prompt+'Computing MH_H2')

         # Compute MH_H2

         vx_shift=(VxHG+2*vxH2)/3
         Tmaxwell=THG+(4./9.)*(TH2-THG +2*mu*mH*(vxH2-VxHG)**2/(6*q))
         mol=1
         #@create_shifted_maxwellian.include
         for k in range(nx):
            MH_H2[:,:,k]=Maxwell[:,:,k]*NHG[k,igen]
            OmegaM[:,:,k]=OmegaM[:,:,k]+Omega_H_H2[k]*MH_H2[:,:,k]
         MH_H2_sum=MH_H2_sum+MH_H2
   #________________________________________________________________________________
   # Compute remaining moments
   #________________________________________________________________________________

   # GammaxH - particle flux in x direction
   for k in range(nx):
      GammaxH[k]=vth*np.sum(Vr2pidVr*np.concatenate((fH(*,*,k),(vx*dVx))))

   # VxH - x velocity
   VxH=GammaxH/nH
   _VxH=VxH/vth

   # magnitude of random velocity at each mesh point
   vr2vx2_ran=np.zeros((nvr,nvx,nx))
   for i in range(nvr):
      for k in range(nx):
         vr2vx2_ran[i,:,k]=vr[i]**2+(vx-_VxH[k])**2

   # pH - pressure 
   #sgbaek check for mu
   for k in range(nx):
      pH[k]=(mu*mH)*vth2*total(Vr2pidVr*np.concatenate(((vr2vx2_ran[:,:,k]*fH[:,:,k]),dVx)))/(3*q)
   if ihe == 1:
      pH[:]=pH[:]*4/mu

   # TH - temperature
   TH=pH/nH

   # piH_xx
   if ihe == 1:
      mu_He=4
      for k in range(nx):
         piH_xx[k]=(mu_He*mH)*vth2*np.sum(Vr2pidVr*np.concatenate((fH[:,:,k],(dVx*(vx-_VxH[k])**2))))/q - pH[k]
   else:
      for k in range(nx):
         piH_xx[k]=(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate((fH[:,:,k],(dVx*(vx-_VxH[k])**2))))/q - pH[k]

   

   # piH_yy
   if ihe == 1:
      mu_He=4
      for k in range(nx):
         piH_yy[k]=(mu_He*mH)*vth2*0.5*np.sum((Vr2pidVr*vr**2)*(fH[:,:,k],dVx))/q - pH[k]
   else:
      for k in range(nx):
         piH_yy[k]=(mu*mH)*vth2*0.5*np.sum((Vr2pidVr*vr**2)*np.concatenate((fH[:,:,k],dVx)))/q - pH[k]


   # piH_zz
   piH_zz=piH_yy

   # qxH
   for k in range(nx):
      qxH[k]=0.5*(mu*mH)*vth3*np.sum(Vr2pidVr*np.concatenate(((vr2vx2_ran[:,:,k]*fH[:,:,k]),(dVx*(vx-_VxH[k])))))
   if ihe == 1:
      qxH[:]=qxH[:]*4/mu
   #________________________________________________________________________________
   # C = RHS of Boltzman equation for total fH
   #________________________________________________________________________________
   for k in range(nx):
      C=vth*(Sn[:,:,k] + Beta_CX_sum[:,:,k] - alpha_c[:,:,k]*fH[:,:,k] + Omega_H_P[k]*MH_P_sum[:,:,k]+Omega_H_H2[k]*MH_H2_sum[:,:,k]+Omega_H_H[k]*MH_H_sum[:,:,k])
      #sgbaek - check for mu
      QH[k]=0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2_ran[:,:,k]*C),dVx)))
      RxH[k]=(mu*mH)*vth*np.sum(Vr2pidVr*np.concatenate((C,(dVx*(vx-_VxH[k])))))
      if ihe == 1:
         QH[k]=QH[k]*4/mu 
         RxH[k]=RxH[k]*4/mu
      NetHSource[k]=np.sum(Vr2pidVr*np.concatenate((C,dVx)))
      Sion[k]=vth*nH[k]*alpha_ion[k]
      SourceH[k]=np.sum(Vr2pidVr*np.concatenate((fSH[:,:,k],dVx)))
      WallH[k]=vth*np.sum(Vr2pidVr*np.concatenate(((gamma_wall[:,:,k]*fH[:,:,k]),dVx)))
      if Recomb:
         SRecomb[k]=vth*ni[k]*Rec[k]
      else:
         SRecomb[k]=0.0
      if H_P_CX:
         CCX=vth*(Beta_CX_sum[:,:,k] - alpha_cx[:,:,k]*fH[:,:,k])
         RxHCX[k]=(mu*mH)*vth*np.sum(Vr2pidVr*np.concatenate((CCX,(dVx*(vx-_VxH[k])))))
         EHCX[k]=0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*CCX),dVx)))
      if H_H2_EL:
         CH_H2=vth*Omega_H_H2[k]*(MH_H2_sum[:,:,k]-fH[:,:,k])
         RxH2_H[k]=(mu*mH)*vth*np.sum(Vr2pidVr*np.concatenate((CH_H2,(dVx*(vx-_VxH[k])))))
         EH2_H[k]=0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*CH_H2),dVx)))
      if H_P_EL:
         CH_P=vth*Omega_H_P[k]*(MH_P_sum[:,:,k]-fH[:,:,k])
         #sgbaek check for mu
         RxP_H[k]=(mu*mH)*vth*np.sum(Vr2pidVr*np.concatenate((CH_P,(dVx*(vx-_VxH[k])))))
         if ihe == 1:
            RxP_H[k]=RxP_H[k]*4/mu
         EP_H[k]=0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate((((vr2vx2[:,:,k]*CH_P),dVx))))
         if ihe == 1:
            EP_H[k]=EP_H[k]*4/mu
         print('k,rxp_h',k,rxp_h[k],ep_h[k])
      CW_H=-vth*(gamma_wall[:,:,k]*fH[:,:,k])
      # sgbaek check for mu
      RxW_H[k]=(mu*mH)*vth*np.sum(Vr2pidVr*np.concatenate((CW_H,(dVx*(vx-_VxH[k])))))
      if ihe == 1:
         RxW_H[k]=RxW_H[k]*4/mu
      EW_H[k]=0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*CW_H),dVx)))
      if ihe == 1:
         EW_H[k]=EW_H[k]*4/mu
      if H_H_EL:
         CH_H=vth*Omega_H_H[k]*(MH_H_sum[:,:,k]-fH[:,:,k])
         for i in range(nvr):
            vr2_2vx_ran2[i,:]=vr[i]**2-2*(vx-_VxH[k])**2
         Epara_PerpH_H[k]=-0.5*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((vr2_2vx_ran2*CH_H),dVx)))
   print('before qxH')



   # qxH_total

   if ihe == 1:
      mu_he=4
      qxH_total=(0.5*nH*(mu_he*mH)*VxH*VxH + 2.5*pH*q)*VxH + q*piH_xx*VxH + qxH
   else:
      qxH_total=(0.5*nH*(mu*mH)*VxH*VxH + 2.5*pH*q)*VxH + q*piH_xx*VxH + qxH


   # QH_total

   if ihe == 1:
      mu_he=4
      QH_total=QH+RxH*VxH + 0.5*(mu_he*mH)*NetHSource*VxH*VxH
   else:
      QH_total=QH+RxH*VxH + 0.5*(mu*mH)*NetHSource*VxH*VxH


   # Albedo

   AlbedoH=0.0
   gammax_plus=vth*np.sum(Vr2pidVr*np.concatenate((fH[:,ip,0],(vx[ip]*dVx[ip]))))
   gammax_minus=vth*np.sum(Vr2pidVr*np.concatenate((fH[:,i_n,0],(vx[i_n]*dVx[i_n]))))
   if np.abs(gammax_plus) > 0:
      AlbedoH=-gammax_minus/gammax_plus

   #________________________________________________________________________________
   # Compute Mesh Errors
   #________________________________________________________________________________
   mesh_error=np.zeros((nvr,nvx,nx))
   max_mesh_error=0.0
   min_mesh_error=0.0
   mtest=5
   moment_error=np.zeros((nx,mtest))
   max_moment_error=np.zeros(mtest)
   C_error=np.zeros(nx)
   CX_error=np.zeros(nx)
   H_H_error=np.zeros((nx,3))
   H_H2_error=np.zeros((nx,3))
   H_P_error=np.zeros((nx,3))
   max_H_H_error=np.zeros(3)
   max_H_H2_error=np.zeros(3)
   max_H_P_error=np.zeros(3)

   if compute_errors:
      if debrief > 1:
         print(prompt+'Computing Collision Operator, Mesh, and Moment Normalized Errors')

      NetHSource2=SourceH+SRecomb-Sion-WallH
      for k in range(nx):
         C_error[k]=np.abs(NetHSource[k]-NetHSource2[k])/np.max(np.abs([NetHSource[k],NetHSource2[k]]))

      # Test conservation of particles for charge exchange operator

      if H_P_CX:
         for k in range(nx):
            CX_A=np.sum(Vr2pidVr*np.concatenate(((alpha_cx[:,:,k]*fH[:,:,k]),dVx)))
            CX_B=np.sum(Vr2pidVr*np.concatenate(((Beta_CX_sum[:,:,k]),dVx)))
            CX_error[k]=np.abs(CX_A-CX_B)/np.max(np.abs([CX_A,CX_B]))

      # Test conservation of particles, x momentum, and total energy of elastic collision operators

      for m in range(3):
         for k in range(nx):
            if m < 2:
               TfH=np.sum(Vr2pidVr*np.concatenate((fH[:,:,k],(dVx*vx**m))))
            else:
               TfH=np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*fH[:,:,k]),dVx)))
            if H_H_EL:
               if m < 2:
                  TH_H=np.sum(Vr2pidVr*np.concatenate((MH_H_sum[:,:,k],(dVx*vx**m))))
               else:
                  TH_H=np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*MH_H_sum[:,:,k]),dVx)))
               H_H_error[k,m]=np.abs(TfH-TH_H)/np.max(np.abs([TfH,TH_H]))
            if H_H2_EL:
               if m < 2:
                  TH_H2=np.sum(Vr2pidVr*np.concatenate((MH_H2_sum[:,:,k],(dVx*vx**m))))
               else:
                  TH_H2=np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*MH_H2_sum[:,:,k]),dVx)))
               H_H2_error[k,m]=np.abs(TfH-TH_H2)/np.max(np.abs([TfH,TH_H2]))
            if H_P_EL:
               if m < 2:
                  TH_P=np.sum(Vr2pidVr*np.concatenate((MH_P_sum[:,:,k],(dVx*vx**m))))
               else:
                  TH_P=np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*MH_P_sum[:,:,k]),dVx)))
               H_P_error[k,m]=np.abs(TfH-TH_P)/np.max(np.abs([TfH,TH_P]))
               print('h_p_error= ',k, m, TfH-TH_p, TfH, TH_p)
         max_H_H_error[m]=np.max(H_H_error[:,m])
         max_H_H2_error[m]= np.max(H_H2_error[:,m])
         max_H_P_error[m]= np.max(H_P_error[:,m])

      if CI_test:

         # Compute Momentum transfer rate via full collision integrals for charge exchange and mixed elastic scattering
         # Then compute error between this and actual momentum transfer resulting from CX and BKG (elastic) models

         if H_P_CX: #P -> H charge exchange momentum transfer via full collision integral
            print(prompt+'Computing P -> H Charge Exchange Momentum Transfer')
            _Sig=np.zeros((nvr*nvx*nvr*nvx,ntheta))
            _Sig[:]=v_v*sigma_cx_h0(v_v2*(0.5*mH*vth2/q))
            SIG_VX_CX=np.zeros((nvr*nvx,nvr*nvx))
            SIG_VX_CX[:]=Vr2pidVrdVx*vx_vx*np.concatenate((_Sig,dtheta))
            alpha_vx_cx=np.zeros((nvr,nvx,nx))
            for k in range(nx):
               Work[:]=fi_hat[:,:,k]*ni[k]
               alpha_vx_cx[:,:,k]=np.concatenate((SIG_VX_CX,Work))
            for k in range(nx):
               RxCI_CX[k]=-(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((Alpha_vx_cx[:,:,k]*fH[:,:,k]),dVx)))
            norm=np.max(np.abs([RxHCX,RxCI_CX]))
            for k in range(nx):
               CI_CX_error[k]=np.abs(RxHCX[k]-RxCI_CX[k])/norm
            print(prompt+'Maximum normalized momentum transfer error in CX collision operator: ',sval(np.max(CI_CX_Error)))
         if H_P_EL: #P -> H momentum transfer via full collision integral
            if ihe == 1: 
                mu_helium=4
                for k in range(nx):
                  RxCI_P_H[k]=-(1.0/2.0)*(mu_helium*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((Alpha_H_P[:,:,k]*fH[:,:,k]),dVx)))
            else:
                for k in range(nx):
                  RxCI_P_H[k]=-(1.0/2.0)*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((Alpha_H_P[:,:,k]*fH[:,:,k]),dVx)))
            norm=np.max(np.abs([RxP_H,RxCI_P_H]))
            for k in range(nx):
               CI_P_H_error[k]=np.abs(RxP_H[k]-RxCI_P_H[k])/norm
            print(prompt+'Maximum normalized momentum transfer error in P -> H elastic BKG collision operator: ',sval(np.max(CI_P_H_Error)))
         if H_H2_EL: #H2 -> H momentum transfer via full collision integral
            for k in range(nx):
               RxCI_H2_H[k]=-(2.0/3.0)*(mu*mH)*vth2*np.sum(Vr2pidVr*np.concatenate(((Alpha_H_H2[:,:,k]*fH[:,:,k]),dVx)))
            norm=np.max(np.abs([RxH2_H,RxCI_H2_H]))
            for k in range(nx):
               CI_H2_H_error[k]=np.abs(RxH2_H[k]-RxCI_H2_H[k])/norm
            print(prompt+'Maximum normalized momentum transfer error in H2 -> H elastic BKG collision operator: ',sval(np.max(CI_H2_H_Error)))
         if H_H_EL: #H -> H perp/parallel energy transfer via full collision integral
            for k in range(nx):
               Work[:]=fH[:,:,k]
               Alpha_H_H[:]=np.concatenate((SIG_H_H,Work))
               Epara_Perp_CI[k]=0.5*(mu*mH)*vth3*np.sum(Vr2pidVr*np.concatenate(((Alpha_H_H*fH[:,:,k]),dVx)))
            norm=np.max(np.abs([Epara_PerpH_H,Epara_Perp_CI]))
            for k in range(nx):
               CI_H_H_error[k]=np.abs(Epara_PerpH_H[k]-Epara_Perp_CI[k])/norm
            print(prompt+'Maximum normalized perp/parallel energy transfer error in H -> H elastic BKG collision operator: ',sval(np.max(CI_H_H_Error)))

      #  Mesh Point Error based on fH satisfying Boltzmann equation

      T1=np.zeros((nvr,nvx,nx))
      T2=T1
      T3=T1
      T4=T1
      T5=T1
      for k in range(nx-1):
         for j in range(nvx):
            T1[:,j,k]=2*vx([j])*(fH[:,j,k+1]-fH[:,j,k])/(x[k+1]-x[k])
         T2[:,:,k]=Sn[:,:,k+1]+Sn[:,:,k]
         T3[:,:,k]=Beta_CX_sum[:,:,k+1]+Beta_CX_sum[:,:,k]
         T4[:,:,k]=alpha_c[:,:,k+1]*fH[:,:,k+1]+alpha_c[:,:,k]*fH[:,:,k]
         T5[:,:,k]=Omega_H_P[k+1]*MH_P_sum[:,:,k+1]+Omega_H_H2[k+1]*MH_H2_sum[:,:,k+1]+Omega_H_H[k+1]*MH_H_sum[:,:,k+1]+ Omega_H_P[k]*  MH_P_sum[:,:,k]+  Omega_H_H2[k]*  MH_H2_sum[:,:,k]+  Omega_H_H[k]*  MH_H_sum[:,:,k]
         Mesh_Error[:,:,k]=np.abs(T1[:,:,k]-T2[:,:,k]-T3[:,:,k]+T4[:,:,k]-T5[:,:,k])/np.max(np.abs([T1[:,:,k],T2[:,:,k],T3[:,:,k],T4[:,:,k],T5[:,:,k]]))
      ave_mesh_error=np.sum(mesh_error)/len(mesh_error)
      max_mesh_error=np.max(mesh_error)
      min_mesh_error=np.min(mesh_error[:,:,0:nx-2])

      #  Moment Error

      for m in range(mtest):
         for k in range(nx-1):
            MT1=np.sum(Vr2pidVr*np.concatenate((T1[:,:,k],(dVx*vx**m))))
            MT2=np.sum(Vr2pidVr*np.concatenate((T2[:,:,k],(dVx*vx**m))))
            MT3=np.sum(Vr2pidVr*np.concatenate((T3[:,:,k],(dVx*vx**m))))
            MT4=np.sum(Vr2pidVr*np.concatenate((T4[:,:,k],(dVx*vx**m))))
            MT5=np.sum(Vr2pidVr*np.concatenate((T5[:,:,k],(dVx*vx**m))))
            moment_error[k,m]=np.abs(MT1-MT2-MT3+MT4-MT5)/np.max(np.abs([MT1,MT2,MT3,MT4,MT5]))
         max_moment_error[m]=np.max(moment_error[:,m])

      # Compute error in qxH_total

      #    qxH_total2 total neutral heat flux profile (watts m^-2)
      #               This is the total heat flux transported by the neutrals
      #               computed in a different way from:

      #               qxH_total2(k)=vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#(vx*dVx)))*0.5*(mu*mH)

      #               This should agree with qxH_total if the definitions of nH, pH, piH_xx,
      #               TH, VxH, and qxH are coded correctly.

      qxH_total2=np.zeros(nx)
      if ihe == 1:
         mu_helium=4
         for k in range(nx):
            qxH_total2[k]=0.5*(mu_helium*mH)*vth3*np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*fH[:,:,k]),(vx*dVx)))) 
      else:
         for k in range(nx):
            qxH_total2[k]=0.5*(mu*mH)*vth3*np.sum(Vr2pidVr*np.concatenate(((vr2vx2[:,:,k]*fH[:,:,k]),(vx*dVx))))
      qxH_total_error=np.abs(qxH_total-qxH_total2)/np.max(np.abs([qxH_total,qxH_total2]))

      # Compute error in QH_total

      Q1=np.zeros(nx)
      Q2=np.zeros(nx)
      QH_total_error=np.zeros(nx)
      for k in range(nx-1):
         Q1[k]=(qxH_total[k+1]-qxH_total[k])/(x[k+1]-x[k])
         Q2[k]=0.5*(QH_total[k+1]+QH_total[k])
         print('k , Q1(k), Q2(k) = ', k, Q1(k), Q2(k))
      QH_total_error=np.abs(Q1-Q2)/np.max(np.abs([Q1,Q2]))

      if debrief > 0:
         print(prompt+'Maximum particle convervation error of total collision operator: ',sval(np.max(C_Error)))
         print(prompt+'Maximum H_P_CX  particle convervation error: ',sval(np.max(CX_Error)))
         print(prompt+'Maximum H_H_EL  particle conservation error: '+sval(max_H_H_error[0]))
         print(prompt+'Maximum H_H_EL  x-momentum conservation error: '+sval(max_H_H_error[1]))
         print(prompt+'Maximum H_H_EL  total energy conservation error: '+sval(max_H_H_error[2]))
         print(prompt+'Maximum H_H2_EL particle conservation error: '+sval(max_H_H2_error[0]))
         print(prompt+'Maximum H_P_EL  particle conservation error: '+sval(max_H_P_error[0]))
         print(prompt+'Average mesh_error =',ave_mesh_error)
         print(prompt+'Maximum mesh_error =',max_mesh_error)
         for m in range(5):
            print(prompt+'Maximum fH vx^'+sval(m)+' moment error: '+sval(max_moment_error[m]))
         print(prompt+'Maximum qxH_total error =',np.max(qxH_total_Error))
         print(prompt+'Maximum QH_total error =',np.max(QH_total_Error))
         if debug > 0:
            press_return()

   if plot > 1:
      fH1d=np.zeros((nvx,nx))
      for k in range(nx):
         fH1d[:,k]=np.concatenate((Vr2pidVr,fH[:,:,k]))
      plot(vx,fH1d[:,0],yrange=[np.min(fH1d),np.max(fH1d)],/nodata,title=_H+' Velocity Distribution Function: fH(vx)',xtitle='vx/Vth')
      for i in range(nx):
         oplot(x,fH1d[:,i],color=(i % 6)+2)
      if pause:
         press_return()

   mid1=np.where(x == 0.7*(np.max(x)+np.min(x))/2)[0]
   mid2=np.where(x == 0.85*(np.max(x)+np.min(x))/2)[0]
   mid3=np.where(x == (np.max(x)+np.min(x))/2)[0]
   mid4=np.where(x == 1.15*(np.max(x)+np.min(x))/2)[0]
   mid5=np.where(x == 1.3*(np.max(x)+np.min(x))/2)[0]
   if plot > 0:
      data=[nH,n,nHP,nH2]
      jp=np.where(data > 0)[0]
      yrange=[np.min(data[jp]),np.max(data[jp])]
      plot(x,nH,/nodata,/ylog,yrange=yrange,title='Density Profiles',xtitle='x (meters)',ytitle='m!U-3!N')
      oplot(x,nH,color=2)
      xyouts(x[mid1],1.2*nH[mid1],_H,color=2)
      oplot(x,n,color=3)
      xyouts(x[mid2],1.2*n[mid2],'e!U-!N',color=3)
      oplot(x,nH2,color=4)
      xyouts(x[mid3],1.2*nH2[mid3],_HH,color=4)
      oplot(x,nHP,color=6)
      xyouts(x[mid4],1.2*nHP[mid4],_Hp,color=6)
      if pause:
         press_return()

   if plot > 0:
      data=[TH,Te,Ti,THP,TH2]
      jp=np.where(data > 0)[0]
      yrange=[np.min(data[jp]),np.max(data[jp])]
      plot(x,TH,/nodata,/ylog,yrange=yrange,title='Temperature Profiles',xtitle='x (meters)',ytitle='eV')
      oplot(x,TH,color=2)
      xyouts(x[mid1],1.2*TH[mid1],_H,color=2)
      oplot(x,Te,color=3)
      xyouts(1.1*x[mid2],1.2*Ti[mid2],_P,color=5)
      oplot(x,Te,color=5)
      xyouts(x[mid3],1.2*Te[mid3],'e!U-!N',color=3)
      oplot(x,TH2,color=4)
      xyouts(x[mid4],1.2*TH2[mid4],_HH,color=4)
      oplot(x,THP,color=6)
      xyouts(x[mid5],1.2*THP[mid5],_Hp,color=6)
      if pause:
         press_return()

   if plot > 0:
      data=[SourceH,SRecomb,Sion]
      jp=np.where(data > 0)[0]
      yrange=[np.min(data[jp]),np.max(data[jp])]
      plot(x,SourceH,/nodata,/ylog,yrange=yrange,title=_H+' Source and Sink Profiles',xtitle='x (meters)',ytitle='m!U-3!N s!U-1!N')
      oplot(x,SourceH,color=2)
      xyouts(x[mid1],1.2*SourceH[mid1],_HH+' Dissociation',color=2)
      oplot(x,SRecomb,color=3)
      xyouts(x[mid2],1.2*SRecomb[mid2],_P+' Recombination',color=3)
      oplot(x,Sion,color=4)
      xyouts(x[mid3],1.2*Sion[mid3],_H+' Ionization',color=4)
      oplot(x,WallH,color=6)
      xyouts(x[mid4],1.2*WallH[mid4],_H+' Side-Wall Loss',color=6)
      if pause:
         press_return()

   if plot > 0:
      gammaxH_plus=np.zeros(nx)
      gammaxH_minus=np.zeros(nx)
      for k in range(nx):
         gammaxH_plus[k]=vth*np.sum(Vr2pidVr*np.concatenate((fH[:,ip,k],(vx[ip]*dVx[ip]))))
         gammaxH_minus[k]=vth*np.sum(Vr2pidVr*np.concatenate((fH[:,i_n,k],(vx[i_n]*dVx[i_n]))))
      data=[gammaxH_plus,gammaxH_minus,gammaxH]
      jp=np.where(data < 1.0e32)[0]
      yrange=[np.min(data[jp]),np.max(data[jp])]
      plot(x,gammaxH,/nodata,yrange=yrange,title=_H+' Fluxes',xtitle='x (meters)',ytitle='m!U-2!N s!U-1!N')
      oplot(x,gammaxH,color=2)
      xyouts(x[mid1],GammaxH[mid1],'!7C!5',color=2)
      oplot(x,gammaxH_plus,color=3)
      xyouts(x[mid2],GammaxH_plus[mid2],'!7C!5!U+!N',color=3)
      oplot(x,GammaxH_minus,color=4)
      xyouts(x[mid3],GammaxH_minus[mid3],'!7C!5!U-!N',color=4)
      if pause:
         press_return()



   xloc=[.15,.3,.45,.6,.75,.95]
   mid=xloc
   midH=xloc
   midH2=xloc

   xH=x

   # Momentum Transfer

   if plot > 0:
      ydata=[RxHCX,RxH2_H,RxP_H,RxW_H]
      yrange=[np.min(ydata),np.max(ydata)]
      plot(x,n,/nodata,yrange=yrange,title='KN1D: x-Momentum Transfer Rates',xtitle='x (meters)',ytitle='nt m!U-3!N')
      #@kn1d.include
      oplot(xH,RxHCX,color=4)
      xyouts(xH(midH[4]),RxHCX(midH[4]),_P+'->'+_H+'(CX)',color=4)
      oplot(xH,RxH2_H,color=3)
      xyouts(xH(midH[1]),1.2*RxH2_H(midH[1]),_HH+'->'+_H,color=3)
      oplot(xH,RxP_H,color=2)
      xyouts(xH(midH[5]),RxP_H(midH[5]),_P+'->'+_H+'(EL)',color=2)
      oplot(xH,-RxW_H,color=3)
      xyouts(xH(midH[0]),-1.2*RxW_H(midH[0]),_H+'->Side Wall',color=3)
      #@kn1d_limiter.include
      if pause:
         press_return()

   print('plot momentum trsnfer')

 

   # Energy Tranfer

   if plot > 0:
      ydata=[EHCX,EH2_H,EP_H,EW_H]
      f=1/1.0e6
      yrange=[np.min(ydata),np.max(ydata)]*f
      plot(x,n*f,/nodata,yrange=yrange,title='KN1D: Energy Transfer Rates',xtitle='x (meters)',ytitle='MW m!U-3!N')
      #@kn1d.include
      oplot(xH,EHCX*f,color=4)
      xyouts(xH(midH[5]),EHCX(midH[5])*f,_P+'->'+_H+'(CX)',color=4,/norm)
      oplot(xH,EH2_H*f,color=3)
      xyouts(xH(midH[4]),EH2_H(midH[4])*f,_HH+'->'+_H,color=3,/norm)
      oplot(xH,EP_H*f,color=2)
      xyouts(xH(midH[4]),EP_H(midH[4])*f,_P+'->'+_H+'(EL)',color=2,/norm)
      oplot(xH,-EW_H*f,color=9)
      xyouts(xH(midH[0]),-EW_H(midH[0])*f,_H+'Side Wall',color=9,/norm)
      #@kn1d_limiter.include
      if pause:
         press_return()


   # Temperature Isotropization

   if plot > 0:
      ydata=[Epara_PerpH_H]
      yrange=[np.min(ydata),np.max(ydata)]
      plot(x,n,/nodata,yrange=yrange,title='KN1D: T!D//!N -> T!D!9x!3!N Isotropization Rates',xtitle='x (meters)',ytitle='watts m!U-3!N')
      #@kn1d.include
      oplot(xH,Epara_PerpH_H,color=4)
      xyouts(xH(midH[3]),Epara_PerpH_H(midH[3]),_H+'<->'+_H+'(EL)',color=4)
      #@kn1d_limiter.include
      if pause:
         press_return()


   # Heat Fluxes

   if plot > 0:
      ydata=[qxH_total]
      f=1/1.0e3
      yrange=[np.min(ydata),np.max(ydata)]*f
      plot(x,n*f,/nodata,yrange=yrange,title='KN1D: Heat Fluxes',xtitle='x (meters)',ytitle='KW m!U-2!N')
      #@kn1d.include
      oplot(xH,qxH_total*f,color=4)
      xyouts(xH(midH[3]),qxH_total(midH[3])*f,_H,color=4)
      #@kn1d_limiter.include
      if pause:
         press_return()



   # Save input parameters in common block

   g.Kinetic_H_input_vx_s = vx
   g.Kinetic_H_input_vr_s = vr
   g.Kinetic_H_input_x_s = x
   g.Kinetic_H_input_Tnorm_s = Tnorm
   g.Kinetic_H_input_mu_s = mu
   g.Kinetic_H_input_Ti_s = Ti
   g.Kinetic_H_input_vxi_s = vxi
   g.Kinetic_H_input_Te_s = Te
   g.Kinetic_H_input_n_s = n
   g.Kinetic_H_input_vxi_s = vxi
   g.Kinetic_H_input_fHBC_s = fHBC
   g.Kinetic_H_input_GammaxHBC_s = GammaxHBC
   g.Kinetic_H_input_PipeDia_s = PipeDia
   g.Kinetic_H_input_fH2_s = fH2
   g.Kinetic_H_input_fSH_s = fSH
   g.Kinetic_H_input_nHP_s = nHP
   g.Kinetic_H_input_THP_s = THP
   g.Kinetic_H_input_fH_s = fH
   g.Kinetic_H_input_Simple_CX_s = Simple_CX
   g.Kinetic_H_input_JH_s = JH
   g.Kinetic_H_input_Recomb_s = Recomb
   g.Kinetic_H_input_H_H_EL_s = H_H_EL
   g.Kinetic_H_input_H_P_EL_s = H_P_EL
   g.Kinetic_H_input_H_H2_EL_s = H_H2_EL
   g.Kinetic_H_input_H_P_CX_s = H_P_CX

# this part was commented out in the original, so it can probably be removed later - gcheyde
# Set output parameters to single precision
 
#   fH=float(fH)
#   nH=float(nH)
#   GammaxH=float(GammaxH)
#   VxH=float(VxH)
#   pH=float(pH)
#   TH=float(TH)
#   qxH=float(qxH)
#   qxH_total=float(qxH_total)
#   NetHSource=float(NetHSource)
#   Sion=float(Sion)
#   QH=float(QH)
#   RxH=float(RxH)
#   QH_total=float(QH_total)
#   AlbedoH=float(AlbedoH)
#   piH_xx=float(piH_xx)
#   piH_yy=float(piH_yy)
#   piH_zz=float(piH_zz)
#   RxH2_H=float(RxH2_H)
#   RxP_H=float(RxP_H)
#   RxHCX=float(RxHCX)
#   EH2_H=float(EH2_H)
#   EP_H=float(EP_H)
#   EHCX=float(EHCX)
#   Epara_PerpH_H=float(Epara_PerpH_H)

#   vr=float(vr)
#   vx=float(vx)
#   x=float(x)

def Return():
   if debug > 0:
      print(prompt+'Finished')
      press_return()
   return fH, nH, GammaxH, VxH, pH, TH, qxH, qxH_total, NetHSource, Sion, QH, RxH, QH_total, AlbedoH, WallH, error