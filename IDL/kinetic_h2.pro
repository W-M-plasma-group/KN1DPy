;+
; Kinetic_H2.pro
;
; This subroutine is part of the "KN1D" atomic and molecular neutal transport code.
;
;   This subroutine solves a 1-D spatial, 2-D velocity kinetic neutral transport
; problem for molecular hydrogen or deuterium (H2) by computing successive generations of
; charge exchange and elastic scattered neutrals. The routine handles electron-impact
; ionization and dissociation, molecular ion charge exchange, and elastic
; collisions with hydrogenic ions, neutral atoms, and molecules.
;
;   The positive vx half of the atomic neutral distribution function is inputted at x(0)
; (with arbitrary normalization). The desired flux on molecules entering the slab geometry at
; x(0) is specified. Background profiles of plasma ions, (e.g., Ti(x), Te(x), n(x), vxi(x),...)
; and atomic distribution function (fH) is inputted. (fH can be computed by procedure 
; "Kinetic_H.pro".) Optionally, a molecular source profile (SH2(x)) is also inputted. 
; The code returns the molecular hydrogen distribution function, fH2(vr,vx,x) for all 
; vx, vr, and x of the specified vr,vx,x grid. The atomic (H) and ionic (P) hydrogen 
; source profiles and the atomic source velocity distribution functions 
; resulting from Franck-Condon reaction product energies of H are also returned.
;
;   Since the problem involves only the x spatial dimension, all distribution functions
; are assumed to have rotational symmetry about the vx axis. Consequently, the distributions
; only depend on x, vx and vr where vr =sqrt(vy^2+vz^2)
;
;  History:
;
;    B. LaBombard   First coding based on Kinetic_Neutrals.pro 22-Dec-2000
;
;    For more information, see write-up: "A 1-D Space, 2-D Velocity, Kinetic
;    Neutral Transport Algorithm for Hydrogen Molecules in an Ionizing Plasma", B. LaBombard
;
; Note: Variable names contain characters to help designate species -
;       atomic neutral (H), molecular neutral (H2), molecular ion (HP), proton (i) or (P)
;
;________________________________________________________________________________
pro Kinetic_H2,vx,vr,x,Tnorm,mu,Ti,Te,n,vxi,fH2BC,GammaxH2BC,NuLoss,PipeDia,fH,SH2,$
       fH2,nH2,GammaxH2,VxH2,pH2,TH2,qxH2,qxH2_total,Sloss,QH2,RxH2,QH2_total,AlbedoH2,WallH2,$
       truncate=truncate,$
       nHP,THP,fSH,SH,SP,SHP,NuE,NuDis,$
       Simple_CX=Simple_CX,Max_Gen=Max_Gen,Compute_H_Source=Compute_H_Source,$
       No_Sawada=No_Sawada,H2_H2_EL=H2_H2_EL,H2_P_EL=H2_P_EL,H2_H_EL=_H2_H_EL,H2_HP_CX=H2_HP_CX,ni_correct=ni_correct,$
       ESH=ESH,Eaxis=Eaxis,error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause

   common Kinetic_H2_Output,piH2_xx,piH2_yy,piH2_zz,RxH2CX,RxH_H2,RxP_H2,RxW_H2,EH2CX,EH_H2,EP_H2,EW_H2,Epara_PerpH2_H2

   common Kinetic_H2_Errors,Max_dx,vbar_error,mesh_error,moment_error,C_Error,CX_Error,Swall_Error,H2_H2_error,Source_Error,$
                            qxH2_total_error,QH2_total_error
;
;  Input:
;		  vx(*)	- fltarr(nvx), normalized x velocity coordinate 
;			  [negative values, positive values],
;			  monotonically increasing. Note: a nonuniform mesh can be used.
;			  Dimensional velocity (note: Vth is based on ATOM mass)
;			  is v = Vth * vx where Vth=sqrt(2 k Tnorm/(mH*mu))
;			  Note: nvx must be even and vx(*) symmetric about 
;			  zero but not contain a zero element
;		  vr(*)	- fltarr(nvr), normalized radial velocity coordinate 
;			  [positive values], monotonically increasing. Note: a non-uniform mesh can be used.
;			  Dimensional velocity is v = Vth * vr where Vth=sqrt(2 k Tnorm/(mH*mu)) 
;			  Note: vr must not contain a zero element
;		   x(*)	- fltarr(nx), spatial coordinate (meters), 
;			  positive, monontonically increasing. Note: a non-uniform mesh can be used.
;		  Tnorm	- Float, temperature corresponding to the thermal speed (see vx and vr above) (eV)
;		     mu	- Float, 1=hydrogen, 2=deuterium
;		     Ti	- fltarr(nx), Ion temperature profile (eV)
;		     Te	- fltarr(nx), electron temperature profile (eV)
;		      n	- fltarr(nx), electron density profile (m^-3)
;		    vxi	- fltarr(nx), x-directed plasma ion and molecular ion flow profile (m s^-1)
;		  fH2BC	- fltarr(nvr,nvx), this is an input boundary condition
;			  specifying the shape of the neutral molecule velocity distribution 
;			  function at location x(0). Normalization is arbitrary.
;		          Only values with positive vx, fH2BC(*,nvx/2:*) are used
;		          by the code.
;	     GammaxH2BC	- float, desired neutral molecule flux density in the +Vx
;			  direction at location x(0) (m^-2 s^-1)
;			  fH2BC is scaled to yield this flux density.
;	         NuLoss	- fltarr(nx), characteristic loss frequency for HP molecules (1/s)
;			  (for open field lines, this is ~Cs/L). If this variable is undefined,
;			  then NuLoss set set to zero.
;	        PipeDia	- fltarr(nx), effective pipe diameter (meters)
;			  This variable allows collisions with the 'side-walls' to be simulated.
;			  If this variable is undefined, then PipeDia set set to zero. Zero values
;			  of PipeDia are ignored (i.e., treated as an infinite diameter).
;                    fH	- fltarr(nvr,nvx,nx), neutral atomic velocity distribution
;                         function. fH is normalized so that the atomic neutral density, nH(k), is 
;			  defined as the velocity space integration:
;                           nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
;			  If this variable is undefined, then no molecule-atom collisions are included.
;			  NOTE: dVx is velocity space differential for Vx axis and Vr2pidVr = Vr*!pi*dVr
;		                with dVr being velocity space differential for Vr axis.
;		    SH2	- fltarr(nx), source profile of wall-temperature (room temperature) H2 molecules (m^-3 s^-1)
;			  If this variable is undefined, it is set equal to zero.
;
;  Input & Output:
;                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
;                         function.
;			  'Seed' values for this may be specified on input. If this parameter
;			  is undefined, then a zero 'seed' value will be used. 
;			  The algorithm outputs a self-consistent fH2.
;			  fH2 is normalized so that the neutral density, nH2(k), is defined as 
;			  the velocity space integration:
;                             nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
;		    nHP	- fltarr(nx), molecular ion density profile (m^-3)
;			     'Seed' values for this may be specified on input. If this parameter
;			     is undefined, then a zero 'seed' value will be used. 
;			     The algorithm outputs a self-consistent profile for nHP.
;		    THP	- fltarr(nx), molecular ion temperature profile (m^-3)
;			     'Seed' values for this may be specified on input. If this parameter
;			     is undefined, then a 'seed' value of 3 eV will be used. 
;			     The algorithm outputs a self-consistent profile for THP.
;
;  Output:
;                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
;                         function. fH2 is normalized so that the neutral density, nH2(k), is 
;			  defined as the velocity space integration:
;                             nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
;                   nH2	- fltarr(nx), neutral molecule density profile (m^-3)
;              GammaxH2	- fltarr(nx), neutral flux profile (# m^-2 s^-1)
;                             computed from GammaxH2(k)=Vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))
;                  VxH2	- fltarr(nx), neutral velocity profile (m s^-1)
;                             computed from GammaxH2/nH2
;
;                       To aid in computing the some of the quantities below, the procedure internally
;                       defines the quantities:
;                       vr2vx2_ran(i,j,k)=vr(i)^2+(vx(j)-VxH2(k))^2
;                                     which is the magnitude of 'random v^2' at each mesh point
;                       vr2vx2(i,j,k)=vr(i)^2+vx(j)^2
;                                     which is the magnitude of 'total v^2' at each mesh point
;                       q=1.602177D-19, mH=1.6726231D-27
;                       C(*,*,*) is the right hand side of the Boltzmann equation, evaluated
;                                using the computed neutral distribution function
;
;                   pH2	- fltarr(nx), neutral pressure (eV m^-2) computed from:
;                         pH2(k)~vth2*total(Vr2pidVr*(vr2vx2_ran(*,*,k)*fH2(*,*,k))#dVx))*(mu*mH)/(3*q)
;                   TH2	- fltarr(nx), neutral temperature profile (eV) computed from: TH2=pH2/nH2
;                  qxH2	- fltarr(nx), neutral random heat flux profile (watts m^-2) computed from:
;                          qxH2(k)~vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH2(*,*,k))#(dVx*(vx-VxH2(k)))))*0.5*(mu*mH)
;            qxH2_total	- fltarr(nx), total neutral heat flux profile (watts m^-2)
;                             This is the total heat flux transported by the neutrals:
;                         qxH2_total=(0.5*nH2*(mu*mH)*VxH2*VxH2 + 2.5*pH2*q)*VxH2 + piH2_xx*VxH2 + qxH2
;                 Sloss	- fltarr(nx), H2 loss rate from ionization and dissociation (SH2) (m^-3 s^-1)
;                   QH2	- fltarr(nx), rate of net thermal energy transfer into neutral molecules (watts m^-3) computed from
;                               QH2(k)~vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C(*,*,k))#dVx))*0.5*(mu*mH)
;                  RxH2	- fltarr(nx), rate of x momentum transfer to neutral molecules (=force, N m^-2).
;                               RxH2(k)~Vth*total(Vr2pidVr*(C(*,*,k)#(dVx*(vx-VxH2(k)))))*(mu*mH)
;             QH2_total	- fltarr(nx), net rate of total energy transfer into neutral molecules
;                          = QH2 + RxH2*VxH2 - 0.5*(mu*mH)*(Sloss-SH2)*VxH2*VxH2 (watts m^-3)
;              AlbedoH2	- float, Ratio of molceular particle flux with Vx < 0 divided by particle flux
;                          with Vx > 0  at x=x(0)
;                          (Note: For SH2 non-zero, the flux with Vx < 0 will include
;                          contributions from molecular hydrogen sources within the 'slab'.
;                          In this case, this parameter does not return the true 'Albedo'.)
;	         WallH2 - fltarr(nx), molecular sink rate and source rate from interation with 'side walls' (m^-3 s^-1)
;			  Note: There is no loss of molecules when they strike the side walls in this model. The 
;		          molecules are just relaunched with a Maxwellian distribution at the wall temperature.
;
;	            nHP	- fltarr(nx), molecular ion density profile (m^-3)
;			     'Seed' values for this may be specified on input.
;		    THP	- fltarr(nx), molecular ion temperature profile (m^-3)
;			     'Seed' values for this may be specified on input.
;		    fSH	- fltarr(nvr,nvx,nx), H source velocity distribution.
;			  fSH is normalized so that the total atomic neutral 
;			  source, SH(k), is defined as the velocity space integration:
;			      SH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
;		     SH	- fltarr(nx), atomic neutral source profile (m^-3 s^-1)
;  		             computed from total(Vr2pidVr*(fSH(*,*,k)#dVx))
;		     SP	- fltarr(nx), proton source profile (m^-3 s^-1)
;		    SHP	- fltarr(nx), molecular ion source profile (m^-3 s^-1)
;		    NuE	- fltarr(nx), energy equilibration frequency of molecular ion with bulk plasma (1/s)
;		  NuDis	- fltarr(nx), molecular ion dissociation frequency (1/s)
;
; KEYWORDS:
;   Output:
;		    ESH	- fltarr(nvr,nx) - returns normalized H source energy distribution
;			     derived from velocity-space distribution:
;  		            ESH(*,*) = 0.5*mu*mH*Vr^2 FSH(*,vx=0,*)*4*pi*Vr/(mu*mH)
;		  Eaxis	- fltarr(nvr) - returns energy coordinate for ESH (eV) = 0.5*mu*mH*Vr^2
;	          error	- Returns error status: 0=no error, solution returned
;					        1=error, no solution returned
;
; COMMON BLOCK KINETIC_H2_OUTPUT
;    Output:
;               piH2_xx	- fltarr(nx), xx element of stress tensor (eV m^-2) computed from:
;                         piH2_xx(k)~vth2*total(Vr2pidVr*(fH2(*,*,k)#(dVx*(vx-VxH2(k))^2)))*(mu*mH)/q - pH2
;               piH2_yy	- fltarr(nx), yy element of stress tensor (eV m^-2) computed from:
;                         piH2_yy(k)~vth2*total((Vr2pidVr*Vr^2)*(fH2(*,*,k)#dVx))*(mu*mH)/q - pH2
;               piH2_zz	- fltarr(nx), zz element of stress tensor (eV m^-2) = piH2_yy
;			   Note: cylindrical system relates r^2 = y^2 + z^2. All other stress tensor elements are zero.
;
;           The following momentum and energy transfer rates are computed from charge-exchange collsions between species:
;                RxH2CX	- fltarr(nx), rate of x momentum transfer from molecular ions to neutral molecules (=force/vol, N m^-3).
;                 EH2CX	- fltarr(nx), rate of energy transfer from molecular ions to neutral molecules (watts m^-3).
;
;           The following momentum and energy transfer rates are computed from elastic collsions between species:
;                RxH_H2	- fltarr(nx), rate of x momentum transfer from neutral atoms to molecules (=force/vol, N m^-3).
;                RxP_H2	- fltarr(nx), rate of x momentum transfer from hydrogen ions to neutral molecules (=force/vol, N m^-3).
;                 EH_H2	- fltarr(nx), rate of energy transfer from neutral atoms to molecules (watts m^-3).
;                 EP_H2	- fltarr(nx), rate of energy transfer from hydrogen ions to neutral molecules (watts m^-3).
;
;           The following momentum and energy transfer rates are computed from collisions with the 'side-walls'
;                RxW_H2	- fltarr(nx), rate of x momentum transfer from wall to molecules (=force/vol, N m^-3).
;                 EW_H2	- fltarr(nx), rate of energy transfer from wall to neutral molecules (watts m^-3).
;
;           The following is the rate of parallel to perpendicular energy transfer computed from elastic collisions
;       Epara_PerpH2_H2	- fltarr(nx), rate of parallel to perp energy transfer within molecular hydrogen species (watts m^-3).
;
; KEYWORDS:
;   Input:
;	       truncate	- float, stop computation when the maximum 
;			  increment of neutral density normalized to 
;			  inputed neutral density is less than this 
;	    		  value in a subsequent generation. Default value is 1.0e-4
;
;             Simple_CX	- if set, then use CX source option (B): Neutral molecules are born
;                         in velocity with a distribution proportional to the local
;                         molecular ion distribution function. Simple_CX=1 is default.
;
;                         if not set, then use CX source option (A): The CX source
;                         neutral molecule distribution function is computed by evaluating the
;                         the CX cross section for each combination of (vr,vx,vr',vx')
;                         and convolving it with the molecule neutral distribution function.
;                         This option requires more CPU time and memory.
;
;      	        Max_gen	- integer, maximum number of collision generations to try including before giving up.
;                         Default is 50.
;      Compute_H_Source	-  if set, then compute fSH, SH, SP, and SHP
;
;	      No_Sawada	- if set, then DO NOT correct reaction rates according to
;			     results from collisional-radiative model of Sawada
;			     [Sawada, K. and Fujimoto, T., Journal of Applied Physics 78 (1995) 2913.]
;	       H2_H2_EL	- if set, then include H2 -> H2 elastic self collisions
;			     Note: if H2_H2_EL is set, then algorithm iterates fH2 until
;	                     self consistent fH2 is achieved.
;	       H2_HP_CX	- if set, then include H2 -> H2(+) charge exchange collisions
;			     Note: if H2_HP_CX is set, then algorithm iterates until
;	                     self consistent nHp is achieved.
;	        H2_P_EL	- if set, then include H2 -> H(+) elastic collisions (does not lead to an iteration loop)
;	        H2_H_EL	- if set, then include H2 -> H elastic collisions (does not lead to an iteration loop)
;            ni_correct	- if set, then algorithm corrects hydrogen ion density
;			     according to quasineutrality: ni=ne-nHp
;			     This is done in an iteration loop.
;
;	 Compute_Errors	- if set, then return error estimates in common block KINETIC_H2_ERRORS below
;
;		   plot	- 0= no plots, 1=summary plots, 2=detail plots, 3=very detailed plots
;		  debug	- 0= do not execute debug code, 1=summary debug, 2=detail debug, 3=very detailed debug
;	        debrief	- 0= do not print, 1=print summary information, 2=print detailed information
;	          pause	- if set, then pause between plots
;
; COMMON BLOCK KINETIC_H2_ERRORS
;
;	if COMPUTE_ERRORS keyword is set then the following is returned in common block KINETIC_H2_ERRORS
;
;	         Max_dx	- float(nx), Max_dx(k) for k=0:nx-2 returns maximum 
;			  allowed x(k+1)-x(k) that avoids unphysical negative 
;			  contributions to fH2
;	     Vbar_error	- float(nx), returns numerical error in computing
;			  the speed of neutrals averged over maxwellian distribution;
;			  over a temperature range spanning the Franck Condon energies
;		          of reactions R2, R3, R4, R5, R6, R7, R8, R10
;			  The average speed should be:
;				 vbar_exact=2*Vth*sqrt(Ti(*)/Tnorm)/sqrt(!pi)
;			  Vbar_error returns: abs(vbar-vbar_exact)/vbar_exact
;			  where vbar is the numerically computed value.
;	     mesh_error	- fltarr(nvr,nvx,nx), normalized error of solution
;			  based on substitution into Boltzmann equation.
;	   moment_error	- fltarr(nx,m), normalized error of solution
;			  based on substitution into velocity space
;			  moments (v^m) of Boltzmann equation, m=[0,1,2,3,4]
;      	        C_error	- fltarr(nx), normalized error in charge exchange and elastic scattering collision 
;				      operator. This is a measure of how well the charge exchange and
;				      elastic scattering portions of the collision operator
;				      conserve particles.
;      	       CX_error	- fltarr(nx), normalized particle conservation error in charge exchange collision operator.
;	    H2_H2_error	- fltarr(nx,[0,1,2]) return normalized errors associated with 
;		          particle [0], x-momentum [1], and total energy [2] convervation of the elastic self-collision operator
;
;      	   Source_error	- fltarr(nx), normalized error in mass continuity equation. This is a measure
;				      of how well atomic plus molecular mass continuity is satisfied.
;      qxH2_total_error	- fltarr(nx), normalized error estimate in computation of qxH2_total
;       QH2_total_error	- fltarr(nx), normalized error estimate in computation of QH2_total
;
; History:
;	22-Dec-2000 - B. LaBombard - first coding.
;	06-Feb-2001 - B. LaBombard - added elastic collisions and molecular sources
;
;______________________________________________________________________
;-
   prompt='Kinetic_H2 => '
;
   common Kinetic_H2_input,vx_s,vr_s,x_s,Tnorm_s,mu_s,Ti_s,Te_s,n_s,vxi_s,fH2BC_s,GammaxH2BC_s,NuLoss_s,PipeDia_s,fH_s,SH2_s,fH2_s,$
               nHP_s,THP_s,Simple_CX_s,Sawada_s,H2_H2_EL_s,H2_P_EL_s,H2_H_EL_s,H2_HP_CX_s,ni_correct_s

   common Kinetic_H2_internal,vr2vx2,vr2vx_vxi2,fw_hat,fi_hat,fHp_hat,EH2_P,sigv,Alpha_Loss,v_v2,v_v,vr2_vx2,vx_vx,$
          Vr2pidVrdVx,SIG_CX,SIG_H2_H2,SIG_H2_H,SIG_H2_P,Alpha_CX,Alpha_H2_H,MH2_H2_sum,Delta_nH2s

   common Kinetic_H2_H_moments,nH,VxH,TH
;
; Internal Debug switches
;
   shifted_Maxwellian_debug=0
   CI_Test=1
   Do_Alpha_CX_Test=0
;
; Internal Tolerances
;
   DeltaVx_tol=.01
   Wpp_tol=.001
;
; Test input parameters
;
   key_default,truncate,1.0e-4
   key_default,Compute_Errors,0
   key_default,Compute_H_Source,0
   key_default,plot,0
   key_default,debug,0
   key_default,pause,0
   key_default,debrief,0
   if debug gt 0 then plot=plot > 1
   if debug gt 0 then debrief=debrief > 1
   if debug gt 0 then pause=1
   key_default,Simple_CX,1
   key_default,Max_Gen,50
   key_default,No_Sawada,0
   if No_Sawada eq 0 then Sawada=1 else Sawada=0
   key_default,H2_H2_EL,0
   key_default,H2_P_EL,0
   key_default,_H2_H_EL,0
   key_default,H2_HP_CX,0
   key_default,ni_correct,0
   error=0
   nvr=n_elements(vr)
   nvx=n_elements(vx)
   nx=n_elements(x)
   vr=double(vr)
   vx=double(vx)
   x=double(x)
   dx=x-shift(x,1) & dx=dx(1:*)
   notpos=where(dx le 0.0,count)
   if count gt 0 then begin & print,prompt+'x(*) must be increasing with index!' & error=1 & goto,return & endif
   if (nvx mod 2) ne 0 then begin & print,prompt+'Number of elements in vx must be even!' & error=1 & goto,return & endif
   if n_elements(Ti) ne nx then begin & print,prompt+'Number of elements in Ti and x do not agree!' & error=1 & goto,return & endif
   if type_of(vxi) eq 0 then vxi=dblarr(nx)
   ; vxi=dblarr(nx)
   if n_elements(vxi) ne nx then begin & print,prompt+'Number of elements in vxi and x do not agree!' & error=1 & goto,return & endif
   if n_elements(Te) ne nx then begin & print,prompt+'Number of elements in Te and x do not agree!' & error=1 & goto,return & endif
   if n_elements(n) ne nx then begin & print,prompt+'Number of elements in n and x do not agree!' & error=1 & goto,return & endif
   if type_of(NuLoss) eq 0 then NuLoss=fltarr(nx)
   ; NuLoss=fltarr(nx)
   if n_elements(NuLoss) ne nx then begin & print,prompt+'Number of elements in NuLoss and x do not agree!' & error=1 & goto,return & endif
   if type_of(PipeDia) eq 0 then PipeDia=dblarr(nx)
   if n_elements(PipeDia) ne nx then begin & print,prompt+'Number of elements in PipeDia and x do not agree!' & error=1 & goto,return & endif
   if type_of(GammaxH2BC) eq 0 then  begin & print,prompt+'GammaxH2BC is not defined!' & error=1 & goto,return & endif
   if type_of(fH) eq 0 then fH=dblarr(nvr,nvx,nx)
   ; fH=dblarr(nvr,nvx,nx)
   if n_elements(fH(*,0,0)) ne nvr then begin & print,prompt+'Number of elements in fH(*,0,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fH(0,*,0)) ne nvx then begin & print,prompt+'Number of elements in fH(0,*,0) and vx do not agree!' & error=1 & goto,return & endif
   if n_elements(fH(0,0,*)) ne nx then begin & print,prompt+'Number of elements in fH(0,0,*) and x do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2BC(*,0)) ne nvr then begin & print,prompt+'Number of elements in fH2BC(*,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2BC(0,*)) ne nvx then begin & print,prompt+'Number of elements in fH2BC(0,*) and vx do not agree!' & error=1 & goto,return & endif
   if type_of(fH2) eq 0 then fH2=dblarr(nvr,nvx,nx)
   ; fH2=dblarr(nvr,nvx,nx)
   if n_elements(fH2(*,0,0)) ne nvr then begin & print,prompt+'Number of elements in fH2(*,0,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2(0,*,0)) ne nvx then begin & print,prompt+'Number of elements in fH2(0,*,0) and vx do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2(0,0,*)) ne nx then begin & print,prompt+'Number of elements in fH2(0,0,*) and x do not agree!' & error=1 & goto,return & endif
   if type_of(SH2) eq 0 then SH2=dblarr(nx)
   if n_elements(SH2) ne nx then begin & print,prompt+'Number of elements in SH2 and x do not agree!' & error=1 & goto,return & endif
   if type_of(nHP) eq 0 then nHP=dblarr(nx)
   ; nHP=dblarr(nx)
   if n_elements(nHP) ne nx then begin & print,prompt+'Number of elements in nHP and x do not agree!' & error=1 & goto,return & endif
   if type_of(THP) eq 0 then THP=dblarr(nx)+3.0
   ; THP=dblarr(nx)+3.0
   if n_elements(THP) ne nx then begin & print,prompt+'Number of elements in THP and x do not agree!' & error=1 & goto,return & endif
   if total(abs(vr)) le 0.0 then begin & print,prompt+'vr is all 0!' & error=1  & goto,return & endif
   ii=where(vr le 0,count)
   if count gt 0 then begin & print,prompt+'vr contains zero or negative element(s)!' & error=1 & goto,return & endif
   if total(abs(vx)) le 0.0 then begin & print,prompt+'vx is all 0!' & error=1 & goto,return & endif
   if total(x) le 0.0 then begin & print,prompt+'Total(x) is less than or equal to 0!' & error=1 & goto,return & endif
   if type_of(Tnorm) eq 0 then begin & print,prompt+'Tnorm is not defined!' & error=1 & goto,return & endif
   if type_of(mu) eq 0 then begin & print,prompt+'mu is not defined!' & error=1 & goto,return & endif
   if mu ne 1 and mu ne 2 then begin & print,prompt+'mu must be 1 or 2!' & error=1 & goto,return & endif

   
   _e='e!U-!N'
   if mu eq 1 then begin
      _p='H!U+!N'
      _H='H!U0!N'
      _H1s='H(1s)'
      _H2s='H!U*!N(2s)'
      _H2p='H!U*!N(2p)'
      _Hn2='H!U*!N(n=2)'
      _Hn3='H!U*!N(n=3)'
      _Hn='H!U*!N(n>=2)'
      _HH='H!D2!N'
      _Hp='H!D2!U+!N'
   endif else begin
      _p='D!U+!N'
      _H='D!U0!N'
      _H1s='D(1s)'
      _H2s='D!U*!N(2s)'
      _H2p='D!U*!N(2p)'
      _Hn2='D!U*!N(n=2)'
      _Hn3='D!U*!N(n=3)'
      _Hn='D!U*!N(n>=2)'
      _HH='D!D2!N'
      _Hp='D!D2!U+!N'
   endelse
   plus=' + '
   arrow=' -> '
   elastic=' (elastic)'
   _R1=_e+plus+_HH+arrow+_e+plus+_Hp+plus+_e
   _R2=_e+plus+_HH+arrow+_e+plus+_H1s+plus+_H1s
   _R3=_e+plus+_HH+arrow+_e+plus+_H1s+plus+_H2s
   _R4=_e+plus+_HH+arrow+_e+plus+_p+plus+_H1s+plus+_e
   _R5=_e+plus+_HH+arrow+_e+plus+_H2p+plus+_H2s
   _R6=_e+plus+_HH+arrow+_e+plus+_H1s+plus+_Hn3
   _R7=_e+plus+_Hp+arrow+_e+plus+_p+plus+_H1s
   _R8=_e+plus+_Hp+arrow+_e+plus+_p+plus+_Hn2
   _R9=_e+plus+_Hp+arrow+_e+plus+_p+plus+_p+plus+_e
   _R10=_e+plus+_Hp+arrow+_H1s+plus+_Hn
   _R11=_HH+plus+_p+arrow+_HH+plus+_p+elastic
   _R12=_HH+plus+_H+arrow+_HH+plus+_H+elastic
   _R13=_HH+plus+_HH+arrow+_HH+plus+_HH+elastic
   _R14=_HH+plus+_Hp+arrow+_Hp+plus+_HH
   _Rn=[' ',_R1,_R2,_R3,_R4,_R5,_R6,_R7,_R8,_R9,_R10,_R11,_R12,_R13,_R14]
   in=where(vx lt 0,count)
   if count lt 1 then begin & print,prompt+'vx contains no negative elements!' & error=1 & goto,return & endif
   ip=where(vx gt 0,count)
   if count lt 1 then begin & print,prompt+'vx contains no positive elements!' & error=1 & goto,return & endif
   iz=where(vx eq 0,count)
   if count gt 0 then begin & print,prompt+'vx contains one or more zero elements!' & error=1 & goto,return & endif
   diff=where(vx(ip) ne -reverse(vx(in)),count)
   if count gt 0 then begin & print,prompt+'vx array elements are not symmetric about zero!' & error=1 & goto,return & endif
   fH2BC_input=fH2BC
   fH2BC_input(*)=0.0
   fH2BC_input(*,ip)=fH2BC(*,ip)
   test=total(fH2BC_input)
   if test le 0.0 then begin & print,prompt+'Values for fH2BC(*,*) with vx > 0 are all zero!' & error=1 & goto,return & endif
;
; Output variables
;
   nH2=dblarr(nx)
   GammaxH2=dblarr(nx)
   VxH2=dblarr(nx)
   pH2=dblarr(nx)
   TH2=dblarr(nx)
   NuDis=dblarr(nx)
   NuE=dblarr(nx)

   qxH2=dblarr(nx)
   qxH2_total=dblarr(nx)
   Sloss=dblarr(nx)
   WallH2=dblarr(nx)
   QH2=dblarr(nx)
   RxH2=dblarr(nx)
   QH2_total=dblarr(nx)
   piH2_xx=dblarr(nx)
   piH2_yy=dblarr(nx)
   piH2_zz=dblarr(nx)
   RxH2CX=dblarr(nx)
   RxH_H2=dblarr(nx)
   RxP_H2=dblarr(nx)
   RxW_H2=dblarr(nx)
   EH2CX=dblarr(nx)
   EH_H2=dblarr(nx)
   EP_H2=dblarr(nx)
   EW_H2=dblarr(nx)
   Epara_PerpH2_H2=dblarr(nx)
   AlbedoH2=0.0D0

   fSH=dblarr(nvr,nvx,nx)
   SH=dblarr(nx)
   SP=dblarr(nx)
   SHP=dblarr(nx)
   ESH=dblarr(nvr,nx)
   Eaxis=dblarr(nx)
;
; Internal variables
;
   mH=1.6726231D-27
   q=1.602177D-19				
   k_boltz=1.380658D-23                 &;Bolzmann's constant, J K^-1
   Twall=293.0*k_boltz/q                &;room temperature (eV)
   Work=dblarr(nvr*nvx)
   fH2G=dblarr(nvr,nvx,nx)
   NH2G=dblarr(nx,max_gen+1)
   Vth=sqrt(2*q*Tnorm/(mu*mH))
   Vth2=vth*vth
   Vth3=Vth2*Vth
   fH2s=dblarr(nx)
   nH2s=dblarr(nx)
   THPs=dblarr(nx)
   nHPs=dblarr(nx)
   Alpha_H2_H2=dblarr(nvr,nvx)
   Omega_H2_P=dblarr(nx)
   Omega_H2_H=dblarr(nx)
   Omega_H2_H2=dblarr(nx)
   VxH2G=dblarr(nx)
   TH2G=dblarr(nx)
   Wperp_paraH2=dblarr(nx)
   vr2vx2_ran2=dblarr(nvr,nvx)
   vr2_2vx_ran2=dblarr(nvr,nvx)
   vr2_2vx2_2D=dblarr(nvr,nvx)
   RxCI_CX=dblarr(nx)
   RxCI_H_H2=dblarr(nx)
   RxCI_P_H2=dblarr(nx)
   Epara_Perp_CI=dblarr(nx)
   CI_CX_error=fltarr(nx)
   CI_H_H2_error=fltarr(nx)
   CI_P_H2_error=fltarr(nx)
   CI_H2_H2_error=fltarr(nx)
   Maxwell=dblarr(nvr,nvx,nx)

   Make_dVr_dVx,vr,vx,Vr2pidVr,VrVr4pidVr,dVx,vrL=vrL,vrR=vrR,vxL=vxL,vxR=vxR,$
                Vol=Vol,Vth_DeltaVx=Vth_DVx,Vx_DeltaVx=Vx_DVx,Vr_DeltaVr=Vr_DVr,Vr2Vx2=Vr2Vx2_2D,$
                jpa=jpa,jpb=jpb,jna=jna,jnb=jnb
;
; Vr^2-2*Vx^2
;
   for i=0,nvr-1 do vr2_2vx2_2D(i,*)=vr(i)^2-2*vx^2
;
; Theta-prime coordinate
;
   ntheta=5     &; use 5 theta mesh points for theta integration
   dTheta=replicate(1.0d0,ntheta)/ntheta
   theta=!pi*(dindgen(ntheta)/ntheta+0.5/ntheta)
   cos_theta=cos(theta)
;
; Determine energy space differentials
;
   Eaxis=vth2*0.5*mu*mH*vr^2/q
   _Eaxis=[Eaxis,2*Eaxis(nvr-1)-Eaxis(nvr-2)]
   Eaxis_mid=[0.0,0.5*(_Eaxis+shift(_Eaxis,-1))]
   dEaxis=shift(Eaxis_mid,-1)-Eaxis_mid
   dEaxis=dEaxis(0:nvr-1)
;
; Scale input molecular distribution function to agree with desired flux
;
   gamma_input=1.0
   if abs(GammaxH2BC) gt 0 then gamma_input=vth*total(Vr2pidVr*(fH2BC_input#(Vx*dVx)))
   ratio=abs(GammaxH2BC)/gamma_input
   fH2BC_input=fH2BC_input*ratio
   if abs(ratio-1) gt 0.01*truncate then begin
      fH2BC=fH2BC_input
   endif
   fH2(*,ip,0)=fH2BC_input(*,ip)
;
; if fH is zero, then turn off elastic H2 <-> H collisions
;
   H2_H_EL=_H2_H_EL
   if total(fH) le 0.0 then H2_H_EL=0
;
; Set iteration scheme
;
   fH2_iterate=0
   if (H2_H2_EL ne 0) or (H2_HP_CX ne 0) or (H2_H_EL ne 0) or (H2_P_EL ne 0) or (ni_correct ne 0) then fH2_iterate=1

   fH2_generations=0
   if (fH2_iterate ne 0) then fH2_generations=1
;
; Set flags to make use of previously computed local parameters 
;
   New_Grid=1
   if type_of(vx_s) ne 0 then begin
      test=0
      ii=where(vx_s ne vx,count) & test=test+count
      ii=where(vr_s ne vr,count) & test=test+count
      ii=where(x_s ne x,count) & test=test+count
      ii=where(Tnorm_s ne Tnorm,count) & test=test+count
      ii=where(mu_s ne mu,count) & test=test+count
      if test le 0 then New_Grid=0
   endif
   New_Protons=1
   if type_of(Ti_s) ne 0 then begin
      test=0
      ii=where(Ti_s ne Ti,count) & test=test+count
      ii=where(n_s ne n,count) & test=test+count
      ii=where(vxi_s ne vxi,count) & test=test+count
      if test le 0 then New_Protons=0
   endif
   New_Electrons=1
   if type_of(Te_s) ne 0 then begin
      test=0
      ii=where(Te_s ne Te,count) & test=test+count
      ii=where(n_s ne n,count) & test=test+count
      if test le 0 then New_Electrons=0
   endif
   New_fH=1
   if type_of(fH_s) ne 0 then begin
      ii=where(fH_s ne fH,count)
      if count le 0 then New_fH=0
   endif
   New_Simple_CX=1
   if type_of(Simple_CX_s) ne 0 then begin
      ii=where(Simple_CX_s ne Simple_CX,count)
      if count le 0 then New_Simple_CX=0
   endif
   New_H2_Seed=1
   if type_of(fH2_s) ne 0 then begin
      ii=where(fH2_s ne fH2,count)
      if count le 0 then New_H2_Seed=0
   endif
   New_HP_Seed=1
   if type_of(nHP_s) ne 0 then begin
      test=0
      ii=where(nHP_s ne nHP,count) & test=test+count
      ii=where(THP_s ne THP,count) & test=test+count
      if test le 0 then New_HP_Seed=0
   endif
   New_ni_correct=1
   if type_of(ni_correct_s) ne 0 then begin
      ii=where(ni_correct_s ne ni_correct,count)
      if count le 0 then New_ni_correct=0
   endif

   Do_sigv=        New_Grid or New_Electrons
   Do_fH_moments= (New_Grid or New_fH) and total(fH) gt 0.0
   Do_Alpha_CX=   (New_Grid or (type_of(Alpha_CX) eq 0) or New_HP_Seed or New_Simple_CX) and H2_HP_CX
;Do_Alpha_CX is updated in fH2_iteration loop
   Do_SIG_CX=     (New_Grid or (type_of(SIG_CX) eq 0) or New_Simple_CX) and (Simple_CX eq 0) and Do_Alpha_CX
   Do_Alpha_H2_H= (New_Grid or (type_of(Alpha_H2_H) eq 0) or New_fH)    and H2_H_EL
   Do_SIG_H2_H=   (New_Grid or (type_of(SIG_H2_H) eq 0))                and Do_Alpha_H2_H
   Do_SIG_H2_H2=  (New_Grid or (type_of(SIG_H2_H2) eq 0))               and H2_H2_EL
   Do_Alpha_H2_P= (New_Grid or (type_of(Alpha_H2_P) eq 0) or New_Protons or New_ni_correct) and H2_P_EL 
; Do_Alpha_H2_P is updated in fH2_iteration loop
   Do_SIG_H2_P=   (New_Grid or (type_of(SIG_H2_P) eq 0))                and Do_Alpha_H2_P
   Do_v_v2=      (New_Grid or (type_of(v_v2) eq 0))                   and (CI_Test or Do_SIG_CX or Do_SIG_H2_H or Do_SIG_H2_H2 or Do_SIG_H2_P)

   nH=dblarr(nx)
   VxH=dblarr(nx)
   TH=dblarr(nx)+1.0
   if Do_fH_moments then begin
      if debrief gt 1 then print,prompt+'Computing vx and T moments of fH'
;
; Compute x flow velocity and temperature of atomic species
;
      for k=0,nx-1 do begin
         nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
         if nH(k) gt 0 then begin
            VxH(k)=vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))/nH(k)
            for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-VxH(k)/vth)^2
            TH(k)=(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fH(*,*,k))#dVx))/(3*q*nH(k))
         endif
      endfor
   endif

   if New_Grid then begin
      if debrief gt 1 then print,prompt+'Computing vr2vx2, vr2vx_vxi2, EH2_P'
;
; Magnitude of total normalized v^2 at each mesh point
;
      vr2vx2=dblarr(nvr,nvx,nx)
      for i=0,nvr-1 do for k=0,nx-1 do vr2vx2(i,*,k)=vr(i)^2+vx^2
;
; Magnitude of total normalized (v-vxi)^2 at each mesh point
;
      vr2vx_vxi2=dblarr(nvr,nvx,nx)
      for i=0,nvr-1 do for k=0,nx-1 do vr2vx_vxi2(i,*,k)=vr(i)^2+(vx-vxi(k)/vth)^2
;
; Molecular hydrogen ion energy in local rest frame of plasma at each mesh point
;
      EH2_P=mH*vr2vx_vxi2*vth2/q
      EH2_P=EH2_P > 0.1    &; sigmav_cx does not handle neutral energies below 0.1 eV
      EH2_P=EH2_P < 2.0E4  &; sigmav_cx does not handle neutral energies above 20 keV
;
; Compute Maxwellian H2 distribution at the wall temperature
;
      fw_Hat=dblarr(nvr,nvx)
;
; NOTE: Molecular ions have 'normalizing temperature' of 2 Tnorm, i.e., in order to
; achieve the same thermal velocity^2, a molecular ion distribution has to have twice the temperature 
; as an atomic ion distribution
;
      if (total(SH2) gt 0) or (total(PipeDia) gt 0) then begin
         if debrief gt 1 then print,prompt+'Computing fw_Hat'
         vx_shift=[0.0]
         Tmaxwell=[Twall]
         mol=2
         create_shifted_Maxwellian,vr,vx,Tmaxwell,vx_shift,mu,mol,Tnorm,_Maxwell
         fw_hat=_Maxwell(*,*,0)
      endif
   endif

   if New_Protons then begin
;
; Compute Fi_hat
;
      if debrief gt 1 then print,prompt+'Computing fi_Hat'
      vx_shift=vxi
      Tmaxwell=Ti
      mol=1
@create_shifted_maxwellian.include
      fi_hat=Maxwell
   endif

   if Do_sigv then begin
      if debrief gt 1 then print,prompt+'Computing sigv'
;
; Compute sigmav rates for each reaction and optionally apply
; CR model corrections of Sawada
;
      sigv=dblarr(nx,11)
;________________________________________________________________________________
; Reaction R1:  e + H2 -> e + H2(+) + e 
;________________________________________________________________________________
      sigv(*,1)=sigmav_ion_HH(Te)
      if Sawada then sigv(*,1)=sigv(*,1)*3.7/2.0
;________________________________________________________________________________
; Reaction R2:  e + H2 -> H(1s) + H(1s)  
;________________________________________________________________________________
      sigv(*,2)=sigmav_H1s_H1s_HH(Te)
      if Sawada then begin
;
;  Construct table
;
         Te_table=alog([5,20,100]) & Ne_Table=alog([1e14,1e17,1e18,1e19,1e20,1e21,1e22])
         fctr_Table=fltarr(7,3)
         fctr_Table(*,0)=[2.2, 2.2, 2.1, 1.9, 1.2,  1.1,  1.05]/5.3
         fctr_Table(*,1)=[5.1, 5.1, 4.3, 3.1, 1.5,  1.25, 1.25]/10.05
         fctr_Table(*,2)=[1.3, 1.3, 1.1, 0.8, 0.38, 0.24, 0.22]/2.1
         _Te=Te
         _Te=_Te > 5
         _Te=_Te < 100
         _n=n > 1e14
         _n=n < 1e22
         fctr=Path_Interp_2D(fctr_Table,Ne_Table,Te_table,alog(_n),alog(_Te))
         sigv(*,2)=(1.0+fctr)*sigv(*,2)
      endif
;________________________________________________________________________________
; Reaction R3:  e + H2 -> e + H(1s) + H*(2s)
;________________________________________________________________________________
      sigv(*,3)=sigmav_H1s_H2s_HH(Te)
;________________________________________________________________________________
; Reaction R4:  e + H2 -> e + p + H(1s)
;________________________________________________________________________________
      sigv(*,4)=sigmav_P_H1s_HH(Te)
      if Sawada then sigv(*,4)=sigv(*,4)*1.0/0.6
;________________________________________________________________________________
; Reaction R5:  e + H2 -> e + H*(2p) + H*(2s)
;________________________________________________________________________________
      sigv(*,5)=sigmav_H2p_H2s_HH(Te)
;________________________________________________________________________________
; Reaction R6:  e + H2 -> e + H(1s) + H*(n=3)
;________________________________________________________________________________
      sigv(*,6)=sigmav_H1s_Hn3_HH(Te)
;________________________________________________________________________________
; Reaction R7:  e + H2(+) -> e + p + H(1s)
;________________________________________________________________________________
      sigv(*,7)=sigmav_P_H1s_HP(Te)
;________________________________________________________________________________
; Reaction R8:  e + H2(+) -> e + p + H*(n=2)
;________________________________________________________________________________
      sigv(*,8)=sigmav_p_Hn2_HP(Te)
;________________________________________________________________________________
; Reaction R9:  e + H2(+) -> e + p + p + e
;________________________________________________________________________________
      sigv(*,9)=sigmav_P_P_HP(Te)
;________________________________________________________________________________
; Reaction R10:  e + H2(+) -> e + H(1s) + H*(n>=2)
;________________________________________________________________________________
      sigv(*,10)=sigmav_H1s_Hn_HP(Te)
;________________________________________________________________________________
;
; Total H2 destruction rate (normalized by vth) = sum of reactions 1-6
;
      alpha_loss=dblarr(nx)
      alpha_loss(*)=n*total(sigv(*,1:6),2)/vth
   endif
;
;________________________________________________________________________________
; Set up arrays for charge exchange and elastic collision computations, if needed
;________________________________________________________________________________
;
   if Do_v_v2 eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing v_v2, v_v, vr2_vx2, and vx_vx'
;
; v_v2=(v-v_prime)^2 at each double velocity space mesh point, including theta angle
;
      v_v2=dblarr(nvr,nvx,nvr,nvx,ntheta)
;
; vr2_vx2=(vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2
;         at each double velocity space mesh point, including theta angle
;
      vr2_vx2=dblarr(nvr,nvx,nvr,nvx,ntheta)
      for m=0,ntheta-1 do begin
         for l=0,nvx-1 do begin
            for k=0,nvr-1 do begin
               for i=0,nvr-1 do begin
                  v_v2(i,*,k,l,m)=Vr(i)^2+Vr(k)^2-2*Vr(i)*Vr(k)*cos_theta(m)+(Vx(*)-Vx(l))^2
                  vr2_vx2(i,*,k,l,m)=Vr(i)^2+Vr(k)^2-2*Vr(i)*Vr(k)*cos_theta(m)-2*(Vx(*)-Vx(l))^2
               endfor
            endfor
         endfor
      endfor
;
; v_v=|v-v_prime| at each double velocity space mesh point, including theta angle
;
      v_v=sqrt(v_v2)
;
; vx_vx=(vx-vx_prime) at each double velocity space mesh point
;
      vx_vx=dblarr(nvr,nvx,nvr,nvx)
      for j=0,nvx-1 do for l=0,nvx-1 do vx_vx(*,j,*,l)=vx(j)-vx(l)
;
; Set Vr'2pidVr'*dVx' for each double velocity space mesh point
;
      Vr2pidVrdVx=dblarr(nvr,nvx,nvr,nvx)
      for k=0,nvr-1 do Vr2pidVrdVx(*,*,k,*)=Vr2pidVr(k)
      for l=0,nvx-1 do Vr2pidVrdVx(*,*,*,l)=Vr2pidVrdVx(*,*,*,l)*dVx(l)
   endif
;
   if Simple_CX eq 0 and Do_SIG_CX eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_CX'
;________________________________________________________________________________
; Option (A) was selected: Compute SigmaV_CX from sigma directly.
; In preparation, compute SIG_CX for present velocity space grid, if it has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_cx * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=v_v*sigma_cx_hh(v_v2*(mH*vth2/q))
;
; Set SIG_CX = vr' x Integral{v_v*sigma_cx} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_CX=dblarr(nvr*nvx,nvr*nvx)
        SIG_CX(*)=Vr2pidVrdVx*(_Sig#dtheta)
;
; SIG_CX is now vr' * sigma_cx(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
   if Do_SIG_H2_H eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H2_H'
;________________________________________________________________________________
; Compute SIG_H2_H for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H2_H * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=v_v*Sigma_EL_H_HH(v_v2*(0.5*mH*vth2/q))
;
; NOTE: using H energy here for cross-sections tabulated as H->H2
;
; Set SIG_H2_H = vr' x vx_vx x Integral{v_v*sigma_H2_H} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_H2_H=dblarr(nvr*nvx,nvr*nvx)
        SIG_H2_H(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)
;
; SIG_H2_H is now vr' * vx_vx * sigma_H2_H(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
   if Do_SIG_H2_P eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H2_P'
;________________________________________________________________________________
; Compute SIG_H2_P for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H2_P * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=v_v*Sigma_EL_P_HH(v_v2*(0.5*mH*vth2/q))
;
; NOTE: using H energy here for cross-sections tabulated as P->H2
;
; Set SIG_H2_P = vr' x vx_vx x Integral{v_v*sigma_H2_P} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_H2_P=dblarr(nvr*nvx,nvr*nvx)
        SIG_H2_P(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)
;
; SIG_H2_P is now vr' * vx_vx * sigma_H2_P(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
   if Do_SIG_H2_H2 eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H2_H2'
;________________________________________________________________________________
; Compute SIG_H2_H2 for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H2_H2 * vr2_vx2 * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=vr2_vx2*v_v*Sigma_EL_HH_HH(v_v2*(mH*mu*vth2/q),/VIS)/8.0
;
; Note: For viscosity, the cross section for D -> D is the same function of
;       center of mass energy as H -> H.
;
; Set SIG_H2_H2 = vr' x Integral{vr2_vx2*v_v*sigma_H2_H2} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_H2_H2=dblarr(nvr*nvx,nvr*nvx)
        SIG_H2_H2(*)=Vr2pidVrdVx*(_Sig#dtheta)
;
; SIG_H2_H2 is now vr' * sigma_H2_H2(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
;________________________________________________________________________________ 
; Compute Alpha_H2_H for inputted fH, if it is needed and has not
; already been computed with the present input parameters
;________________________________________________________________________________ 
;
   if Do_Alpha_H2_H eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_H2_H'
      Alpha_H2_H=dblarr(nvr,nvx,nx)
      for k=0,nx-1 do begin
         Work(*)=fH(*,*,k)
         Alpha_H2_H(*,*,k)=SIG_H2_H#Work
      endfor
   endif
;
;________________________________________________________________________________
; Compute nH2
;________________________________________________________________________________
   for k=0,nx-1 do begin
      nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
   endfor
;
   if New_H2_Seed then begin
      MH2_H2_sum=dblarr(nvr,nvx,nx)
      Delta_nH2s=1.0
   endif
;
; Compute Side-Wall collision rate
;
   gamma_wall=dblarr(nvr,nvx,nx)
   for k=0,nx-1 do begin
      if PipeDia(k) gt 0.0 then begin
         for j=0,nvx-1 do gamma_wall(*,j,k)=2*vr/PipeDia(k)
      endif
   endfor
;

fH2_Iterate:
;
;  This is the iteration entry point for fH2, THP and nHP iteration.
;  Save 'seed' values for comparison later
;
   fH2s=fH2
   nH2s=nH2
   THPs=THP
   nHPs=nHP
;
;________________________________________________________________________________ 
; Compute Alpha_CX for present THP and nHP, if it is needed and has not
; already been computed with the present parameters
;________________________________________________________________________________ 
;
   if Do_Alpha_CX eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_CX'
;
; Set Maxwellian Molecular Ion Distribution Function (assumed to be drifting with ion velocity, vxi)
;
      vx_shift=vxi
      Tmaxwell=THP
      mol=2

@create_shifted_maxwellian.include
      fHp_hat=Maxwell
;
      if Simple_CX then begin
;________________________________________________________________________________
; Option (B): Use maxwellian weighted <sigma v>
;________________________________________________________________________________
;
; THp/mu at each mesh point
;
         THp_mu=dblarr(nvr,nvx,nx)
         for k=0,nx-1 do THp_mu(*,*,k)=THp(k)/mu
;
; Molecular Charge Exchange sink rate
;
         alpha_cx=sigmav_cx_HH(THp_mu,EH2_P)/vth
         for k=0,nx-1 do alpha_cx(*,*,k)=alpha_cx(*,*,k)*nHp(k)
;________________________________________________________________________________
;
      endif else begin
;________________________________________________________________________________
; Option (A): Compute SigmaV_CX from sigma directly via SIG_CX
;________________________________________________________________________________
;
         alpha_cx=dblarr(nvr,nvx,nx)
         for k=0,nx-1 do begin
            Work(*)=fHp_hat(*,*,k)*nHp(k)
            alpha_cx(*,*,k)=SIG_CX#Work
         endfor
         if do_alpha_cx_test then begin
            alpha_cx_test=sigmav_cx_HH(THp_mu,EH2_P)/vth
            for k=0,nx-1 do alpha_cx_test(*,*,k)=alpha_cx_test(*,*,k)*nHp(k)
            print,'Compare alpha_cx and alpha_cx_test'
            press_return
         endif
      endelse
   endif
;________________________________________________________________________________ 
; Compute Alpha_H2_P for present Ti and ni (optionally correcting for nHP), 
; if it is needed and has not already been computed with the present parameters
;________________________________________________________________________________ 
   if Do_Alpha_H2_P eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_H2_P'
      Alpha_H2_P=dblarr(nvr,nvx,nx)
      ni=n
      if ni_correct then ni=(n-nHp) > 0.0
      for k=0,nx-1 do begin
         Work(*)=fi_hat(*,*,k)*ni(k)
         Alpha_H2_P(*,*,k)=SIG_H2_P#Work
      endfor
   endif
;________________________________________________________________________________
; Compute Omega values if nH2 is non-zero
;________________________________________________________________________________
   ii=where(nH2 le 0,count)
   if count le 0 then begin
;
; Compute VxH2
;
      if H2_P_EL or H2_H_EL or H2_H2_EL then begin
         for k=0,nx-1 do begin
            VxH2(k)=vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))/nH2(k)
         endfor
      endif
      ; if H2_P_EL or H2_H_EL or H2_H2_EL then for k=0,nx-1 do VxH2(k)=vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))/nH2(k) - ZZZZZ used to be all in one line
;________________________________________________________________________________
; Compute Omega_H2_P for present fH2 and Alpha_H2_P if H2_P elastic collisions are included
;________________________________________________________________________________
      if H2_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H2_P'
         for k=0,nx-1 do begin
            DeltaVx=(VxH2(k)-vxi(k))/vth
            MagDeltaVx=abs(DeltaVx) > DeltaVx_tol
            DeltaVx=sign(DeltaVx)*MagDeltaVx
            Omega_H2_P(k)=total(Vr2pidVr*((Alpha_H2_P(*,*,k)*fH2(*,*,k))#dVx))/(nH2(k)*DeltaVx)
         endfor
         Omega_H2_P=Omega_H2_P > 0.0
      endif
;________________________________________________________________________________
; Compute Omega_H2_H for present fH2 and Alpha_H2_H if H2_H elastic collisions are included
;________________________________________________________________________________
      if H2_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H2_H'
         for k=0,nx-1 do begin
            DeltaVx=(VxH2(k)-vxH(k))/vth
            MagDeltaVx=abs(DeltaVx) > DeltaVx_tol
            DeltaVx=sign(DeltaVx)*MagDeltaVx
            Omega_H2_H(k)=total(Vr2pidVr*((Alpha_H2_H(*,*,k)*fH2(*,*,k))#dVx))/(nH2(k)*DeltaVx)
         endfor
         Omega_H2_H=Omega_H2_H > 0.0
      endif
;________________________________________________________________________________
; Compute Omega_H2_H2 for present fH2 if H2_H2 elastic collisions are included
;________________________________________________________________________________
      if H2_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H2_H2'
         if total(MH2_H2_sum) le 0 then begin
            for k=0,nx-1 do begin
               for i=0,nvr-1 do vr2_2vx_ran2(i,*)=vr(i)^2-2*(vx-VxH2(k)/vth)^2
               Wperp_paraH2(k)=total(Vr2pidVr*((vr2_2vx_ran2*fH2(*,*,k))#dVx))/nH2(k)
            endfor
         endif else begin
            for k=0,nx-1 do begin
               M_fH2=MH2_H2_sum(*,*,k)-fH2(*,*,k)
               Wperp_paraH2(k)=-total(Vr2pidVr*((vr2_2vx2_2D*M_fH2)#dVx))/nH2(k)
            endfor
         endelse
         for k=0,nx-1 do begin
            Work(*)=fH2(*,*,k)
            Alpha_H2_H2(*)=SIG_H2_H2#Work
            Wpp=Wperp_paraH2(k)
            MagWpp=abs(Wpp) > Wpp_tol
            Wpp=sign(Wpp)*MagWpp
            Omega_H2_H2(k)=total(Vr2pidVr*((Alpha_H2_H2*Work)#dVx))/(nH2(k)*Wpp)
         endfor
         Omega_H2_H2=Omega_H2_H2 > 0
      endif
   endif
;
; Total Elastic scattering frequency
;
   Omega_EL=Omega_H2_P+Omega_H2_H+Omega_H2_H2

   ;stop

;
; Total collision frequency
;
   alpha_c=dblarr(nvr,nvx,nx)
   if H2_HP_CX then begin
      for k=0,nx-1 do alpha_c(*,*,k)=alpha_cx(*,*,k)+alpha_loss(k)+Omega_EL(k)+gamma_wall(*,*,k)
   endif else begin
      for k=0,nx-1 do alpha_c(*,*,k)=alpha_loss(k)+Omega_EL(k)+gamma_wall(*,*,k)
   endelse
;   
; Test x grid spacing based on Eq.(27) in notes
;
   if debrief gt 1 then print,prompt+'Testing x grid spacing'
   Max_dx=fltarr(nx)  & Max_dx(*)=1.0e32
   for k=0,nx-1 do begin
      for j=ip(0),nvx-1 do begin 
         denom=alpha_c(*,j,k)
         Max_dx(k)= Max_dx(k) < min(2*vx(j)/denom)
      endfor
   endfor
   dx=shift(x,-1)-x
   Max_dxL=Max_dx(0:nx-2)
   Max_dxR=Max_dx(1:nx-1)
   Max_dx=Max_dxL < Max_dxR
   ilarge=where(Max_dx lt dx(0:nx-2),count)
   if count gt 0 then begin
      print,prompt+'x mesh spacing is too large!'
      debug=1
      out=''
      jj=0
      print,'   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)   x(k+1)-x(k)  Max_dx(k)'
      for ii=0,n_elements(ilarge)-1 do begin
         jj=jj+1
         out=out+string(format="('(',I3,')')",ilarge(ii))+' '+string(format="(E9.2)",x(ilarge(ii)+1)-x(ilarge(ii)))+$
                 ' '+string(format="(E9.2)",Max_dx(ilarge(ii)))+' '
         if jj gt 4 then begin
            print,out
            jj=0
            out=''
         endif
      endfor
      if jj gt 0 then print,out
      error=1
      goto,return
      ; I might need to scrap the above logic for some test cases (where it doesn't seem to work) 
   endif
;
; Define parameters Ak, Bk, Ck, Dk, Fk, Gk
;
   Ak=dblarr(nvr,nvx,nx)
   Bk=dblarr(nvr,nvx,nx)
   Ck=dblarr(nvr,nvx,nx)
   Dk=dblarr(nvr,nvx,nx)
   Fk=dblarr(nvr,nvx,nx)
   Gk=dblarr(nvr,nvx,nx)

   for k=0,nx-2 do begin
      for j=ip(0),nvx-1 do begin
         denom=2*vx(j)+(x(k+1)-x(k))*alpha_c(*,j,k+1)
         Ak(*,j,k)=(2*vx(j)-(x(k+1)-x(k))*alpha_c(*,j,k))/denom
         Bk(*,j,k)=(x(k+1)-x(k))/denom
         Fk(*,j,k)=(x(k+1)-x(k))*fw_hat(*,j)*(SH2(k+1)+SH2(k))/(vth*denom)
      endfor
   endfor
   for k=1,nx-1 do begin
      for j=0,ip(0)-1 do begin
         denom=-2*vx(j)+(x(k)-x(k-1))*alpha_c(*,j,k-1)
         Ck(*,j,k)=(-2*vx(j)-(x(k)-x(k-1))*alpha_c(*,j,k))/denom
         Dk(*,j,k)=(x(k)-x(k-1))/denom
         Gk(*,j,k)=(x(k)-x(k-1))*fw_hat(*,j)*(SH2(k)+SH2(k-1))/(vth*denom)
      endfor
   endfor


;
; Compute first-flight (0th generation) neutral distribution function
;

   
   Swall_sum=dblarr(nvr,nvx,nx)
   Beta_CX_sum=dblarr(nvr,nvx,nx)
   MH2_P_sum=dblarr(nvr,nvx,nx)
   MH2_H_sum=dblarr(nvr,nvx,nx)
   MH2_H2_sum=dblarr(nvr,nvx,nx)
   igen=0
   if debrief gt 0 then print,prompt+'Computing molecular neutral generation#'+sval(igen)
   fH2G(*,ip,0)=fH2(*,ip,0)


   for k=0,nx-2 do begin
      fH2G(*,ip,k+1)=fH2G(*,ip,k)*Ak(*,ip,k)+Fk(*,ip,k)
   endfor

   for k=nx-1,1,-1 do begin
      fH2G(*,in,k-1)=fH2G(*,in,k)*Ck(*,in,k)+Gk(*,in,k)
   endfor
;
; Compute first-flight neutral density profile
;
   for k=0,nx-1 do NH2G(k,igen)=total(Vr2pidVr*(fH2G(*,*,k)#dVx))

;
   if plot gt 1 then begin
      fH21d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH21d(*,k)=Vr2pidVr#fH2G(*,*,k)
      endfor
      plot,vx,fH21d(*,0),/nodata,yrange=[0,max(fH21d)],title='First Generation '+_HH
      for i=0,nx-1 do oplot,vx,fH21d(*,i),color=(i mod 8)+2
      if debug gt 0 then press_return
   endif
;
; Set total molecular neutral distribution function to first flight generation
;
   fH2=fH2G
   nH2=NH2G(*,0)

;
   if fH2_generations eq 0 then goto,fH2_done
   fH2_done=0

next_generation:
   if igen+1 gt max_gen then begin
      if debrief gt 0 then print,prompt+'Completed '+sval(max_gen)+' generations. Returning present solution...'
      goto,fH2_done
   endif
   igen=igen+1
   if debrief gt 0 then print,prompt+'Computing molecular neutral generation#'+sval(igen)
;________________________________________________________________________________
; Compute Swall from previous generation
;________________________________________________________________________________
   Swall=dblarr(nvr,nvx,nx)
   if total(gamma_wall) gt 0 then begin
      if debrief gt 1 then print,prompt+'Computing Swall'
      for k=0, nx-1 do begin
         Swall(*,*,k)=fw_hat*total(Vr2pidVr*((gamma_wall(*,*,k)*fH2G(*,*,k))#dVx))
      endfor
      ;for k=0,nx-1 do Swall(*,*,k)=fw_hat*total(Vr2pidVr*((gamma_wall(*,*,k)*fH2G(*,*,k))#dVx)) - ZZZ remember to put this back in because it's one line
;
; Sum wall collision source over all generations
;
      Swall_Sum=Swall_Sum+Swall
   endif

;________________________________________________________________________________
; Compute Beta_CX from previous generation
;________________________________________________________________________________
   Beta_CX=dblarr(nvr,nvx,nx)
   if H2_HP_CX then begin
      if debrief gt 1 then print,prompt+'Computing Beta_CX'
      if Simple_CX then begin
;
; Option (B): Compute charge exchange source with assumption that CX source neutrals have
;             molecular ion distribution function
;
         for k=0,nx-1 do Beta_CX(*,*,k)=fHp_hat(*,*,k)*total(Vr2pidVr*((alpha_cx(*,*,k)*fH2G(*,*,k))#dVx))
;
      endif else begin
;
; Option (A): Compute charge exchange source using fH2 and vr x sigma x v_v at each velocity mesh point
;
         for k=0,nx-1 do begin
            Work(*)=fH2G(*,*,k)
            Beta_CX(*,*,k)=nHp(k)*fHp_hat(*,*,k)*(SIG_CX#Work)
         endfor
      endelse
;
; Sum charge exchange source over all generations
;
      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   endif
;________________________________________________________________________________
; Compute MH2 from previous generation
;________________________________________________________________________________
;

   MH2_H2=dblarr(nvr,nvx,nx)
   MH2_P=dblarr(nvr,nvx,nx)
   MH2_H=dblarr(nvr,nvx,nx)
   OmegaM=dblarr(nvr,nvx,nx)
   if H2_H2_EL or H2_P_EL or H2_H_EL then begin
;
; Compute VxH2G, TH2G
;
      for k=0,nx-1 do begin
         VxH2G(k)=vth*total(Vr2pidVr*(fH2G(*,*,k)#(Vx*dVx)))/NH2G(k,igen-1)
         for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-VxH2G(k)/vth)^2
         TH2G(k)=(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fH2G(*,*,k))#dVx))/(3*q*NH2G(k,igen-1))
      endfor
      if H2_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_H2'
;
; Compute MH2_H2  
;
         vx_shift=VxH2G
         Tmaxwell=TH2G
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_H2(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_H2(k)*MH2_H2(*,*,k)
         endfor
         MH2_H2_sum=MH2_H2_sum+MH2_H2
      endif
      if H2_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_P'
;
; Compute MH2_P  
;
         vx_shift=(2*VxH2G+vxi)/3
         Tmaxwell=TH2G+(4./9.)*(Ti-TH2G +mu*mH*(vxi-VxH2G)^2/(6*q))
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_P(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_P(k)*MH2_P(*,*,k)
         endfor
         MH2_P_sum=MH2_P_sum+MH2_P
      endif
      if H2_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_H'
;
; Compute MH2_H
;
         vx_shift=(2*VxH2G+VxH)/3
         Tmaxwell=TH2G+(4./9.)*(TH-TH2G +mu*mH*(VxH-VxH2G)^2/(6*q))
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_H(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_H(k)*MH2_H(*,*,k)
         endfor
         MH2_H_sum=MH2_H_sum+MH2_H
      endif
   endif
;________________________________________________________________________________
; Compute next generation molecular distribution
;________________________________________________________________________________
;
   fH2G(*)=0.0D0
   for k=0,nx-2 do fH2G(*,ip,k+1)=Ak(*,ip,k)*fH2G(*,ip,k) $
               + Bk(*,ip,k)*(Swall(*,ip,k+1)+Beta_CX(*,ip,k+1)+OmegaM(*,ip,k+1)+Swall(*,ip,k)+Beta_CX(*,ip,k)+OmegaM(*,ip,k))

   for k=nx-1,1,-1 do fH2G(*,in,k-1)=Ck(*,in,k)*fH2G(*,in,k) $
               + Dk(*,in,k)*(Swall(*,in,k-1)+Beta_CX(*,in,k-1)+OmegaM(*,in,k-1)+Swall(*,in,k)+Beta_CX(*,in,k)+OmegaM(*,in,k))

   for k=0,nx-1 do nH2G(k,igen)=total(Vr2pidVr*(fH2G(*,*,k)#dVx))

   if plot gt 1 then begin
      fH21d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH21d(*,k)=Vr2pidVr#fH2G(*,*,k)
      endfor
      plot,vx,fH21d(*,0),/nodata,yrange=[0,max(fH21d)],title=sval(igen)+' Generation '+_HH
      for i=0,nx-1 do oplot,vx,(fH21d(*,i) > 0.9),color=(i mod 8)+2
      if debug gt 0 then press_return
   endif
;
; Add result to total neutral distribution function
;
   fH2=fH2+fH2G
   nH2=nH2+nH2G(*,igen)
;________________________________________________________________________________
; Compute 'generation error': Delta_nH2G=max(NH2G(*,igen)/max(nH2))
; and decide if another generation should be computed
;________________________________________________________________________________
   Delta_nH2G=max(NH2G(*,igen)/max(nH2))
   if fH2_Iterate then begin
      ;print, 'generation error breakpoint: ', igen
      ;STOP
;
; If fH2 'seed' is being iterated, then do another generation until the 'generation error'
; is less than 0.003 times the 'seed error' or is less than TRUNCATE
;
      if (Delta_nH2G lt 0.003*Delta_nH2s) or (Delta_nH2G lt truncate) then goto,fH2_done
   endif else begin
;
; If fH2 'seed' is NOT being iterated, then do another generation unitl the 'generation error'
; is less than parameter TRUNCATE
;
      if Delta_nH2G lt truncate then goto,fH2_done
   endelse
;________________________________________________________________________________
;
   goto,next_generation

fH2_done:
;
   if plot gt 0 then begin
      plot,x,NH2G(*,0),/nodata,yrange=[max(NH2G)*truncate,max(NH2G)],title=_HH+' Density by Generation',xtitle='x (m)',$
           ytitle='Density (m!U-3!N)',/ylog
      for i=0,igen do oplot,x,NH2G(*,i),color=(i mod 8)+2
      if pause then press_return
   endif
;
; Compute H2 density profile
   for k=0,nx-1 do nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
;
; GammaxH2 - particle flux in x direction
   for k=0,nx-1 do GammaxH2(k)=vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))
;
; VxH2 - x velocity
   VxH2=GammaxH2/nH2
   _VxH2=VxH2/vth
;
; magnitude of random velocity at each mesh point
   vr2vx2_ran=dblarr(nvr,nvx,nx)
   for i=0,nvr-1 do for k=0,nx-1 do vr2vx2_ran(i,*,k)=vr(i)^2+(vx-_VxH2(k))^2
;
; pH2 - pressure 
   for k=0,nx-1 do pH2(k)=(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH2(*,*,k))#dVx))/(3*q)
;
; TH2 - temperature
   TH2=pH2/nH2
;
; Compute NuDis - Dissociation frequency
   NuDis=n*total(sigv(*,7:10),2)
;
; Compute NuE (assume np=ne) - Energy equilibration frequency H(+) <-> H2(+)
   NuE=7.7e-7*n*1.0e-6/(sqrt(mu)*Ti^1.5)
;
; Compute H2(+) density profile
   nHP=nH2*n*sigv(*,1)/(NuDis + NuLoss)
;
; Compute THP - temperature of molecular ions
   THP=Ti*NuE/(NuE + NuDis + NuLoss)

   if fH2_Iterate then begin
;________________________________________________________________________________
; Compute 'seed error': Delta_nH2s=(|nH2s-nH2|)/max(nH2) 
; If Delta_nH2s is greater than 10*truncate then iterate fH2
;________________________________________________________________________________
      Delta_nH2s=max(abs(nH2s-nH2))/max(nH2)
      if Delta_nH2s gt 10*truncate then begin
         goto,fH2_iterate
      endif
   endif
;________________________________________________________________________________
; Update Swall_sum using last generation
;________________________________________________________________________________
;
   Swall=dblarr(nvr,nvx,nx)
   if total(gamma_wall) gt 0 then begin
      for k=0,nx-1 do Swall(*,*,k)=fw_hat*total(Vr2pidVr*((gamma_wall(*,*,k)*fH2G(*,*,k))#dVx))
      Swall_Sum=Swall_Sum+Swall
   endif

;________________________________________________________________________________
; Update Beta_CX_sum using last generation
;________________________________________________________________________________
   Beta_CX=dblarr(nvr,nvx,nx)
   if H2_HP_CX then begin
      if debrief gt 1 then print,prompt+'Computing Beta_CX'
      if Simple_CX then begin
;
; Option (B): Compute charge exchange source with assumption that CX source neutrals have
;             molecular ion distribution function
;
         for k=0,nx-1 do Beta_CX(*,*,k)=fHp_hat(*,*,k)*total(Vr2pidVr*((alpha_cx(*,*,k)*fH2G(*,*,k))#dVx))
;
      endif else begin
;
; Option (A): Compute charge exchange source using fH2 and vr x sigma x v_v at each velocity mesh point
;
         for k=0,nx-1 do begin
            Work(*)=fH2G(*,*,k)
            Beta_CX(*,*,k)=nHp(k)*fHp_hat(*,*,k)*(SIG_CX#Work)
         endfor
      endelse
      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   endif
;________________________________________________________________________________
; Update MH2_*_sum using last generation
;________________________________________________________________________________
;
   MH2_H2=dblarr(nvr,nvx,nx)
   MH2_P=dblarr(nvr,nvx,nx)
   MH2_H=dblarr(nvr,nvx,nx)
   OmegaM=dblarr(nvr,nvx,nx)
   if H2_H2_EL or H2_P_EL or H2_H_EL then begin
;
; Compute VxH2G, TH2G
;
      for k=0,nx-1 do begin
         VxH2G(k)=vth*total(Vr2pidVr*(fH2G(*,*,k)#(Vx*dVx)))/NH2G(k,igen)
         for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-VxH2G(k)/vth)^2
         TH2G(k)=(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fH2G(*,*,k))#dVx))/(3*q*NH2G(k,igen))
      endfor
      if H2_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_H2'
;
; Compute MH2_H2  
;
         vx_shift=VxH2G
         Tmaxwell=TH2G
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_H2(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_H2(k)*MH2_H2(*,*,k)
         endfor
         MH2_H2_sum=MH2_H2_sum+MH2_H2
      endif
      if H2_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_P'
;
; Compute MH2_P  
;
         vx_shift=(2*VxH2G+vxi)/3
         Tmaxwell=TH2G+(4./9.)*(Ti-TH2G +mu*mH*(vxi-VxH2G)^2/(6*q))
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_P(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_P(k)*MH2_P(*,*,k)
         endfor
         MH2_P_sum=MH2_P_sum+MH2_P
      endif
      if H2_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH2_H'
;
; Compute MH2_H
;
         vx_shift=(2*VxH2G+VxH)/3
         Tmaxwell=TH2G+(4./9.)*(TH-TH2G +mu*mH*(VxH-VxH2G)^2/(6*q))
         mol=2
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH2_H(*,*,k)=Maxwell(*,*,k)*NH2G(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H2_H(k)*MH2_H(*,*,k)
         endfor
         MH2_H_sum=MH2_H_sum+MH2_H
      endif
   endif
;________________________________________________________________________________
; Compute remaining moments
;________________________________________________________________________________
; piH2_xx
   for k=0,nx-1 do piH2_xx(k)=(2*mu*mH)*vth2*total(Vr2pidVr*(fH2(*,*,k)#(dVx*(vx-_VxH2(k))^2)))/q - pH2(k)
;
; piH2_yy
   for k=0,nx-1 do piH2_yy(k)=(2*mu*mH)*vth2*0.5*total((Vr2pidVr*Vr^2)*(fH2(*,*,k)#dVx))/q - pH2(k)
;
; piH2_zz
   piH2_zz=piH2_yy
;
; qxH2
   for k=0,nx-1 do qxH2(k)=0.5*(2*mu*mH)*vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH2(*,*,k))#(dVx*(vx-_VxH2(k)))))
;________________________________________________________________________________
; C = RHS of Boltzman equation for total fH2
;________________________________________________________________________________
   for k=0,nx-1 do begin
      C=vth*(fw_hat(*,*)*SH2(k)/vth + Swall_sum(*,*,k) + Beta_CX_sum(*,*,k) - alpha_c(*,*,k)*fH2(*,*,k) + $
             Omega_H2_P(k)*MH2_P_sum(*,*,k)+Omega_H2_H(k)*MH2_H_sum(*,*,k)+Omega_H2_H2(k)*MH2_H2_sum(*,*,k))
      QH2(k)=0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C)#dVx))
      RxH2(k)=(2*mu*mH)*vth*total(Vr2pidVr*(C#(dVx*(vx-_VxH2(k)))))
      Sloss(k)=-total(Vr2pidVr*(C#dVx))+SH2(k)
      WallH2(k)=total(Vr2pidVr*((gamma_wall(*,*,k)*fH2(*,*,k))#dVx))
      if H2_H_EL then begin
         CH2_H=vth*Omega_H2_H(k)*(MH2_H_sum(*,*,k)-fH2(*,*,k))
         RxH_H2(k)=(2*mu*mH)*vth*total(Vr2pidVr*(CH2_H#(dVx*(vx-_VxH2(k)))))
         EH_H2(k)=0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CH2_H)#dVx))
      endif
      if H2_P_EL then begin
         CH2_P=vth*Omega_H2_P(k)*(MH2_P_sum(*,*,k)-fH2(*,*,k))
         RxP_H2(k)=(2*mu*mH)*vth*total(Vr2pidVr*(CH2_P#(dVx*(vx-_VxH2(k)))))
         EP_H2(k)=0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CH2_P)#dVx))
      endif
      if H2_HP_CX then begin
         CH2_HP_CX=vth*(Beta_CX_sum(*,*,k) - alpha_cx(*,*,k)*fH2(*,*,k))
         RxH2CX(k)=(2*mu*mH)*vth*total(Vr2pidVr*(CH2_HP_CX#(dVx*(vx-_VxH2(k)))))
         EH2CX(k)=0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CH2_HP_CX)#dVx))
      endif
      CW_H2=vth*(Swall_sum(*,*,k) - gamma_wall(*,*,k)*fH2(*,*,k))
      RxW_H2(k)=(2*mu*mH)*vth*total(Vr2pidVr*(CW_H2#(dVx*(vx-_VxH2(k)))))
      EW_H2(k)=0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CW_H2)#dVx))
      if H2_H2_EL then begin
         CH2_H2=vth*Omega_H2_H2(k)*(MH2_H2_sum(*,*,k)-fH2(*,*,k))
         for i=0,nvr-1 do vr2_2vx_ran2(i,*)=vr(i)^2-2*(vx-_VxH2(k))^2
         Epara_PerpH2_H2(k)=-0.5*(2*mu*mH)*vth2*total(Vr2pidVr*((vr2_2vx_ran2*CH2_H2)#dVx))
      endif
   endfor
;
; qxH2_total
;
   qxH2_total=(0.5*nH2*(2*mu*mH)*VxH2*VxH2 + 2.5*pH2*q)*VxH2 + q*piH2_xx*VxH2 + qxH2
;
; QH2_total
;
   QH2_total=QH2+RxH2*VxH2 - 0.5*(2*mu*mH)*(Sloss-SH2)*VxH2*VxH2
;
; Albedo
;
   AlbedoH2=0.0
   gammax_plus=vth*total(Vr2pidVr*(fH2(*,ip,0)#(Vx(ip)*dVx(ip))))
   gammax_minus=vth*total(Vr2pidVr*(fH2(*,in,0)#(Vx(in)*dVx(in))))
   if abs(gammax_plus) gt 0 then AlbedoH2=-gammax_minus/gammax_plus
;
;________________________________________________________________________________
; Compute Mesh Errors
;________________________________________________________________________________
   mesh_error=fltarr(nvr,nvx,nx)
   max_mesh_error=0.0
   min_mesh_error=0.0
   mtest=5
   moment_error=fltarr(nx,mtest)
   max_moment_error=fltarr(mtest)
   C_error=fltarr(nx)
   CX_error=fltarr(nx)
   Wall_error=fltarr(nx)
   H2_H2_error=fltarr(nx,3)
   H2_H_error=fltarr(nx,3)
   H2_P_error=fltarr(nx,3)
   max_H2_H2_error=fltarr(3)
   max_H2_H_error=fltarr(3)
   max_H2_P_error=fltarr(3)
;
   if compute_errors then begin
      if debrief gt 1 then print,prompt+'Computing Collision Operator, Mesh, and Moment Normalized Errors'
;
      Sloss2=vth*Alpha_Loss*nH2
      for k=0,nx-1 do C_error(k)=abs(Sloss(k)-Sloss2(k))/max(abs([Sloss(k),Sloss2(k)]))
;
; Test conservation of particles for charge exchange operator
;
      if H2_HP_CX then begin
         for k=0,nx-1 do begin
            CX_A=total(Vr2pidVr*((alpha_cx(*,*,k)*fH2(*,*,k))#dVx))
            CX_B=total(Vr2pidVr*((Beta_CX_sum(*,*,k))#dVx))
            CX_error(k)=abs(CX_A-CX_B)/max(abs([CX_A,CX_B]))
         endfor
      endif
;
; Test conservation of particles for wall collision operator
;
      if total(PipeDia) gt 0 then begin
         for k=0,nx-1 do begin
            Wall_A=WallH2(k)
            Wall_B=total(Vr2pidVr*((Swall_sum(*,*,k))#dVx))
            if max(abs([Wall_A,Wall_B])) gt 0 then Wall_error(k)=abs(Wall_A-Wall_B)/max(abs([Wall_A,Wall_B]))
         endfor
      endif
;
; Test conservation of particles, x momentum, and total energy of elastic collision operators
;
      for m=0,2 do begin
         for k=0,nx-1 do begin
            if m lt 2 then begin
               TfH2=total(Vr2pidVr*(fH2(*,*,k)#(dVx*Vx^m)))
            endif else begin
               TfH2=total(Vr2pidVr*((vr2vx2(*,*,k)*fH2(*,*,k))#dVx))
            endelse
            if H2_H2_EL then begin
               if m lt 2 then begin
                  TH2_H2=total(Vr2pidVr*(MH2_H2_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH2_H2=total(Vr2pidVr*((vr2vx2(*,*,k)*MH2_H2_sum(*,*,k))#dVx))
               endelse
               H2_H2_error(k,m)=abs(TfH2-TH2_H2)/max(abs([TfH2,TH2_H2]))
            endif
            if H2_H_EL then begin
               if m lt 2 then begin
                  TH2_H=total(Vr2pidVr*(MH2_H_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH2_H=total(Vr2pidVr*((vr2vx2(*,*,k)*MH2_H_sum(*,*,k))#dVx))
               endelse
               H2_H_error(k,m)=abs(TfH2-TH2_H)/max(abs([TfH2,TH2_H]))
            endif
            if H2_P_EL then begin
               if m lt 2 then begin
                  TH2_P=total(Vr2pidVr*(MH2_P_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH2_P=total(Vr2pidVr*((vr2vx2(*,*,k)*MH2_P_sum(*,*,k))#dVx))
               endelse
               H2_P_error(k,m)=abs(TfH2-TH2_P)/max(abs([TfH2,TH2_P]))
            endif
         endfor
         max_H2_H2_error(m)=max(H2_H2_error(*,m))
         max_H2_H_error(m)= max(H2_H_error(*,m))
         max_H2_P_error(m)= max(H2_P_error(*,m))
      endfor
;
      if CI_test then begin
          minRx=1.0e-6
          minEpara_perp=1.0e-6
;
; Compute Momentum transfer rate via full collision integrals for charge exchange and mixed elastic scattering
; Then compute error between this and actual momentum transfer resulting from CX and BKG (elastic) models
;
         if H2_HP_CX then begin &; H2(+) -> H2 charge exchange momentum transfer via full collision integral
            print,prompt+'Computing H2(+) -> H2 Charge Exchange Momentum Transfer'
            _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
            _Sig(*)=v_v*sigma_cx_hh(v_v2*(mH*vth2/q))
            SIG_VX_CX=dblarr(nvr*nvx,nvr*nvx)
            SIG_VX_CX(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)
            alpha_vx_cx=dblarr(nvr,nvx,nx)
            for k=0,nx-1 do begin
               Work(*)=nHp(k)*fHp_hat(*,*,k)
               alpha_vx_cx(*,*,k)=SIG_VX_CX#Work
            endfor
            for k=0,nx-1 do RxCI_CX(k)=-(2*mu*mH)*vth2*total(Vr2pidVr*((Alpha_vx_cx(*,*,k)*fH2(*,*,k))#dVx))
            norm=max(abs([RxH2CX,RxCI_CX]))
            for k=0,nx-1 do CI_CX_error(k)=abs(RxH2CX(k)-RxCI_CX(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in CX collision operator: ',sval(max(CI_CX_Error))
         endif
         if H2_P_EL then begin &; P -> H2 momentum transfer via full collision integral
            for k=0,nx-1 do RxCI_P_H2(k)=-(1.0/3.0)*(2*mu*mH)*vth2*total(Vr2pidVr*((Alpha_H2_P(*,*,k)*fH2(*,*,k))#dVx))
            norm=max(abs([RxP_H2,RxCI_P_H2]))
            for k=0,nx-1 do CI_P_H2_error(k)=abs(RxP_H2(k)-RxCI_P_H2(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in P -> H2 elastic BKG collision operator: ',sval(max(CI_P_H2_Error))
         endif
         if H2_H_EL then begin &; H -> H2 momentum transfer via full collision integral
            for k=0,nx-1 do RxCI_H_H2(k)=-(1.0/3.0)*(2*mu*mH)*vth2*total(Vr2pidVr*((Alpha_H2_H(*,*,k)*fH2(*,*,k))#dVx))
            norm=max(abs([RxH_H2,RxCI_H_H2]))
            for k=0,nx-1 do CI_H_H2_error(k)=abs(RxH_H2(k)-RxCI_H_H2(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in H -> H2 elastic BKG collision operator: ',sval(max(CI_H_H2_Error))
         endif
         if H2_H2_EL then begin &; H2 -> H2 perp/parallel energy transfer via full collision integral
            for k=0,nx-1 do begin
               Work(*)=fH2(*,*,k)
               Alpha_H2_H2(*)=SIG_H2_H2#Work
               Epara_Perp_CI(k)=0.5*(2*mu*mH)*vth3*total(Vr2pidVr*((Alpha_H2_H2*fH2(*,*,k))#dVx))
            endfor
            norm=max(abs([Epara_PerpH2_H2,Epara_Perp_CI]))
            for k=0,nx-1 do CI_H2_H2_error(k)=abs(Epara_PerpH2_H2(k)-Epara_Perp_CI(k))/norm
            print,prompt+'Maximum normalized perp/parallel energy transfer error in H2 -> H2 elastic BKG collision operator: ',$
                  sval(max(CI_H2_H2_Error))
         endif
      endif
;
;  Mesh Point Error based on fH2 satisfying Boltzmann equation
;
      T1=dblarr(nvr,nvx,nx) & T2=T1 & T3=T1 & T4=T1 & T5=T1 & T6=T1
      for k=0,nx-2 do begin
         for j=0,nvx-1 do T1(*,j,k)=2*vx(j)*(fH2(*,j,k+1)-fH2(*,j,k))/(x(k+1)-x(k))
         T2(*,*,k)=fw_hat(*,*)*(SH2(k+1)+SH2(k))/vth
         T3(*,*,k)=Beta_CX_sum(*,*,k+1)+Beta_CX_sum(*,*,k)
         T4(*,*,k)=alpha_c(*,*,k+1)*fH2(*,*,k+1)+alpha_c(*,*,k)*fH2(*,*,k)
         T5(*,*,k)=Omega_H2_P(k+1)*MH2_P_sum(*,*,k+1)+Omega_H2_H(k+1)*MH2_H_sum(*,*,k+1)+Omega_H2_H2(k+1)*MH2_H2_sum(*,*,k+1)+$
                   Omega_H2_P(k)*MH2_P_sum(*,*,k)+Omega_H2_H(k)*MH2_H_sum(*,*,k)+Omega_H2_H2(k)*MH2_H2_sum(*,*,k)
         T6(*,*,k)=Swall_sum(*,*,k+1)+Swall_sum(*,*,k)
         Mesh_Error(*,*,k)=abs(T1(*,*,k)-T2(*,*,k)-T3(*,*,k)+T4(*,*,k)-T5(*,*,k)-T6(*,*,k))/$
                              max(abs([T1(*,*,k),T2(*,*,k),T3(*,*,k),T4(*,*,k),T5(*,*,k),T6(*,*,k)]))
      endfor
      ave_mesh_error=total(mesh_error)/n_elements(mesh_error)
      max_mesh_error=max(mesh_error)
      min_mesh_error=min(mesh_error(*,*,0:nx-2))
;
;  Moment Error
;
      for m=0,mtest-1 do begin
         for k=0,nx-2 do begin
            MT1=total(Vr2pidVr*(T1(*,*,k)#(dVx*Vx^m)))
            MT2=total(Vr2pidVr*(T2(*,*,k)#(dVx*Vx^m)))
            MT3=total(Vr2pidVr*(T3(*,*,k)#(dVx*Vx^m)))
            MT4=total(Vr2pidVr*(T4(*,*,k)#(dVx*Vx^m)))
            MT5=total(Vr2pidVr*(T5(*,*,k)#(dVx*Vx^m)))
            MT6=total(Vr2pidVr*(T6(*,*,k)#(dVx*Vx^m)))
            moment_error(k,m)=abs(MT1-MT2-MT3+MT4-MT5-MT6)/max(abs([MT1,MT2,MT3,MT4,MT5,MT6]))
         endfor
         max_moment_error(m)=max(moment_error(*,m))
      endfor
;
; Compute error in qxH2_total
;
;    qxH2_total2 total neutral heat flux profile (watts m^-2)
;               This is the total heat flux transported by the neutrals
;               computed in a different way from:
;
;               qxH2_total2(k)=vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH2(*,*,k))#(Vx*dVx)))*0.5*(2*mu*mH)
;
;               This should agree with qxH2_total if the definitions of nH2, pH2, piH2_xx,
;               TH2, VxH2, and qxH2 are coded correctly.
;
      qxH2_total2=dblarr(nx)
      for k=0,nx-1 do qxH2_total2(k)=0.5*(2*mu*mH)*vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH2(*,*,k))#(Vx*dVx)))
      qxH2_total_error=abs(qxH2_total-qxH2_total2)/max(abs([qxH2_total,qxH2_total2]))
;
; Compute error in QH2_total
;
      Q1=dblarr(nx)
      Q2=dblarr(nx)
      QH2_total_error=fltarr(nx)
      for k=0,nx-2 do begin
         Q1(k)=(qxH2_total(k+1)-qxH2_total(k))/(x(k+1)-x(k))
         Q2(k)=0.5*(QH2_total(k+1)+QH2_total(k))
      endfor
      QH2_total_error=abs(Q1-Q2)/max(abs([Q1,Q2]))
;
      if debrief gt 0 then begin
         print,prompt+'Maximum particle convervation error of total collision operator: ',sval(max(C_Error))
         print,prompt+'Maximum H2_HP_CX particle convervation error: ',sval(max(CX_Error))
         print,prompt+'Maximum H2_Wall particle convervation error: ',sval(max(Wall_Error))
         print,prompt+'Maximum H2_H2_EL particle conservation error: '+sval(max_H2_H2_error(0))
         print,prompt+'Maximum H2_H2_EL x-momentum conservation error: '+sval(max_H2_H2_error(1))
         print,prompt+'Maximum H2_H2_EL total energy conservation error: '+sval(max_H2_H2_error(2))
         print,prompt+'Maximum H2_H_EL  particle conservation error: '+sval(max_H2_H_error(0))
         print,prompt+'Maximum H2_P_EL  particle conservation error: '+sval(max_H2_P_error(0))
         print,prompt+'Average mesh_error =',ave_mesh_error
         print,prompt+'Maximum mesh_error =',max_mesh_error
         for m=0,4 do print,prompt+'Maximum fH2 vx^'+sval(m)+' moment error: '+sval(max_moment_error(m))
         print,prompt+'Maximum qxH2_total error =',max(qxH2_total_Error)
         print,prompt+'Maximum QH2_total error =',max(QH2_total_Error)
         if debug gt 0 then press_return
      endif
   endif

   mid1=locate(x,0.7*(max(x)+min(x))/2)
   mid2=locate(x,0.85*(max(x)+min(x))/2)
   mid3=locate(x,(max(x)+min(x))/2)
   mid4=locate(x,1.15*(max(x)+min(x))/2)
   mid5=locate(x,1.3*(max(x)+min(x))/2)
   mid6=locate(x,1.45*(max(x)+min(x))/2)

   if plot gt 1 then begin
      fH21d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH21d(*,k)=Vr2pidVr#fH2(*,*,k)
      endfor
      plot,vx,fH21d(*,0),/nodata,title=_HH+' Velocity Distribution Function: fH2(Vx)',xtitle='Vx/Vth'
      for i=0,nx-1 do oplot,vx,fH21d(*,i),color=(i mod 6)+2
      if pause then press_return
   endif

   if plot gt 0 then begin
      data=[nH,n,nHP,nH2]
      jp=where(data gt 0)
      yrange=[min(data(jp)),max(data(jp))]
      plot,x,nH,/nodata,/ylog,yrange=yrange,title='Density Profiles',xtitle='x (meters)',ytitle='m!U-3!N'
      oplot,x,nH,color=2
      xyouts,x(mid1),1.2*nH(mid1),_H,color=2
      oplot,x,n,color=3
      xyouts,x(mid2),1.2*n(mid2),'e!U-!N',color=3
      oplot,x,nH2,color=4
      xyouts,x(mid3),1.2*nH2(mid3),_HH,color=4
      oplot,x,nHP,color=6
      xyouts,x(mid4),1.2*nHP(mid4),_Hp,color=6
      if pause then press_return
   endif
;
   if plot gt 0 then begin
      data=[TH,Te,THP,TH2]
      jp=where(data gt 0)
      yrange=[min(data(jp)),max(data(jp))]
      plot,x,TH,/nodata,/ylog,yrange=yrange,title='Temperature Profiles',xtitle='x (meters)',ytitle='eV'
      oplot,x,TH,color=2
      xyouts,x(mid1),1.2*TH(mid1),_H,color=2
      oplot,x,Te,color=3
      xyouts,x(mid2),1.2*Te(mid2),'e!U-!N',color=3
      oplot,x,TH2,color=4
      xyouts,x(mid3),1.2*TH2(mid3),_HH,color=4
      oplot,x,THP,color=6
      xyouts,x(mid4),1.2*THP(mid4),_Hp,color=6
      if pause then press_return
   endif

   if Compute_H_Source then begin
      if debrief gt 1 then print,prompt+'Computing Velocity Distributions of H products...'
;
; Set Normalized Franck-Condon Velocity Distributions for reactions R2, R3, R4, R5, R6, R7, R8, R10
;
; Make lookup table to select reaction Rn in SFCn
;        Rn=2 3 4 5 6 7 8   10
      nFC=[0,0,0,1,2,3,4,5,6,0,7]
;
      SFCn=dblarr(nvr,nvx,nx,8)
      Eave=dblarr(nx,8)
      Emax=dblarr(nx,8)
      Emin=dblarr(nx,8)
;
;  Reaction R2: e + H2 -> e + H(1s) + H(1s)
;
      ii=nFC(2)
      Eave(*,ii)=3.0
      Emax(*,ii)=4.25
      Emin(*,ii)=2
;
;  Reaction R3: e + H2 -> e + H(1s) + H*(2s)
;
      ii=nFC(3)
      Eave(*,ii)=0.3
      Emax(*,ii)=0.55
      Emin(*,ii)=0.0
;
; Reaction R4:  e + H2 -> e + H(+) + H(1s) + e
;
      ii=nFC(4)
      Ee=3*Te/2			&;Note the FC energy depends on electron energy
      kk=where(Ee le 26.0,count)
      if count gt 0 then Eave(kk,ii)=0.25
      kk=where(Ee gt 26.0 and Ee le 41.6,count)
      if count gt 0 then begin
         Eave(kk,ii)=0.5*(Ee(kk)-26)
         Eave(kk,ii)=Eave(kk,ii) > 0.25
      endif
      kk=where(Ee gt 41.6,count)
      if count gt 0 then Eave(kk,ii)=7.8
      Emax(*,ii)=1.5*Eave(*,ii)	&;Note the max/min values here are a guess
      Emin(*,ii)=0.5*Eave(*,ii)	&;Note the max/min values here are a guess
;
;  Reaction R5: e + H2 -> e + H*(2p) + H*(2s)
;
      ii=nFC(5)
      Eave(*,ii)=4.85
      Emax(*,ii)=5.85
      Emin(*,ii)=2.85
;
;  Reaction R6: e + H2 -> e + H(1s) + H*(n=3)
;
      ii=nFC(6)
      Eave(*,ii)=2.5
      Emax(*,ii)=3.75
      Emin(*,ii)=1.25
;
;  Reaction R7: e + H2(+) -> e + H(+) + H(1s)
;
      ii=nFC(7)
      Eave(*,ii)=4.3
      Emax(*,ii)=4.3+2.1	  &;Note the max/min values here are a guess
      Emin(*,ii)=4.3-2.1     &;Note the max/min values here are a guess
;
;  Reaction R8: e + H2(+) -> e + H(+) + H*(n=2)
;
      ii=nFC(8)
      Eave(*,ii)=1.5
      Emax(*,ii)=1.5+0.75	  &;Note the max/min values here are a guess
      Emin(*,ii)=1.5-0.75    &;Note the max/min values here are a guess
;
;  Reaction R10: e + H2(+) -> H(1s) + H*(n>=2)
;
      ii=nFC(10)
;
; Compute relative cross-sections for populating a specific n level for reaction R10
; (see page 62 in Janev, "Elementary Processes in Hydrogen-Helium Plasmas", Springer-Verlag, 1987)
;         n=2   3    4    5    6
      R10rel=[0.1,0.45,0.22,0.12,0.069]
      for k=7,10 do R10rel=[R10rel,10.0/k^3]
      En=13.58/(2+Findgen(9))^2	&;energy of levels
;
      for k=0,nx-1 do begin
         EHn=0.5*(Ee-En)*R10rel/total(R10rel)
         EHn=EHn > 0.0
         Eave(k,ii)=total(EHn)
         Eave(k,ii)=Eave(k,ii) > 0.25
         Emax(k,ii)=1.5*Eave(k,ii)  &;Note the max/min values here are a guess
         Emin(k,ii)=0.5*Eave(k,ii)  &;Note the max/min values here are a guess
      endfor
;
; Set SFCn values for reactions R2, R3, R4, R5, R6, R7, R8, R10
;
      Vfc=dblarr(nvr,nvx,nx)
      Tfc=dblarr(nvr,nvx,nx)
      magV=sqrt(vr2vx2)
      _THP=dblarr(nvr,nvx,nx)
      _TH2=dblarr(nvr,nvx,nx)
      for k=0,nx-1 do begin
         _THP(*,*,k)=THP(k)/Tnorm
         _TH2(*,*,k)=TH2(k)/Tnorm
      endfor
;
; The following function is choosen to represent the velocity distribution of the
; hydrogen products for a given reaction, accounting for the Franck-Condon energy 
; distribution and accounting for additional velocity spreading due to the finite 
; temperature of the molcules (neutral and ionic) prior to breakup:
;
;     f(Vr,Vx) = exp( -0.5*mH*mu*(|v|-Vfc+0.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )
;
;       	|v|=sqrt(Vr^2+Vx^2)
;	        Tfc= Franck-Condon 'temperature' = (Emax-Emin)/4
;	        Vfc= Franck-Condon  velocity = sqrt(2 Eave/mH/mu)
;		Tmol= temperature of H2 molecule (neutral or ionic)
;
;    This function is isotropic in velocity space and can be written in terms
;  of a distribution in particle speed, |v|, 
;
;     f(|v|) = exp( -(|v|-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )
;
; with velocities normalized by vth and T normalized by Tnorm.
;
;  Recognizing the the distribution in energy, f(E), and particle speed, f(|v|),
;  are related by  f(E) dE = f(|v|) 2 pi v^2 dv, and that dE/dv = mH mu v,
;  f(E) can be written as
;
;     f(E) = f(|v|) 2 pi |v|/(mH mu) = const. |v| exp( -(|v|-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*Tmol) )
;
; The function f(Vr,Vx) was chosen because it has has the following characteristics:
;
; (1) For Tmol/2 << Tfc,  the peak in the v^2 times the energy distribution, can be found
;    by finding the |v| where df(E)/d|v| =0
;
;    df(E)/d|v|= 0 = 3v^2 exp() - 2(|v|-Vfc+1.5*Tfc/Vfc)/Tfc v^3 exp() 
;                    2(|v|-Vfc+1.5*Tfc/Vfc)/Tfc |v| = 3
;
;    which is satisfied when |v|=Vfc. Thus the energy-weighted energy distribution peaks
;    at the velocity corresponding to the average Franck-Condon energy.
;
; (2) for Tmol/2 >> Tfc ~ Vfc^2, the velocity distribution becomes
;
;	f(|v|) = exp( -2(|v|-Vfc+1.5*Tfc/Vfc)^2/Tmol )
;
;    which leads to a velocity distribution that approaches the molecular velocity
;    distribution with the magnitude of the average velocity divided by 2. This
;    is the appropriate situation for when the Franck-Condon energies are negligible
;    relative to the thermal speed of the molecules.
;
      Rn=[2,3,4,5,6,7,8,10]
      for jRn=0,n_elements(Rn)-1 do begin
         ii=nFC(Rn(jRn))
         Tfc(0,0,*)=0.25*(Emax(*,ii)-Emin(*,ii))/Tnorm	&;Franck-Condon 'effective temperature'
         Vfc(0,0,*)=sqrt(Eave(*,ii)/Tnorm)     		&;Velocity corresponding to Franck-Condon 'mean evergy'
         for k=0,nx-1 do begin
            Vfc(*,*,k)=Vfc(0,0,k)
            Tfc(*,*,k)=Tfc(0,0,k)
         endfor
         if Rn(jRn) le 6 then begin
;
;   For R2-R6, the Franck-Condon 'mean energy' is taken equal to Eave
;	   and the 'temperature' corresponds to the sum of the Franck-Condon 'temperature', Tfc,
;          and the temperature of the H2 molecules, TH2. (Note: directed neutral molecule velocity
;	   is not included and assumed to be small)
;          
            arg=-(magV-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*_TH2)
            SFCn(*,*,*,ii)=exp( arg > (-80) )
         endif else begin
;
;   For R7, R8 and R10, the Franck-Condon 'mean energy' is taken equal to Eave
;	   and the 'temperature' corresponds to the sum of the Franck-Condon 'temperature', Tfc,
;          and the temperature of the H2(+) molecular ions, THP. (Note: directed molecular ion velocity
;	   is not included and assumed to be small)
;          
            arg=-(magV-Vfc+1.5*Tfc/Vfc)^2/(Tfc+0.5*_THP)
            SFCn(*,*,*,ii)=exp( arg > (-80) )
         endelse
         for k=0,nx-1 do begin
            integral = total(Vr2pidVr*(SFCn(*,*,k,ii)#dVx))
            SFCn(*,*,k,ii)=SFCn(*,*,k,ii) / integral
             ; SFCn(*,*,k,ii)=SFCn(*,*,k,ii)/(total(Vr2pidVr*(SFCn(*,*,k,ii)#dVx))) - ZZZZ this is the original line
         endfor
      endfor
      nm=3
      if plot gt 3 then begin
         ii=nFC(Rn(0))
         m=indgen(nm)*(nx-1)/(nm-1)
         kp=0
         contour,SFCn(*,*,m(0),ii)
         for mm=0,nm-1 do begin
            for jRn=0,n_elements(Rn)-1 do begin
               ii=nFC(Rn(jRn))
               color=(kp mod 6)+1
               contour,SFCn(*,*,m(mm),ii),/noerase,color=color
               xyouts,.7,.8-kp*.03,/norm,'R'+sval(Rn(jRn))+' Te:'+sval(Te(m(mm))),color=color
               kp=kp+1
               if pause then press_return
            endfor
            ii=nFC(Rn(0))
            contour,SFCn(*,*,0,ii),/noerase,title='SFCn - Franck-Condon Velocity Distributions'
         endfor
      endif

      if plot gt 2 then begin
         m=indgen(nm)*(nx-1)/(nm-1)
         kp=0
         plot,/nodata,/xlog,Eaxis,Eaxis,yrange=[0,1],title='Energy Distribution of '+_H+' Products',xtitle='Energy (eV)'
         for mm=0,2 do begin
            kp=kp+1
            xyouts,.85,0.92-kp*.02,/norm,charsize=0.8,'Te:'+sval(Te(m(mm)),l=4)
            kp=kp+1
            for jRn=0,n_elements(Rn)-1 do begin
               ii=nFC(Rn(jRn))
               color=(jRn mod n_elements(Rn))+1
               EFC=Eaxis*SFCn(*,ip(0),m(mm),ii)*VrVr4pidVr/dEaxis
               EFC=EFC/Max(EFC)
               oplot,Eaxis,EFC,color=color
               xyouts,.87,.92-kp*.02,/norm,charsize=0.8,'R'+sval(Rn(jRn)),color=color
               kp=kp+1
            endfor
         endfor
         if pause then press_return
      endif

      if plot gt 0 then begin
         kp=0
         for jRn=0,n_elements(Rn)-1 do begin
            ii=nFC(Rn(jRn))
            Ebar=dblarr(nx)
            for k=0,nx-1 do Ebar(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*SFCn(*,*,k,ii))#dVx))/q
            color=(kp mod n_elements(Rn))+1
            if kp eq 0 then plot,x,Ebar,title='Average Energy of '+_H+', '+_p+' Products',xtitle='x (meters)',ytitle='eV',$
                                 yrange=[0.1,100],/ylog
            oplot,x,Ebar,color=color
            xyouts,.4,.90-.03*kp,/norm,'R'+sval(Rn(jRn))+': '+_rn(Rn(jRn)),color=color,charsize=0.8
            kp=kp+1
         endfor
         if pause then press_return
      endif
   
      Vbar_Error=fltarr(nx)
      if compute_errors then begin
;
; Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti(x)/Tnorm)/sqrt(!pi)
;
         
         TFC=min(Eave(0,*))+(max(Eave(0,*))-min(Eave(0,*)))*findgen(nx)/(nx-1)
         vx_shift=TFC & vx_shift(*)=0.0
         Tmaxwell=TFC
         mol=1
@create_shifted_maxwellian.include
         vbar_test=vth*sqrt(vr2vx2(*,*,0))
         for k=0,nx-1 do begin 
            vbar=total(Vr2pidVr*((vbar_test*Maxwell(*,*,k))#dVx))
            vbar_exact=2*Vth*sqrt(TFC(k)/Tnorm)/sqrt(!pi) 
            Vbar_Error(k)=abs(vbar-vbar_exact)/vbar_exact
         endfor
         if debrief gt 0 then begin
            print,prompt+'Maximum Vbar error over FC energy range =',max(Vbar_Error)
         endif
      endif
;
; Compute atomic hydrogen source distribution function
; using normalized FC source distributions SFCn
;
      for k=0,nx-1 do begin
         fSH(*,*,k)=n(k)*nH2(k)*( 2*sigv(k,2)*SFCn(*,*,k,nFC(2))+$
                                   2*sigv(k,3)*SFCn(*,*,k,nFC(3))+$
                                     sigv(k,4)*SFCn(*,*,k,nFC(4))+$
                                   2*sigv(k,5)*SFCn(*,*,k,nFC(5))+$
                                   2*sigv(k,6)*SFCn(*,*,k,nFC(6)) )
         fSH(*,*,k)=fSH(*,*,k)+$
                     n(k)*nHP(k)*(   sigv(k,7)* SFCn(*,*,k,nFC(7))+$
                                     sigv(k,8)* SFCn(*,*,k,nFC(8))+$
                                   2*sigv(k,10)*SFCn(*,*,k,nFC(10)) )
      endfor
;
; Compute total H and H(+) sources
;
      for k=0,nx-1 do begin
         SH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
         SP(k)=n(k)*nH2(k)*sigv(k,4)+$
               n(k)*nHP(k)*( sigv(k,7) + sigv(k,8) + 2*sigv(k,9) )
      endfor
;
; Compute total HP source
;
      SHP=n*nH2*sigv(*,1)
;
; Compute energy distribution of H source
;
      for k=0,nx-1 do begin
         ESH(*,k)=Eaxis*fSH(*,ip(0),k)*VrVr4pidVr/dEaxis
         ESH(*,k)=ESH(*,k)/max(ESH(*,k))
      endfor
;
      if plot gt 2 then begin
         fH21d=fltarr(nvx,nx)
         for k=0,nx-1 do begin
            fH21d(*,k)=Vr2pidVr#fSH(*,*,k)
         endfor
         plot,vx,fH21d(*,0),/nodata,title=_H+' Source Velocity Distribution Function: fSH(Vx)',$
                   ytitle='m!U-3!N s!U-1!N dVx!U-1!N',xtitle='Vx/Vth'
         for i=0,nx-1 do oplot,vx,fH21d(*,i),color=(i mod 6)+2
         if pause then press_return
      endif

      if plot gt 2 then begin
         plot,/xlog,Eaxis,ESH(*,0),/nodata,title=_H+' Source Energy Distribution: ESH(E)= E fSH(E)',xtitle='E (eV)'
         for k=0,nx-1 do oplot,Eaxis,ESH(*,k),color=(k mod 6)+2
         if pause then press_return
      endif
      
      if plot gt 1 then begin
         Ebar=dblarr(nx)
         for k=0,nx-1 do Ebar(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*fSH(*,*,k))#dVx))/q
         for k=0,nx-1 do Ebar(k)=Ebar(k)/(total(Vr2pidVr*(fSH(*,*,k)#dVx)))
         plot,x,Ebar,title='Average Energy of '+_H+' Source Distribution',xtitle='x (meters)',ytitle='eV',/ylog,/nodata
         oplot,x,Ebar,color=2
         if pause then press_return
      endif

      if plot gt 0 then begin
         data=[SH,SP,SHP,SH2,NuLoss*nHP,NuDis*nHP]
         jp=where(data gt 0)
         yrange=[min(data(jp)),max(data(jp))]
         plot,x,SH,/nodata,/ylog,yrange=yrange,title='Source and Sink Profiles',$
              ytitle='m!U-3!N s!U-1!N',xtitle='x (meters)'
         oplot,x,SH,color=2
         xyouts,x(mid1),1.2*SH(mid1),_H+' source',color=2
         oplot,x,SHP,color=3
         xyouts,x(mid2),1.2*SHP(mid2),_Hp+' source',color=3
         oplot,x,SP,color=4
         xyouts,x(mid3),0.6*SP(mid3),_p+' source',color=4
         oplot,x,NuLoss*nHp,color=5
         xyouts,x(mid4),0.6*NuLoss(mid4)*nHp(mid4),_Hp+' loss',color=5
         oplot,x,NuDis*nHp,color=6
         xyouts,x(mid5),0.6*NuDis(mid5)*nHp(mid5),_Hp+' dissoc.',color=6
         oplot,x,SH2,color=1
         xyouts,x(mid6),0.8*SH2(mid6),_HH+' source',color=1
         if pause then press_return
      endif
      if plot gt 0 then begin
         gammaxH2_plus=dblarr(nx)
         gammaxH2_minus=dblarr(nx)
         for k=0,nx-1 do begin
            gammaxH2_plus(k)=vth*total(Vr2pidVr*(fH2(*,ip,k)#(Vx(ip)*dVx(ip))))
            gammaxH2_minus(k)=vth*total(Vr2pidVr*(fH2(*,in,k)#(Vx(in)*dVx(in))))
         endfor
         data=[gammaxH2_plus,gammaxH2_minus,gammaxH2]
         jp=where(data lt 1.0e32)
         yrange=[min(data(jp)),max(data(jp))]
         plot,x,gammaxH2,/nodata,yrange=yrange,title=_HH+' Fluxes',xtitle='x (meters)',ytitle='m!U-2!N s!U-1!N'
         oplot,x,gammaxH2,color=2
         xyouts,x(mid1),GammaxH2(mid1),'!7C!5',color=2
         oplot,x,gammaxH2_plus,color=3
         xyouts,x(mid2),GammaxH2_plus(mid2),'!7C!5!U+!N',color=3
         oplot,x,GammaxH2_minus,color=4
         xyouts,x(mid3),GammaxH2_minus(mid3),'!7C!5!U-!N',color=4
         if pause then press_return
      endif
;     
      Source_Error=fltarr(nx)
;
; Compute Source error
;     
      if compute_errors then begin
         if debrief gt 1 then print,prompt+'Computing Source Error'
;
; Test Mass Balance
;
; The relationship, 2 dGammaxH2/dx - 2 SH2 + SH + SP + 2 nHp x Nuloss = 0, should be satisfied.
;
         dGammaxH2dx=dblarr(nx-1)
         SH_p=dblarr(nx-1)
         for k=0,nx-2 do dGammaxH2dx(k)=(GammaxH2(k+1)-GammaxH2(k))/(x(k+1)-x(k))
         for k=0,nx-2 do SH_p(k)=0.5*(SH(k+1)+SP(k+1)+2*NuLoss(k+1)*nHp(k+1)-2*SH2(k+1) $
                                     +SH(k)+SP(k)+2*NuLoss(k)*nHp(k)-2*SH2(k))
         max_source=max([SH,2*SH2])
         for k=0,nx-2 do Source_Error(k)=abs(2*dGammaxH2dx(k)+SH_p(k))/max(abs([2*dGammaxH2dx(k),SH_p(k),max_source]))
         if debrief gt 0 then print,prompt+'Maximum Normalized Source_error =',max(Source_Error)
      endif
   endif
;
; Save input parameters in common block
;
   vx_s=vx
   vr_s=vr
   x_s=x
   Tnorm_s=Tnorm
   mu_s=mu
   Ti_s=Ti
   vxi_s=vxi
   Te_s=Te
   n_s=n
   vxi_s=vxi
   fH2BC_s=fH2BC
   GammaxH2BC_s=GammaxH2BC
   NuLoss_s=NuLoss
   PipeDia_s=PipeDia
   fH_s=fH
   SH2_s=SH2
   fH2_s=fH2
   nHP_s=nHP
   THP_s=THP
   Simple_CX_s=Simple_CX
   Sawada_s=Sawada
   H2_H2_EL_s=H2_H2_EL
   H2_P_EL_s=H2_P_EL
   H2_H_EL_s=H2_H_EL
   H2_HP_CX_s=H2_HP_CX
   ni_correct_s=ni_correct
;
; Set output parameters to single precision
; 
;   fH2=float(fH2)
;   nH2=float(nH2)
;   GammaxH2=float(GammaxH2)
;   VxH2=float(VxH2)
;   pH2=float(pH2)
;   TH2=float(TH2)
;   qxH2=float(qxH2)
;   qxH2_total=float(qxH2_total)
;   Sloss=float(Sloss)
;   QH2=float(QH2)
;   RxH2=float(RxH2)
;   QH2_total=float(QH2_total)
;   AlbedoH2=float(AlbedoH2)
;   nHP=float(nHP)
;   THP=float(THP)
;   fSH=float(fSH)
;   SH=float(SH)
;   SP=float(SP)
;   SHP=float(SHP)
;   NuE=float(NuE)
;   NuDis=float(NuDis)
;   piH2_xx=float(piH2_xx)
;   piH2_yy=float(piH2_yy)
;   piH2_zz=float(piH2_zz)
;   RxH_H2=float(RxH_H2)
;   RxHp_H2=float(RxHp_H2)
;   RxP_H2=float(RxP_H2)
;   RxP_H2_CX=float(RxP_H2_CX)
;   RxH_H2=float(RxH_H2)
;   EHp_H2=float(EHp_H2)
;   EP_H2=float(EP_H2)
;   EP_H2_CX=float(EP_H2_CX)
;   Epara_PerpH2_H2=float(Epara_PerpH2_H2)
;   ESH=float(ESH)
;   Eaxis=float(Eaxis)
;
;   vr=float(vr)
;   vx=float(vx)
;   x=float(x)

Return:
   if debug gt 0 then begin
      print,prompt+'Finished'
      press_return
   endif
   return
   end
