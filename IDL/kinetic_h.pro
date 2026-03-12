;+
; Kinetic_H.pro
;
; This subroutine is part of the "KN1D" atomic and molecular neutral transport code.
;
;   This subroutine solves a 1-D spatial, 2-D velocity kinetic neutral transport 
; problem for atomic hydrogen (H) or deuterium by computing successive generations of 
; charge exchange and elastic scattered neutrals. The routine handles electron-impact 
; ionization, proton-atom charge exchange, radiative recombination, and elastic
; collisions with hydrogenic ions, neutral atoms, and molecules.
;
;   The positive vx half of the atomic neutral distribution function is inputted at x(0) 
; (with arbitrary normalization) and the desired flux of hydrogen atoms entering the slab,
; at x(0) is specified. Background profiles of plasma ions, (e.g., Ti(x), Te(x), n(x), vxi(x),...)
; molecular ions, (nHP(x), THP(x)), and molecular distribution function (fH) are inputted.
;
; Optionally, the hydrogen source velocity distribution function is also inputted.
; (The H source and fH2 distribution functions can be computed using procedure 
; "Kinetic_H2.pro".) The code returns the atomic hydrogen distribution function, fH(vr,vx,x) 
; for all vx, vr, and x of the specified vr,vx,x grid.
;
;   Since the problem involves only the x spatial dimension, all distribution functions
; are assumed to have rotational symmetry about the vx axis. Consequently, the distributions
; only depend on x, vx and vr where vr =sqrt(vy^2+vz^2)
;
;  History:
;
;    B. LaBombard   First coding based on Kinetic_Neutrals.pro 		22-Dec-2000
;
;    For more information, see write-up: "A 1-D Space, 2-D Velocity, Kinetic 
;    Neutral Transport Algorithm for Hydrogen Atoms in an Ionizing Plasma", B. LaBombard
;
; Note: Variable names contain characters to help designate species -
;	atomic neutral (H), molecular neutral (H2), molecular ion (HP), proton (i) or (P) 
;
;________________________________________________________________________________
pro Kinetic_H,vx,vr,x,Tnorm,mu,Ti,Te,n,vxi,fHBC,GammaxHBC,PipeDia,fH2,fSH,nHP,THP,$
       fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,WallH,$
       truncate=truncate,Simple_CX=Simple_CX,Max_Gen=Max_Gen,$
       No_Johnson_Hinnov=No_Johnson_Hinnov,No_Recomb=No_Recomb,$
       H_H_EL=H_H_EL,H_P_EL=H_P_EL,H_H2_EL=_H_H2_EL,H_P_CX=H_P_CX,ni_correct=ni_correct,$
       error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause

   common Kinetic_H_Output,piH_xx,piH_yy,piH_zz,RxHCX,RxH2_H,RxP_H,RxW_H,EHCX,EH2_H,EP_H,EW_H,Epara_PerpH_H,SourceH,SRecomb

   common Kinetic_H_Errors,Max_dx,vbar_error,mesh_error,moment_error,C_Error,CX_Error,H_H_error,$
                           qxH_total_error,QH_total_error
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
;		   fHBC	- fltarr(nvr,nvx), this is an input boundary condition
;			  specifying the shape of the neutral atom velocity distribution 
;			  function at location x(0). Normalization is arbitrary.
;		          Only values with positive vx, fHBC(*,nvx/2:*) are used
;		          by the code.
;	      GammaxHBC	- float, desired neutral atom flux density in the +Vx
;			  direction at location x(0) (m^-2 s^-1)
;			  fHBC is scaled to yield this flux density.
;	        PipeDia	- fltarr(nx), effective pipe diameter (meters)
;			  This variable allows collisions with the 'side-walls' to be simulated.
;			  If this variable is undefined, then PipeDia set set to zero. Zero values
;			  of PipeDia are ignored (i.e., treated as an infinite diameter).
;                   fH2	- fltarr(nvr,nvx,nx), neutral molecule velocity distribution
;                         function. fH2 is normalized so that the molecular neutral density, nH2(k), is 
;			  defined as the velocity space integration: nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
;                         If this variable is undefined, then it is set equal to zero and
;                         no molecule-atom collisions are included.
;			  NOTE: dVx is velocity space differential for Vx axis and Vr2pidVr = Vr*!pi*dVr
;		                with dVr being velocity space differential for Vr axis.
;                   fSH	- fltarr(nvr,nvx,nx), atomic hydrogen source velocity distribution.
;                         fSH must be normalized so that the total atomic neutral
;                         source, SourceH(k), is defined as the velocity space integration:
;                             SourceH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
;			  fSH can be computed from IDL procedure Kinetic_H2.pro
;                         If this variable is undefined, then it is set equal to zero.
;                   nHP	- fltarr(nx), molecular ion density profile (m^-3)
;                         If this parameter is undefined, then it is set equal to zero.
;			  nHP can be computed from IDL procedure Kinetic_H2.pro
;                   THP	- fltarr(nx), molecular ion temperature profile (m^-3)
;                         If this parameter is undefined, then it is set equal to 3 eV at each grid point.
;			  THP can be computed from IDL procedure Kinetic_H2.pro
;
;  Input & Output:
;                    fH	- fltarr(nvr,nvx,nx), neutral atom velocity distribution
;                         function. 'Seed' values for this may be specified on input. 
;		          If this parameter is undefined on input, then a zero 'seed' value will be used. 
;			  The algorithm outputs a self-consistent fH.
;			  fH is normalized so that the neutral density, nH(k), is defined as 
;			  the velocity space integration: nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
;
;  Output:
;                    nH	- fltarr(nx), neutral atom density profile (m^-3)
;               GammaxH	- fltarr(nx), neutral atom flux profile (# m^-2 s^-1)
;                             computed from GammaxH(k)=Vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))
;                   VxH	- fltarr(nx), neutral atom velocity profile (m s^-1)
;                             computed from GammaxH/nH
;
;                       To aid in computing the some of the quantities below, the procedure internally
;                       defines the quantities:
;                       vr2vx2_ran(i,j,k)=vr(i)^2+(vx(j)-VxH(k))^2
;                                     which is the magnitude of 'random v^2' at each mesh point
;                       vr2vx2(i,j,k)=vr(i)^2+vx(j)^2
;                                     which is the magnitude of 'total v^2' at each mesh point
;                       q=1.602177D-19, mH=1.6726231D-27
;                       C(*,*,*) is the right hand side of the Boltzmann equation, evaluated
;                                using the computed neutral distribution function
;
;                    pH	- fltarr(nx), neutral atom pressure (eV m^-2) computed from:
;                         pH(k)~vth2*total(Vr2pidVr*(vr2vx2_ran(*,*,k)*fH(*,*,k))#dVx))*(mu*mH)/(3*q)
;                    TH	- fltarr(nx), neutral atom temperature profile (eV) computed from: TH=pH/nH
;                   qxH	- fltarr(nx), neutral atom random heat flux profile (watts m^-2) computed from:
;                          qxH(k)~vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH(*,*,k))#(dVx*(vx-VxH(k)))))*0.5*(mu*mH)
;             qxH_total	- fltarr(nx), total neutral atom heat flux profile (watts m^-2)
;                             This is the total heat flux transported by the neutrals:
;                         qxH_total=(0.5*nH*(mu*mH)*VxH*VxH + 2.5*pH*q)*VxH + piH_xx*VxH + qxH
;	     NetHSource	- fltarr(nx), net H0 source [H0 source - ionization sink - wall sink] (m^-3 s^-1) computed from
;				NetHSource(k)=total(Vr2pidVr*(C(*,*,k)#dVx))
;		   Sion	- fltarr(nx), H ionization rate (m^-3 s^-1) 
;                    QH	- fltarr(nx), rate of net thermal energy transfer into neutral atoms (watts m^-3) computed from
;                               QH(k)~vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C(*,*,k))#dVx))*0.5*(mu*mH)
;                   RxH	- fltarr(nx), rate of x momentum transfer to neutral atoms (=force, N m^-2).
;                               RxH(k)~Vth*total(Vr2pidVr*(C(*,*,k)#(dVx*(vx-VxH(k)))))*(mu*mH)
;              QH_total	- fltarr(nx), net rate of total energy transfer into neutral atoms
;                          = QH + RxH*VxH - 0.5*(mu*mH)*(Sloss-SourceH)*VxH*VxH (watts m^-3)
;               AlbedoH	- float, Ratio of atomic neutral particle flux with Vx < 0 divided by particle flux
;                          with Vx > 0  at x=x(0)
;                          (Note: For fSH non-zero, the flux with Vx < 0 will include
;                          contributions from molecular hydrogen sources within the 'slab'.
;                          In this case, this parameter does not return the true 'Albedo'.)
;	          WallH	- fltarr(nx), atomic neutral sink rate arising from hitting the 'side walls' (m^-3 s^-1)
;			   Unlike the molecules in Kinetic_H2, wall collisions result in the destruction of atoms.
;			   This parameter can be used to specify a resulting source of molecular
;			   neutrals in Kinetic_H2. (molecular source = 2 times WallH)
;
; KEYWORDS:
;   Output:
;	          error	- Returns error status: 0=no error, solution returned
;					        1=error, no solution returned
;
; COMMON BLOCK Kinetic_H_OUTPUT
;    Output:
;                piH_xx	- fltarr(nx), xx element of stress tensor (eV m^-2) computed from:
;                         piH_xx(k)~vth2*total(Vr2pidVr*(fH(*,*,k)#(dVx*(vx-VxH(k))^2)))*(mu*mH)/q - pH
;                piH_yy	- fltarr(nx), yy element of stress tensor (eV m^-2) computed from:
;                         piH_yy(k)~vth2*total((Vr2pidVr*Vr^2)*(fH(*,*,k)#dVx))*(mu*mH)/q - pH
;                piH_zz	- fltarr(nx), zz element of stress tensor (eV m^-2) = piH_yy
;			   Note: cylindrical system relates r^2 = y^2 + z^2. All other stress tensor elements are zero.
;
;           The following momentum and energy transfer rates are computed from charge-exchange collsions between species:
;                 RxHCX	- fltarr(nx), rate of x momentum transfer from hydrogren ions to atoms (=force/vol, N m^-3).
;                  EHCX	- fltarr(nx), rate of energy transfer from hydrogren ions to atoms (watts m^-3).
;               
;           The following momentum and energy transfer rates are computed from elastic collsions between species:
;                RxH2_H	- fltarr(nx), rate of x momentum transfer from neutral molecules to atoms (=force/vol, N m^-3).
;                 RxP_H	- fltarr(nx), rate of x momentum transfer from hydrogen ions to neutral atoms (=force/vol, N m^-3).
;                 EH2_H	- fltarr(nx), rate of energy transfer from neutral molecules to atoms (watts m^-3).
;                  EP_H	- fltarr(nx), rate of energy transfer from hydrogen ions to neutral atoms (watts m^-3).
;
;           The following momentum and energy transfer rates are computed from collisions with the 'side-walls'
;                 RxW_H	- fltarr(nx), rate of x momentum transfer from wall to neutral atoms (=force/vol, N m^-3).
;                  EW_H	- fltarr(nx), rate of energy transfer from wall to neutral atoms (watts m^-3).
;
;           The following is the rate of parallel to perpendicular energy transfer computed from elastic collisions
;         Epara_PerpH_H	- fltarr(nx), rate of parallel to perp energy transfer within atomic hydrogen species (watts m^-3).
;
;           Source/Sink info:
;               SourceH	- fltarr(nx), source rate of neutral atoms from H2 dissociation (from integral of inputted fSH) (m^-3 s^-1).
;                SRecom	- fltarr(nx), source rate of neutral atoms from recombination (m^-3 s^-1).
;
; KEYWORDS:
;   Input:
;	       truncate	- float, stop computation when the maximum 
;			  increment of neutral density normalized to 
;			  inputed neutral density is less than this 
;	    		  value in a subsequent generation. Default value is 1.0e-4
;
;             Simple_CX	- if set, then use CX source option (B): Neutrals are born
;                         in velocity with a distribution proportional to the local
;                         ion distribution function. Simple_CX=1 is default.
;
;                         if not set, then use CX source option (A): The CX source
;                         neutral atom distribution function is computed by evaluating the
;                         the CX cross section for each combination of (vr,vx,vr',vx')
;                         and convolving it with the neutral atom distribution function.
;                         This option requires more CPU time and memory.
;
;      	  	Max_gen	- integer, maximum number of collision generations to try including before giving up.
;                         Default is 50.
;
;     No_Johnson_Hinnov	- if set, then compute ionization and recombination rates
;			  directly from reaction rates published by Janev* for
;			  ground state hydrogen
;
;			      Ionization:    e + H(1s) -> p + e 
;			      Recombination: e + p -> H(1s) + hv
;
;			  *Janev, R.K., et al, "Elementary processes in hydrogen-helium plasmas",
;			   (Springer-Verlag, Berlin ; New York, 1987)
;
;			  Otherwise, compute ionization and recombination rates using
;		          results from the collisional-radiative model published by Johnson
;			  and Hinnov [L.C.Johnson and E. Hinnov, J. Quant. Spectrosc. Radiat.
; 			  Transfer. vol. 13 pp.333-358]. This is the default.
;			  Note: charge exchange is always computed using the ground state reaction
;		          rates published by Janev:
;
;			      Charge Exchange: p + H(1s) -> H(1s) + p
;			  
;	      No_Recomb	- if set, then DO NOT include recombination as a source of atomic neutrals
;		          in the algorithm
;
;	 	 H_H_EL	- if set, then include H -> H elastic self collisions
;			     Note: if H_H_EL is set, then algorithm iterates fH until
;	                     self consistent fH is achieved.
;	 	 H_P_CX	- if set, then include H -> H(+) charge exchange collisions 
;	         H_P_EL	- if set, then include H -> H(+) elastic collisions 
;	        H_H2_EL	- if set, then include H -> H2 elastic collisions 
;            ni_correct	- if set, then algorithm corrects hydrogen ion density
;			     according to quasineutrality: ni=ne-nHP. Otherwise, nHP is assumed to be small.
;
;	 Compute_Errors	- if set, then return error estimates in common block Kinetic_H_ERRORS below
;
;		   plot	- 0= no plots, 1=summary plots, 2=detail plots, 3=very detailed plots
;		  debug	- 0= do not execute debug code, 1=summary debug, 2=detail debug, 3=very detailed debug
;	        debrief	- 0= do not print, 1=print summary information, 2=print detailed information
;	          pause	- if set, then pause between plots
;
; COMMON BLOCK Kinetic_H_ERRORS
;
;	if COMPUTE_ERRORS keyword is set then the following is returned in common block Kinetic_H_ERRORS
;
;	         Max_dx	- float(nx), Max_dx(k) for k=0:nx-2 returns maximum 
;			  allowed x(k+1)-x(k) that avoids unphysical negative 
;			  contributions to fH
;	     Vbar_error	- float(nx), returns numerical error in computing
;			  the speed of ions averged over maxwellian distribution.
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
;	      H_H_error	-  fltarr(nx,[0,1,2]) return normalized errors associated with 
;		           particle [0], x-momentum [1], and total energy [2] convervation of the elastic self-collision operator
;
;       qxH_total_error	- fltarr(nx), normalized error estimate in computation of qxH_total
;        QH_total_error	- fltarr(nx), normalized error estimate in computation of QH_total
;
; History:
;	22-Dec-2000 - B. LaBombard - first coding.
;	11-Feb-2001 - B. LaBombard - added elastic collisions 
;
;______________________________________________________________________
;-
   prompt='Kinetic_H => '
;
   common Kinetic_H_input,vx_s,vr_s,x_s,Tnorm_s,mu_s,Ti_s,Te_s,n_s,vxi_s,fHBC_s,GammaxHBC_s,PipeDia_s,fH2_s,fSH_s,nHP_s,THP_s,$
                          fH_s,Simple_CX_s,JH_s,Recomb_s,H_H_EL_s,H_P_EL_s,H_H2_EL_s,H_P_CX_s

   common Kinetic_H_internal,vr2vx2,vr2vx_vxi2,fi_hat,ErelH_P,Ti_mu,ni,sigv,alpha_ion,v_v2,v_v,vr2_vx2,vx_vx,$
          Vr2pidVrdVx,SIG_CX,SIG_H_H,SIG_H_H2,SIG_H_P,Alpha_CX,Alpha_H_H2,Alpha_H_P,MH_H_sum,Delta_nHs,Sn,Rec
   common Kinetic_H_H2_Moments,nH2,VxH2,TH2
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
   key_default,plot,0
   key_default,debug,0
   key_default,pause,0
   key_default,debrief,0
   if debug gt 0 then plot=plot > 1
   if debug gt 0 then debrief=debrief > 1
   if debug gt 0 then pause=1
   key_default,Simple_CX,1
   key_default,Max_Gen,50
   key_default,No_Johnson_Hinnov,0
   JH=1 & if No_Johnson_Hinnov then JH=0
   key_default,No_Recomb,0
   Recomb=1 & if No_Recomb then Recomb=0

   key_default,H_H_EL,0
   key_default,H_P_EL,0
   key_default,_H_H2_EL,0
   key_default,H_P_CX,0
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
   if n_elements(vxi) ne nx then begin & print,prompt+'Number of elements in vxi and x do not agree!' & error=1 & goto,return & endif
   if n_elements(Te) ne nx then begin & print,prompt+'Number of elements in Te and x do not agree!' & error=1 & goto,return & endif
   if n_elements(n) ne nx then begin & print,prompt+'Number of elements in n and x do not agree!' & error=1 & goto,return & endif
   if type_of(GammaxHBC) eq 0 then  begin & print,prompt+'GammaxHBC is not defined!' & error=1 & goto,return & endif
   if type_of(PipeDia) eq 0 then PipeDia=dblarr(nx)

   if n_elements(PipeDia) ne nx then begin & print,prompt+'Number of elements in PipeDia and x do not agree!' & error=1 & goto,return & endif
   if n_elements(fHBC(*,0)) ne nvr then begin & print,prompt+'Number of elements in fHBC(*,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fHBC(0,*)) ne nvx then begin & print,prompt+'Number of elements in fHBC(0,*) and vx do not agree!' & error=1 & goto,return & endif
   if type_of(fH2) eq 0 then fH2=dblarr(nvr,nvx,nx)
   ; fH2=dblarr(nvr,nvx,nx)
   if n_elements(fH2(*,0,0)) ne nvr then begin & print,prompt+'Number of elements in fH2(*,0,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2(0,*,0)) ne nvx then begin & print,prompt+'Number of elements in fH2(0,*,0) and vx do not agree!' & error=1 & goto,return & endif
   if n_elements(fH2(0,0,*)) ne nx then begin & print,prompt+'Number of elements in fH2(0,0,*) and x do not agree!' & error=1 & goto,return & endif
   if type_of(fSH) eq 0 then fSH=dblarr(nvr,nvx,nx)
   ;fSH=dblarr(nvr,nvx,nx)
   if n_elements(fSH(*,0,0)) ne nvr then begin & print,prompt+'Number of elements in fSH(*,0,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fSH(0,*,0)) ne nvx then begin & print,prompt+'Number of elements in fSH(0,*,0) and vx do not agree!' & error=1 & goto,return & endif
   if n_elements(fSH(0,0,*)) ne nx then begin & print,prompt+'Number of elements in fSH(0,0,*) and x do not agree!' & error=1 & goto,return & endif
   if type_of(nHP) eq 0 then nHP=dblarr(nx)
   ;nHP=dblarr(nx)
   if n_elements(nHP) ne nx then begin & print,prompt+'Number of elements in nHP and x do not agree!' & error=1 & goto,return & endif
   if type_of(THP) eq 0 then THP=dblarr(nx)+1.0
   ;THP=dblarr(nx)+1.0
   if n_elements(THP) ne nx then begin & print,prompt+'Number of elements in THP and x do not agree!' & error=1 & goto,return & endif
   if type_of(fH) eq 0 then fH=dblarr(nvr,nvx,nx)
   ;fH=dblarr(nvr,nvx,nx)
   if n_elements(fH(*,0,0)) ne nvr then begin & print,prompt+'Number of elements in fH(*,0,0) and vr do not agree!' & error=1 & goto,return & endif
   if n_elements(fH(0,*,0)) ne nvx then begin & print,prompt+'Number of elements in fH(0,*,0) and vx do not agree!' & error=1 & goto,return & endif
   if n_elements(fH(0,0,*)) ne nx then begin & print,prompt+'Number of elements in fH(0,0,*) and x do not agree!' & error=1 & goto,return & endif
   if total(abs(vr)) le 0.0 then begin & print,prompt+'vr is all 0!' & error=1  & goto,return & endif
   ii=where(vr le 0,count)
   if count gt 0 then begin & print,prompt+'vr contains zero or negative element(s)!' & error=1 & goto,return & endif
   if total(abs(vx)) le 0.0 then begin & print,prompt+'vx is all 0!' & error=1 & goto,return & endif
   if total(x) le 0.0 then begin & print,prompt+'Total(x) is less than or equal to 0!' & error=1 & goto,return & endif
   if type_of(Tnorm) eq 0 then begin & print,prompt+'Tnorm is not defined!' & error=1 & goto,return & endif
   if type_of(mu) eq 0 then begin & print,prompt+'mu is not defined!' & error=1 & goto,return & endif
   if mu ne 1 and mu ne 2 then begin & print,prompt+'mu must be 1 or 2!' & error=1 & goto,return & endif
   _e='e!U-!N'
   _hv='hv'
   if mu eq 1 then begin
      _p='H!U+!N'
      _H='H!U0!N'
      _H1s='H(1s)'
      _Hs='H!U*!N(2s)'
      _Hp='H!U*!N(2p)'
      _Hn2='H!U*!N(n=2)'
      _Hn3='H!U*!N(n=3)'
      _Hn='H!U*!N(n>=2)'
      _HH='H!D2!N'
      _Hp='H!D2!U+!N'
   endif else begin
      _p='D!U+!N'
      _H='D!U0!N'
      _H1s='D(1s)'
      _Hs='D!U*!N(2s)'
      _Hp='D!U*!N(2p)'
      _Hn2='D!U*!N(n=2)'
      _Hn3='D!U*!N(n=3)'
      _Hn='D!U*!N(n>=2)'
      _HH='D!D2!N'
      _Hp='D!D2!U+!N'
   endelse
   plus=' + '
   arrow=' -> '
   elastic=' (elastic)'
   _R1=_e+plus+_H1s+arrow+_e+plus+_p+plus+_e
   _R2=_e+plus+_p+arrow+plus+_H1s+plus+_hv
   _R3=_p+plus+_H1s+arrow+_H1s+plus+_p
   _R4=_H+plus+_p+arrow+_H+plus+_p+elastic
   _R5=_H+plus+_HH+arrow+_H+plus+_HH+elastic
   _R6=_H+plus+_H+arrow+_H+plus+_H+elastic
   _Rn=[' ',_R1,_R2,_R3,_R4,_R5,_R6]
   in=where(vx lt 0,count)
   if count lt 1 then begin & print,prompt+'vx contains no negative elements!' & error=1 & goto,return & endif
   ip=where(vx gt 0,count)
   if count lt 1 then begin & print,prompt+'vx contains no positive elements!' & error=1 & goto,return & endif
   iz=where(vx eq 0,count)
   if count gt 0 then begin & print,prompt+'vx contains one or more zero elements!' & error=1 & goto,return & endif
   diff=where(vx(ip) ne -reverse(vx(in)),count)
   if count gt 0 then begin & print,prompt+'vx array elements are not symmetric about zero!' & error=1 & goto,return & endif
   fHBC_input=fHBC
   fHBC_input(*)=0.0
   fHBC_input(*,ip)=fHBC(*,ip)
   test=total(fHBC_input)
   if test le 0.0 and abs(GammaxHBC) gt 0 then begin & print,prompt+'Values for fHBC(*,*) with vx > 0 are all zero!' & error=1 & goto,return & endif
;
; Output variables
;
   nH=dblarr(nx)
   GammaxH=dblarr(nx)
   VxH=dblarr(nx)
   pH=dblarr(nx)
   TH=dblarr(nx)
   qxH=dblarr(nx)
   qxH_total=dblarr(nx)
   NetHSource=dblarr(nx)
   WallH=dblarr(nx)
   Sion=dblarr(nx)
   QH=dblarr(nx)
   RxH=dblarr(nx)
   QH_total=dblarr(nx)
   piH_xx=dblarr(nx)
   piH_yy=dblarr(nx)
   piH_zz=dblarr(nx)
   RxHCX=dblarr(nx)
   RxH2_H=dblarr(nx)
   RxP_H=dblarr(nx)
   RxW_H=dblarr(nx)
   EHCX=dblarr(nx)
   EH2_H=dblarr(nx)
   EP_H=dblarr(nx)
   EW_H=dblarr(nx)
   Epara_PerpH_H=dblarr(nx)
   AlbedoH=0.0D0
   SourceH=dblarr(nx)
   SRecomb=dblarr(nx)
;
; Internal variables
;
   mH=1.6726231D-27
   q=1.602177D-19				
   k_boltz=1.380658D-23                 &;Bolzmann's constant, J K^-1
   Twall=293.0*k_boltz/q                &;room temperature (eV)
   Work=dblarr(nvr*nvx)
   fHG=dblarr(nvr,nvx,nx)
   NHG=dblarr(nx,max_gen+1)
   Vth=sqrt(2*q*Tnorm/(mu*mH))
   Vth2=vth*vth
   Vth3=Vth2*Vth
   fHs=dblarr(nx)
   nHs=dblarr(nx)
   Alpha_H_H=dblarr(nvr,nvx)
   Omega_H_P=dblarr(nx)
   Omega_H_H2=dblarr(nx)
   Omega_H_H=dblarr(nx)
   VxHG=dblarr(nx)
   THG=dblarr(nx)
   Wperp_paraH=dblarr(nx)
   vr2vx2_ran2=dblarr(nvr,nvx)
   vr2_2vx_ran2=dblarr(nvr,nvx)
   vr2_2vx2_2D=dblarr(nvr,nvx)
   RxCI_CX=dblarr(nx)
   RxCI_H2_H=dblarr(nx)
   RxCI_P_H=dblarr(nx)
   Epara_Perp_CI=dblarr(nx)
   CI_CX_error=fltarr(nx)
   CI_H2_H_error=fltarr(nx)
   CI_P_H_error=fltarr(nx)
   CI_H_H_error=fltarr(nx)
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
; Scale input molecular distribution function to agree with desired flux
;
   gamma_input=1.0
   if abs(GammaxHBC) gt 0 then gamma_input=vth*total(Vr2pidVr*(fHBC_input#(Vx*dVx)))
   ratio=abs(GammaxHBC)/gamma_input
   fHBC_input=fHBC_input*ratio
   if abs(ratio-1) gt 0.01*truncate then begin
      fHBC=fHBC_input
   endif
   print, ip
   fH(*,ip,0)=fHBC_input(*,ip)


;
; if fH2 is zero, then turn off elastic H2 <-> H collisions
;
   H_H2_EL=_H_H2_EL
   if total(fH2) le 0.0 then H_H2_EL=0
;
; Set iteration scheme
;
   fH_iterate=0
   if (H_H_EL ne 0) or (H_P_EL ne 0) or (H_H2_EL ne 0) then fH_iterate=1

   fH_generations=0
   if (fH_iterate ne 0) or (H_P_CX ne 0) then fH_generations=1

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
   New_Molecular_Ions=1
   if type_of(nHP_s) ne 0 then begin
      test=0
      ii=where(nHP_s ne nHP,count) & test=test+count
      ii=where(THP_s ne THP,count) & test=test+count
      if test le 0 then New_Molecular_Ions=0
   endif
   New_Electrons=1
   if type_of(Te_s) ne 0 then begin
      test=0
      ii=where(Te_s ne Te,count) & test=test+count
      ii=where(n_s ne n,count) & test=test+count
      if test le 0 then New_Electrons=0
   endif
   New_fH2=1
   if type_of(fH2_s) ne 0 then begin
      ii=where(fH2_s ne fH2,count)
      if count le 0 then New_fH2=0
   endif
   New_fSH=1
   if type_of(fSH_s) ne 0 then begin
      ii=where(fSH_s ne fSH,count)
      if count le 0 then New_fSH=0
   endif
   New_Simple_CX=1
   if type_of(Simple_CX_s) ne 0 then begin
      ii=where(Simple_CX_s ne Simple_CX,count)
      if count le 0 then New_Simple_CX=0
   endif
   New_H_Seed=1
   if type_of(fH_s) ne 0 then begin
      ii=where(fH_s ne fH,count)
      if count le 0 then New_H_Seed=0
   endif


   Do_sigv=        New_Grid or New_Electrons
   Do_ni=          New_Grid or New_Electrons or New_Protons or New_Molecular_Ions
   Do_fH2_moments=(New_Grid or New_fH2) and total(fH2) gt 0.0
   Do_Alpha_CX=   (New_Grid or (type_of(Alpha_CX) eq 0) or Do_ni or New_Simple_CX)        and H_P_CX
   Do_SIG_CX=     (New_Grid or (type_of(SIG_CX) eq 0) or New_Simple_CX)  and (Simple_CX eq 0) and Do_Alpha_CX
   Do_Alpha_H_H2= (New_Grid or (type_of(Alpha_H_H2) eq 0) or New_fH2)    and H_H2_EL
   Do_SIG_H_H2=   (New_Grid or (type_of(SIG_H_H2) eq 0))                 and Do_Alpha_H_H2
   Do_SIG_H_H=    (New_Grid or (type_of(SIG_H_H) eq 0))                  and H_H_EL
   Do_Alpha_H_P=  (New_Grid or (type_of(Alpha_H_P) eq 0) or Do_ni) 	 and H_P_EL 
   Do_SIG_H_P=   (New_Grid or (type_of(SIG_H_P) eq 0))                 and Do_Alpha_H_P
   Do_v_v2=      (New_Grid or (type_of(v_v2) eq 0))                    and (CI_Test or Do_SIG_CX or Do_SIG_H_H2 or Do_SIG_H_H or Do_SIG_H_P)

   nH2=dblarr(nx)
   vxH2=dblarr(nx)
   TH2=dblarr(nx)+1.0
   if Do_fH2_moments then begin
      if debrief gt 1 then print,prompt+'Computing vx and T moments of fH2'
;
; Compute x flow velocity and temperature of molecular species
;
      for k=0,nx-1 do begin
         nH2(k)=total(Vr2pidVr*(fH2(*,*,k)#dVx))
         if nH2(k) gt 0 then begin
            vxH2(k)=vth*total(Vr2pidVr*(fH2(*,*,k)#(Vx*dVx)))/nH2(k)
            for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-vxH2(k)/vth)^2
            TH2(k)=(2*mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fH2(*,*,k))#dVx))/(3*q*nH2(k))
         endif
      endfor
   endif

   if New_Grid then begin
      if debrief gt 1 then print,prompt+'Computing vr2vx2, vr2vx_vxi2, ErelH_P'
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
; Atomic hydrogen ion energy in local rest frame of plasma at each mesh point
;
      ErelH_P=0.5*mH*vr2vx_vxi2*vth2/q
      ErelH_P=ErelH_P > 0.1    &; sigmav_cx does not handle neutral energies below 0.1 eV
      ErelH_P=ErelH_P < 2.0E4  &; sigmav_cx does not handle neutral energies above 20 keV
   endif
;
   if New_Protons then begin
      if debrief gt 1 then print,prompt+'Computing Ti/mu at each mesh point'
;
; Ti/mu at each mesh point
;
      Ti_mu=dblarr(nvr,nvx,nx)
      for k=0,nx-1 do Ti_mu(*,*,k)=Ti(k)/mu
;
; Compute Fi_hat
;
      if debrief gt 1 then print,prompt+'Computing fi_Hat'
      vx_shift=vxi
      Tmaxwell=Ti
      mol=1
@create_shifted_maxwellian.include
;      Create_Shifted_Maxwellian_Core, vr, vx, vx_shift, Tmaxwell, Maxwell, $
;         vth, Tnorm, Vr2pidVr, dVx, vol, vth_Dvx,   $
;         vx_Dvx, vr_Dvr, vr2vx2, vr2vx2_2D,          $
;         jpa, jpb, jna, jnb, mol, mu, mH, q, shifted_Maxwellian_debug
      
      fi_hat=Maxwell

      SAVE, fi_hat, FILENAME='fi_hat_NEW.dat'
   endif

   if Compute_errors then begin
      if debrief gt 1 then print,prompt+'Computing Vbar_Error'
;
; Test: The average speed of a non-shifted maxwellian should be 2*Vth*sqrt(Ti(x)/Tnorm)/sqrt(!pi)
;
      vx_shift=dblarr(nx)
      Tmaxwell=Ti
      mol=1
@create_shifted_maxwellian.include

      vbar_test=dblarr(nvr,nvx,ntheta)
      Vbar_Error=dblarr(nx)
      for m=0,ntheta-1 do vbar_test(*,*,m)=vr2vx2(*,*,0)
      _vbar_test=dblarr(nvr*nvx,ntheta)
      _vbar_test(*)=vth*sqrt(vbar_test)
      vbar_test=dblarr(nvr,nvx)
      vbar_test(*)=_vbar_test#dtheta
      for k=0,nx-1 do begin
         vbar=total(Vr2pidVr*((vbar_test*Maxwell(*,*,k))#dVx))
         vbar_exact=2*Vth*sqrt(Ti(k)/Tnorm)/sqrt(!pi) 
         Vbar_Error(k)=abs(vbar-vbar_exact)/vbar_exact
      endfor
      if debrief gt 0 then begin
         print,prompt+'Maximum Vbar error =',max(Vbar_Error)
      endif
   endif

   if Do_ni then begin
      if debrief gt 1 then print,prompt+'Computing ni profile'
      ni=n
      if ni_correct then ni=n-nHP
      ni=ni > 0.01*n
   endif


   if Do_sigv then begin
      if debrief gt 1 then print,prompt+'Computing sigv'
;
; Compute sigmav rates for each reaction with option to use rates
; from CR model of Johnson-Hinnov
;
      sigv=dblarr(nx,3)
;________________________________________________________________________________
; Reaction R1:  e + H -> e + H(+) + e 
;________________________________________________________________________________
      if JH then begin
         sigv(*,1)=JHS_Coef(n,Te,/no_null)
      endif else begin
         sigv(*,1)=sigmav_ion_h0(Te)
      endelse
;________________________________________________________________________________
; Reaction R2:  e + H(+) -> H(1s) + hv
;________________________________________________________________________________
      if JH then begin
         sigv(*,2)=JHAlpha_Coef(n,Te,/no_null)
      endif else begin
         sigv(*,2)=SigmaV_rec_H1s(Te)
      endelse
;
; H ionization rate (normalized by vth) = reaction 1
;
      alpha_ion=n*sigv(*,1)/vth
;
; Recombination rate (normalized by vth) = reaction 2
;
      Rec=n*sigv(*,2)/vth

   endif
;
; Compute Total Atomic Hydrogen Source
;
   Sn=dblarr(nvr,nvx,nx)
;
; Add Recombination (optionally) and User-Supplied Hydrogen Source (velocity space distribution)
;
   for k=0,nx-1 do begin
      Sn(*,*,k)=fSH(*,*,k)/vth
      if Recomb then Sn(*,*,k)=Sn(*,*,k)+fi_hat(*,*,k)*ni(k)*Rec(k)
   endfor
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
; vr2_vx2=0.125* [ vr2 + vr2_prime - 2*vr*vr_prime*cos(theta) - 2*(vx-vx_prime)^2 ]
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
        _Sig(*)=v_v*sigma_cx_h0(v_v2*(0.5*mH*vth2/q))
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
   if Do_SIG_H_H eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H_H'
;________________________________________________________________________________
; Compute SIG_H_H for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H_H * vr2_vx2 * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=vr2_vx2*v_v*sigma_EL_H_H(v_v2*(0.5*mH*mu*vth2/q),/VIS)/8.0
;
; Note: For viscosity, the cross section for D -> D is the same function of
;       center of mass energy as H -> H.
;
; Set SIG_H_H = vr' x Integral{vr2_vx2*v_v*sigma_H_H} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_H_H=dblarr(nvr*nvx,nvr*nvx)
        SIG_H_H(*)=Vr2pidVrdVx*(_Sig#dtheta)

;
; SIG_H_H is now vr' * sigma_H_H(v_v) * vr2_vx2 * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
   if Do_SIG_H_H2 eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H_H2'
;________________________________________________________________________________
; Compute SIG_H_H2 for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H_H2 * v_v at all possible relative velocities
;
        _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
        _Sig(*)=v_v*sigma_EL_H_HH(v_v2*(0.5*mH*vth2/q))
;
; NOTE: using H energy here for cross-sections tabulated as H->H2
;
; Set SIG_H_H2 = vr' x vx_vx x Integral{v_v*sigma_H_H2} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
        SIG_H_H2=dblarr(nvr*nvx,nvr*nvx)
        SIG_H_H2(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)
;
; SIG_H_H2 is now vr' *vx_vx * sigma_H_H2(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;
   if Do_SIG_H_P eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing SIG_H_P'
;________________________________________________________________________________
; Compute SIG_H_P for present velocity space grid, if it is needed and has not 
; already been computed with the present input parameters
;________________________________________________________________________________
;
; Compute sigma_H_P * v_v at all possible relative velocities
;
      _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
      _Sig(*)=v_v*sigma_EL_P_H(v_v2*(0.5*mH*vth2/q))

;
; Set SIG_H_P = vr' x vx_vx x Integral{v_v*sigma_H_P} over theta=0,2pi times differential velocity space element Vr'2pidVr'*dVx'
;
      SIG_H_P=dblarr(nvr*nvx,nvr*nvx)
      SIG_H_P(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)


;
; SIG_H_P is now vr' *vx_vx * sigma_H_P(v_v) * v_v (intergated over theta) for all possible ([vr,vx],[vr',vx'])
;
   endif
;________________________________________________________________________________ 
; Compute Alpha_CX for present Ti and ni, if it is needed and has not
; already been computed with the present parameters
;________________________________________________________________________________ 
   if Do_Alpha_CX eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_CX'

      if Simple_CX then begin
;________________________________________________________________________________
; Option (B): Use maxwellian weighted <sigma v>
;________________________________________________________________________________
;
; Charge Exchange sink rate
;
         alpha_cx=sigmav_cx_H0(Ti_mu,ErelH_P)/vth
         for k=0,nx-1 do alpha_cx(*,*,k)=alpha_cx(*,*,k)*ni(k)

;________________________________________________________________________________
;
      endif else begin
;________________________________________________________________________________
; Option (A): Compute SigmaV_CX from sigma directly via SIG_CX
;________________________________________________________________________________
;
         alpha_cx=dblarr(nvr,nvx,nx)
         for k=0,nx-1 do begin
            Work(*)=fi_hat(*,*,k)*ni(k)
            alpha_cx(*,*,k)=SIG_CX#Work
         endfor
         if do_alpha_cx_test then begin
            alpha_cx_test=sigmav_cx_H0(Ti_mu,ErelH_P)/vth
            for k=0,nx-1 do alpha_cx_test(*,*,k)=alpha_cx_test(*,*,k)*ni(k)
            print,'Compare alpha_cx and alpha_cx_test'
            press_return
         endif
      endelse
   endif
;________________________________________________________________________________ 
; Compute Alpha_H_H2 for inputted fH, if it is needed and has not
; already been computed with the present input parameters
;________________________________________________________________________________ 
   if Do_Alpha_H_H2 eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_H_H2'
      Alpha_H_H2=dblarr(nvr,nvx,nx)
      for k=0,nx-1 do begin
         Work(*)=fH2(*,*,k)
         Alpha_H_H2(*,*,k)=SIG_H_H2#Work
      endfor
   endif
;________________________________________________________________________________ 
; Compute Alpha_H_P for present Ti and ni 
; if it is needed and has not already been computed with the present parameters
;________________________________________________________________________________ 
   if Do_Alpha_H_P eq 1 then begin
      if debrief gt 1 then print,prompt+'Computing Alpha_H_P'
      Alpha_H_P=dblarr(nvr,nvx,nx)
      for k=0,nx-1 do begin
         Work(*)=fi_hat(*,*,k)*ni(k)
         HELP, Work
         HELP, SIG_H_P
         Alpha_H_P(*,*,k)=SIG_H_P#Work
      endfor
   endif

;
;________________________________________________________________________________
; Compute nH
;________________________________________________________________________________
   for k=0,nx-1 do nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))

   if New_H_Seed then begin
      MH_H_sum=dblarr(nvr,nvx,nx)
      Delta_nHs=1.0
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
fH_Iterate:
;
;  This is the entry point for fH iteration.
;  Save 'seed' values for comparison later
;
   fHs=fH
   nHs=nH
;________________________________________________________________________________
; Compute Omega values if nH is non-zero
;________________________________________________________________________________
   ii=where(nH le 0,count)
   if count le 0 then begin
;
; Compute VxH
;
      if H_P_EL or H_H2_EL or H_H_EL then for k=0,nx-1 do VxH(k)=vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))/nH(k)
;________________________________________________________________________________
; Compute Omega_H_P for present fH and Alpha_H_P if H_P elastic collisions are included
;________________________________________________________________________________
      print, 'H_P_EL'
      print, H_P_EL
      if H_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H_P'
         for k=0,nx-1 do begin
            DeltaVx=(VxH(k)-vxi(k))/vth
            MagDeltaVx=abs(DeltaVx) > DeltaVx_tol
            DeltaVx=sign(DeltaVx)*MagDeltaVx
            Omega_H_P(k)=total(Vr2pidVr*((Alpha_H_P(*,*,k)*fH(*,*,k))#dVx))/(nH(k)*DeltaVx)
         endfor
         Omega_H_P=Omega_H_P > 0.0
      endif
;________________________________________________________________________________
; Compute Omega_H_H2 for present fH and Alpha_H_H2 if H_H2 elastic collisions are included
;________________________________________________________________________________
      if H_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H_H2'
         for k=0,nx-1 do begin
            DeltaVx=(VxH(k)-vxH2(k))/vth
            MagDeltaVx=abs(DeltaVx) > DeltaVx_tol
            DeltaVx=sign(DeltaVx)*MagDeltaVx
            Omega_H_H2(k)=total(Vr2pidVr*((Alpha_H_H2(*,*,k)*fH(*,*,k))#dVx))/(nH(k)*DeltaVx)
         endfor
         Omega_H_H2=Omega_H_H2 > 0.0
      endif
;________________________________________________________________________________
; Compute Omega_H_H for present fH if H_H elastic collisions are included
;________________________________________________________________________________
      if H_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing Omega_H_H'
         if total(MH_H_sum) le 0 then begin
            for k=0,nx-1 do begin
               for i=0,nvr-1 do vr2_2vx_ran2(i,*)=vr(i)^2-2*(vx-VxH(k)/vth)^2
               Wperp_paraH(k)=total(Vr2pidVr*((vr2_2vx_ran2*fH(*,*,k))#dVx))/nH(k)
            endfor
         endif else begin
            for k=0,nx-1 do begin
               M_fH=MH_H_sum(*,*,k)-fH(*,*,k)
               Wperp_paraH(k)=-total(Vr2pidVr*((vr2_2vx2_2D*M_fH)#dVx))/nH(k)
            endfor
         endelse
         for k=0,nx-1 do begin
            Work(*)=fH(*,*,k)
            Alpha_H_H(*)=SIG_H_H#Work
            Wpp=Wperp_paraH(k)
            MagWpp=abs(Wpp) > Wpp_tol
            Wpp=sign(Wpp)*MagWpp
            Omega_H_H(k)=total(Vr2pidVr*((Alpha_H_H*Work)#dVx))/(nH(k)*Wpp)
         endfor
         Omega_H_H=Omega_H_H > 0.0
      endif
   endif
;
; Total Elastic scattering frequency
;
   Omega_EL=Omega_H_P+Omega_H_H2+Omega_H_H
;
; Total collision frequency
;
   alpha_c=dblarr(nvr,nvx,nx)
   if H_P_CX then begin
      for k=0,nx-1 do alpha_c(*,*,k)=alpha_cx(*,*,k)+alpha_ion(k)+Omega_EL(k)+gamma_wall(*,*,k)
   endif else begin
      for k=0,nx-1 do alpha_c(*,*,k)=alpha_ion(k)+Omega_EL(k)+gamma_wall(*,*,k)
   endelse
;   
; Test x grid spacing based on Eq.(27) in notes
;
   if debrief gt 1 then print,prompt+'Testing x grid spacing'
   Max_dx=fltarr(nx)  & Max_dx(*)=1.0e32
   for k=0,nx-1 do begin
      for j=ip(0),nvx-1 do begin 
         Max_dx(k)= Max_dx(k) < min(2*vx(j)/alpha_c(*,j,k))
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
         Fk(*,j,k)=(x(k+1)-x(k))*(Sn(*,j,k+1)+Sn(*,j,k))/denom
      endfor
   endfor
   for k=1,nx-1 do begin
      for j=0,ip(0)-1 do begin
         denom=-2*vx(j)+(x(k)-x(k-1))*alpha_c(*,j,k-1)
         Ck(*,j,k)=(-2*vx(j)-(x(k)-x(k-1))*alpha_c(*,j,k))/denom
         Dk(*,j,k)=(x(k)-x(k-1))/denom
         Gk(*,j,k)=(x(k)-x(k-1))*(Sn(*,j,k)+Sn(*,j,k-1))/denom
      endfor
   endfor
;
; Compute first-flight (0th generation) neutral distribution function
;
   Beta_CX_sum=dblarr(nvr,nvx,nx)
   MH_P_sum=dblarr(nvr,nvx,nx)
   MH_H2_sum=dblarr(nvr,nvx,nx)
   MH_H_sum=dblarr(nvr,nvx,nx)
   igen=0
   if debrief gt 0 then print,prompt+'Computing atomic neutral generation#'+sval(igen)

   fHG(*,ip,0)=fH(*,ip,0)
   for k=0,nx-2 do fHG(*,ip,k+1)=fHG(*,ip,k)*Ak(*,ip,k)+Fk(*,ip,k)
   for k=nx-1,1,-1 do fHG(*,in,k-1)=fHG(*,in,k)*Ck(*,in,k)+Gk(*,in,k)
;
; Compute first-flight neutral density profile
;
   for k=0,nx-1 do NHG(k,igen)=total(Vr2pidVr*(fHG(*,*,k)#dVx))

;
   if plot gt 1 then begin
      fH1d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH1d(*,k)=Vr2pidVr#fHG(*,*,k)
      endfor
      plot,vx,fH1d(*,0),/nodata,yrange=[0,max(fH1d)],title='First Generation '+_H
      for i=0,nx-1 do oplot,vx,fH1d(*,i),color=(i mod 8)+2
      if debug gt 0 then press_return
   endif
;
; Set total atomic neutral distribution function to first flight generation
;
   fH=fHG
   nH=NHG(*,0)

;
   if fH_generations eq 0 then goto,fH_done

next_generation:
   if igen+1 gt max_gen then begin
      if debrief gt 0 then print,prompt+'Completed '+sval(max_gen)+' generations. Returning present solution...'
      goto,fH_done
   endif
   igen=igen+1
   if debrief gt 0 then print,prompt+'Computing atomic neutral generation#'+sval(igen)
;________________________________________________________________________________
; Compute Beta_CX from previous generation
;________________________________________________________________________________
   Beta_CX=dblarr(nvr,nvx,nx)
   if H_P_CX then begin
      if debrief gt 1 then print,prompt+'Computing Beta_CX'
      if Simple_CX then begin
;
; Option (B): Compute charge exchange source with assumption that CX source neutrals have
;             ion distribution function
;
         for k=0,nx-1 do Beta_CX(*,*,k)=fi_hat(*,*,k)*total(Vr2pidVr*((alpha_cx(*,*,k)*fHG(*,*,k))#dVx))

;
      endif else begin
;
; Option (A): Compute charge exchange source using fH and vr x sigma x v_v at each velocity mesh point
;
         for k=0,nx-1 do begin
            Work(*)=fHG(*,*,k)
            Beta_CX(*,*,k)=ni(k)*fi_hat(*,*,k)*(SIG_CX#Work)
         endfor
      endelse
;
; Sum charge exchange source over all generations
;
      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   endif
;________________________________________________________________________________
; Compute MH from previous generation
;________________________________________________________________________________
;
   MH_H=dblarr(nvr,nvx,nx)
   MH_P=dblarr(nvr,nvx,nx)
   MH_H2=dblarr(nvr,nvx,nx)
   OmegaM=dblarr(nvr,nvx,nx)
   if H_H_EL or H_P_EL or H_H2_EL then begin
;
; Compute VxHG, THG
;
      for k=0,nx-1 do begin
         VxHG(k)=vth*total(Vr2pidVr*(fHG(*,*,k)#(Vx*dVx)))/NHG(k,igen-1)
         for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-VxHG(k)/vth)^2
         THG(k)=(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fHG(*,*,k))#dVx))/(3*q*NHG(k,igen-1))
      endfor
      if H_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_H'
;
; Compute MH_H  
;
         vx_shift=VxHG
         Tmaxwell=THG
         mol=1
@create_shifted_maxwellian.include
;         Create_Shifted_Maxwellian_Core, vr, vx, vx_shift, Tmaxwell, Maxwell, $
;            vth, Tnorm, Vr2pidVr, dVx, vol, vth_Dvx,   $
;            vx_Dvx, vr_Dvr, vr2vx2, vr2vx2_2D,          $
;            jpa, jpb, jna, jnb, mol, mu, mH, q, shifted_Maxwellian_debug
         for k=0,nx-1 do begin
            MH_H(*,*,k)=Maxwell(*,*,k)*NHG(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_H(k)*MH_H(*,*,k)
         endfor
         MH_H_sum=MH_H_sum+MH_H
      endif
      if H_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_P'
;
; Compute MH_P  
;
         vx_shift=(VxHG+vxi)/2
         Tmaxwell=THG+(2./4.)*(Ti-THG +mu*mH*(vxi-VxHG)^2/(6*q))
         mol=1
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH_P(*,*,k)=Maxwell(*,*,k)*NHG(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_P(k)*MH_P(*,*,k)
         endfor
         MH_P_sum=MH_P_sum+MH_P
      endif
      if H_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_H2'
;
; Compute MH_H2
;
         vx_shift=(VxHG+2*vxH2)/3
         Tmaxwell=THG+(4./9.)*(TH2-THG +2*mu*mH*(vxH2-VxHG)^2/(6*q))
         mol=1
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH_H2(*,*,k)=Maxwell(*,*,k)*NHG(k,igen-1)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_H2(k)*MH_H2(*,*,k)
         endfor
         MH_H2_sum=MH_H2_sum+MH_H2
      endif
   endif
;________________________________________________________________________________
; Compute next generation atomic distribution
;________________________________________________________________________________
;
   fHG(*)=0.0
   for k=0,nx-2 do fHG(*,ip,k+1)=Ak(*,ip,k)*fHG(*,ip,k) $
                               + Bk(*,ip,k)*(Beta_CX(*,ip,k+1)+OmegaM(*,ip,k+1)+Beta_CX(*,ip,k)+OmegaM(*,ip,k))

   for k=nx-1,1,-1 do fHG(*,in,k-1)=Ck(*,in,k)*fHG(*,in,k) $
                               + Dk(*,in,k)*(Beta_CX(*,in,k-1)+OmegaM(*,in,k-1)+Beta_CX(*,in,k)+OmegaM(*,in,k))

   print
   for k=0,nx-1 do nHG(k,igen)=total(Vr2pidVr*(fHG(*,*,k)#dVx))


   if plot gt 1 then begin
      fH1d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH1d(*,k)=Vr2pidVr#fHG(*,*,k)
      endfor
      plot,vx,fH1d(*,0),/nodata,yrange=[0,max(fH1d)],title=sval(igen)+' Generation '+_H
      for i=0,nx-1 do oplot,vx,(fH1d(*,i) > 0.9),color=(i mod 8)+2
      if debug gt 0 then press_return
   endif
;
; Add result to total neutral distribution function
;
   fH=fH+fHG
   nH=nH+nHG(*,igen)
;________________________________________________________________________________
; Compute 'generation error': Delta_nHG=max(NHG(*,igen)/max(nH))
; and decide if another generation should be computed
;________________________________________________________________________________
   Delta_nHG=max(NHG(*,igen)/max(nH))
   if fH_Iterate then begin
      print, 'fH iteration breakpoint'

;
; If fH 'seed' is being iterated, then do another generation until the 'generation error'
; is less than 0.003 times the 'seed error' or is less than TRUNCATE
;
      if (Delta_nHG lt 0.003*Delta_nHs) or (Delta_nHG lt truncate) then begin
         goto, fH_done
      endif
      if (Delta_nHG lt 0.003*Delta_nHs) or (Delta_nHG lt truncate) then goto,fH_done
   endif else begin
;
; If fH 'seed' is NOT being iterated, then do another generation unitl the 'generation error'
; is less than parameter TRUNCATE
;
      if Delta_nHG lt truncate then goto,fH_done
   endelse
;________________________________________________________________________________
;
   print, 'GOING TO NEXT GENERATION'
   goto,next_generation

fH_done:
;
   if plot gt 0 then begin
      plot,x,NHG(*,0),/nodata,yrange=[max(NHG)*truncate,max(NHG)],title=_H+' Density by Generation',xtitle='x (m)',$
           ytitle='Density (m!U-3!N)',/ylog
      for i=0,igen do oplot,x,NHG(*,i),color=(i mod 8)+2
      if pause then press_return
   endif

   ;stop
   ;print, NHG(*,0)
;
; Compute H density profile
   for k=0,nx-1 do nH(k)=total(Vr2pidVr*(fH(*,*,k)#dVx))
;
   if fH_Iterate then begin
;________________________________________________________________________________
; Compute 'seed error': Delta_nHs=(|nHs-nH|)/max(nH) 
; If Delta_nHs is greater than 10*truncate then iterate fH
;________________________________________________________________________________
      Delta_nHs=max(abs(nHs-nH))/max(nH)
      if Delta_nHs gt 10*truncate then begin
         goto,fH_iterate
      endif
   endif
;________________________________________________________________________________
; Update Beta_CX_sum using last generation
;________________________________________________________________________________
   Beta_CX=dblarr(nvr,nvx,nx)
   if H_P_CX then begin
      if debrief gt 1 then print,prompt+'Computing Beta_CX'
      if Simple_CX then begin
;
; Option (B): Compute charge exchange source with assumption that CX source neutrals have
;             ion distribution function
;
         for k=0,nx-1 do Beta_CX(*,*,k)=fi_hat(*,*,k)*total(Vr2pidVr*((alpha_cx(*,*,k)*fHG(*,*,k))#dVx))
;        
      endif else begin
;
; Option (A): Compute charge exchange source using fH and vr x sigma x v_v at each velocity mesh point
;
         for k=0,nx-1 do begin
            Work(*)=fHG(*,*,k)
            Beta_CX(*,*,k)=ni(k)*fi_hat(*,*,k)*(SIG_CX#Work)
         endfor
      endelse
;
; Sum charge exchange source over all generations
;
      Beta_CX_Sum=Beta_CX_Sum+Beta_CX
   endif
;________________________________________________________________________________
; Update MH_*_sum using last generation
;________________________________________________________________________________
;
   MH_H=dblarr(nvr,nvx,nx)
   MH_P=dblarr(nvr,nvx,nx)
   MH_H2=dblarr(nvr,nvx,nx)
   OmegaM=dblarr(nvr,nvx,nx)
   if H_H_EL or H_P_EL or H_H2_EL then begin
;
; Compute VxHG, THG
;
      for k=0,nx-1 do begin
         VxHG(k)=vth*total(Vr2pidVr*(fHG(*,*,k)#(Vx*dVx)))/NHG(k,igen)
         for i=0,nvr-1 do vr2vx2_ran2(i,*)=vr(i)^2+(vx-VxHG(k)/vth)^2
         THG(k)=(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran2*fHG(*,*,k))#dVx))/(3*q*NHG(k,igen))
      endfor
      if H_H_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_H'
;
; Compute MH_H  
;
         vx_shift=VxHG
         Tmaxwell=THG
         mol=1
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH_H(*,*,k)=Maxwell(*,*,k)*NHG(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_H(k)*MH_H(*,*,k)
         endfor
         MH_H_sum=MH_H_sum+MH_H
      endif
      if H_P_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_P'
;
; Compute MH_P  
;
         vx_shift=(VxHG+vxi)/2
         Tmaxwell=THG+(2./4.)*(Ti-THG +mu*mH*(vxi-VxHG)^2/(6*q))
         mol=1
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH_P(*,*,k)=Maxwell(*,*,k)*NHG(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_P(k)*MH_P(*,*,k)
         endfor
         MH_P_sum=MH_P_sum+MH_P
      endif
      if H_H2_EL then begin
         if debrief gt 1 then print,prompt+'Computing MH_H2'
;
; Compute MH_H2
;
         vx_shift=(VxHG+2*vxH2)/3
         Tmaxwell=THG+(4./9.)*(TH2-THG +2*mu*mH*(vxH2-VxHG)^2/(6*q))
         mol=1
@create_shifted_maxwellian.include
         for k=0,nx-1 do begin
            MH_H2(*,*,k)=Maxwell(*,*,k)*NHG(k,igen)
            OmegaM(*,*,k)=OmegaM(*,*,k)+Omega_H_H2(k)*MH_H2(*,*,k)
         endfor
         MH_H2_sum=MH_H2_sum+MH_H2
      endif
   endif
;________________________________________________________________________________
; Compute remaining moments
;________________________________________________________________________________
;
; GammaxH - particle flux in x direction
   for k=0,nx-1 do GammaxH(k)=vth*total(Vr2pidVr*(fH(*,*,k)#(Vx*dVx)))
;
; VxH - x velocity
   VxH=GammaxH/nH
   _VxH=VxH/vth
;
; magnitude of random velocity at each mesh point
   vr2vx2_ran=dblarr(nvr,nvx,nx)
   for i=0,nvr-1 do for k=0,nx-1 do vr2vx2_ran(i,*,k)=vr(i)^2+(vx-_VxH(k))^2
;
; pH - pressure 
   for k=0,nx-1 do pH(k)=(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH(*,*,k))#dVx))/(3*q)
;
; TH - temperature
   TH=pH/nH
;
; piH_xx
   for k=0,nx-1 do piH_xx(k)=(mu*mH)*vth2*total(Vr2pidVr*(fH(*,*,k)#(dVx*(vx-_VxH(k))^2)))/q - pH(k)
;
; piH_yy
   for k=0,nx-1 do piH_yy(k)=(mu*mH)*vth2*0.5*total((Vr2pidVr*Vr^2)*(fH(*,*,k)#dVx))/q - pH(k)
;
; piH_zz
   piH_zz=piH_yy
;
; qxH
   for k=0,nx-1 do qxH(k)=0.5*(mu*mH)*vth3*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*fH(*,*,k))#(dVx*(vx-_VxH(k)))))
;________________________________________________________________________________
; C = RHS of Boltzman equation for total fH
;________________________________________________________________________________
   for k=0,nx-1 do begin
      C=vth*(Sn(*,*,k) + Beta_CX_sum(*,*,k) - alpha_c(*,*,k)*fH(*,*,k) + $
             Omega_H_P(k)*MH_P_sum(*,*,k)+Omega_H_H2(k)*MH_H2_sum(*,*,k)+Omega_H_H(k)*MH_H_sum(*,*,k))
      QH(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2_ran(*,*,k)*C)#dVx))
      RxH(k)=(mu*mH)*vth*total(Vr2pidVr*(C#(dVx*(vx-_VxH(k)))))
      NetHSource(k)=total(Vr2pidVr*(C#dVx))
      Sion(k)=vth*nH(k)*alpha_ion(k)
      SourceH(k)=total(Vr2pidVr*(fSH(*,*,k)#dVx))
      WallH(k)=vth*total(Vr2pidVr*((gamma_wall(*,*,k)*fH(*,*,k))#dVx))
      if Recomb then begin
         SRecomb(k)=vth*ni(k)*Rec(k)
      endif else begin
         SRecomb(k)=0.0
      endelse
      if H_P_CX then begin
         CCX=vth*(Beta_CX_sum(*,*,k) - alpha_cx(*,*,k)*fH(*,*,k))
         RxHCX(k)=(mu*mH)*vth*total(Vr2pidVr*(CCX#(dVx*(vx-_VxH(k)))))
         EHCX(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CCX)#dVx))
      endif
      if H_H2_EL then begin
         CH_H2=vth*Omega_H_H2(k)*(MH_H2_sum(*,*,k)-fH(*,*,k))
         RxH2_H(k)=(mu*mH)*vth*total(Vr2pidVr*(CH_H2#(dVx*(vx-_VxH(k)))))
         EH2_H(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CH_H2)#dVx))
      endif
      if H_P_EL then begin
         CH_P=vth*Omega_H_P(k)*(MH_P_sum(*,*,k)-fH(*,*,k))
         RxP_H(k)=(mu*mH)*vth*total(Vr2pidVr*(CH_P#(dVx*(vx-_VxH(k)))))
         EP_H(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CH_P)#dVx))
      endif
      CW_H=-vth*(gamma_wall(*,*,k)*fH(*,*,k))
      RxW_H(k)=(mu*mH)*vth*total(Vr2pidVr*(CW_H#(dVx*(vx-_VxH(k)))))
      EW_H(k)=0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2vx2(*,*,k)*CW_H)#dVx))
      if H_H_EL then begin
         CH_H=vth*Omega_H_H(k)*(MH_H_sum(*,*,k)-fH(*,*,k))
         for i=0,nvr-1 do vr2_2vx_ran2(i,*)=vr(i)^2-2*(vx-_VxH(k))^2
         Epara_PerpH_H(k)=-0.5*(mu*mH)*vth2*total(Vr2pidVr*((vr2_2vx_ran2*CH_H)#dVx))
      endif
   endfor
;
; qxH_total
;
   qxH_total=(0.5*nH*(mu*mH)*VxH*VxH + 2.5*pH*q)*VxH + q*piH_xx*VxH + qxH
;
; QH_total
;
   QH_total=QH+RxH*VxH + 0.5*(mu*mH)*NetHSource*VxH*VxH
;
; Albedo
;
   AlbedoH=0.0
   gammax_plus=vth*total(Vr2pidVr*(fH(*,ip,0)#(Vx(ip)*dVx(ip))))
   gammax_minus=vth*total(Vr2pidVr*(fH(*,in,0)#(Vx(in)*dVx(in))))
   if abs(gammax_plus) gt 0 then AlbedoH=-gammax_minus/gammax_plus
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
   H_H_error=fltarr(nx,3)
   H_H2_error=fltarr(nx,3)
   H_P_error=fltarr(nx,3)
   max_H_H_error=fltarr(3)
   max_H_H2_error=fltarr(3)
   max_H_P_error=fltarr(3)

;
   if compute_errors then begin
      if debrief gt 1 then print,prompt+'Computing Collision Operator, Mesh, and Moment Normalized Errors'
;
      NetHSource2=SourceH+SRecomb-Sion-WallH
      for k=0,nx-1 do C_error(k)=abs(NetHSource(k)-NetHSource2(k))/max(abs([NetHSource(k),NetHSource2(k)]))

;
; Test conservation of particles for charge exchange operator
;
      if H_P_CX then begin
         for k=0,nx-1 do begin
            CX_A=total(Vr2pidVr*((alpha_cx(*,*,k)*fH(*,*,k))#dVx))
            CX_B=total(Vr2pidVr*((Beta_CX_sum(*,*,k))#dVx))
            CX_error(k)=abs(CX_A-CX_B)/max(abs([CX_A,CX_B]))
         endfor
      endif
;
; Test conservation of particles, x momentum, and total energy of elastic collision operators
;
      for m=0,2 do begin
         for k=0,nx-1 do begin
            if m lt 2 then begin
               TfH=total(Vr2pidVr*(fH(*,*,k)#(dVx*Vx^m)))
            endif else begin
               TfH=total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#dVx))
            endelse
            if H_H_EL then begin
               if m lt 2 then begin
                  TH_H=total(Vr2pidVr*(MH_H_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH_H=total(Vr2pidVr*((vr2vx2(*,*,k)*MH_H_sum(*,*,k))#dVx))
               endelse
               H_H_error(k,m)=abs(TfH-TH_H)/max(abs([TfH,TH_H]))
            endif
            if H_H2_EL then begin
               if m lt 2 then begin
                  TH_H2=total(Vr2pidVr*(MH_H2_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH_H2=total(Vr2pidVr*((vr2vx2(*,*,k)*MH_H2_sum(*,*,k))#dVx))
               endelse
               H_H2_error(k,m)=abs(TfH-TH_H2)/max(abs([TfH,TH_H2]))
            endif
            if H_P_EL then begin
               if m lt 2 then begin
                  TH_P=total(Vr2pidVr*(MH_P_sum(*,*,k)#(dVx*Vx^m)))
               endif else begin
                  TH_P=total(Vr2pidVr*((vr2vx2(*,*,k)*MH_P_sum(*,*,k))#dVx))
               endelse
               H_P_error(k,m)=abs(TfH-TH_P)/max(abs([TfH,TH_P]))
            endif
         endfor
         max_H_H_error(m)=max(H_H_error(*,m))
         max_H_H2_error(m)= max(H_H2_error(*,m))
         max_H_P_error(m)= max(H_P_error(*,m))
      endfor
;
      if CI_test then begin
;
; Compute Momentum transfer rate via full collision integrals for charge exchange and mixed elastic scattering
; Then compute error between this and actual momentum transfer resulting from CX and BKG (elastic) models
;
         if H_P_CX then begin &; P -> H charge exchange momentum transfer via full collision integral
            print,prompt+'Computing P -> H Charge Exchange Momentum Transfer'
            _Sig=dblarr(nvr*nvx*nvr*nvx,ntheta)
            _Sig(*)=v_v*sigma_cx_h0(v_v2*(0.5*mH*vth2/q))
            SIG_VX_CX=dblarr(nvr*nvx,nvr*nvx)
            SIG_VX_CX(*)=Vr2pidVrdVx*vx_vx*(_Sig#dtheta)
            alpha_vx_cx=dblarr(nvr,nvx,nx)
            for k=0,nx-1 do begin
               Work(*)=fi_hat(*,*,k)*ni(k)
               alpha_vx_cx(*,*,k)=SIG_VX_CX#Work
            endfor
            for k=0,nx-1 do RxCI_CX(k)=-(mu*mH)*vth2*total(Vr2pidVr*((Alpha_vx_cx(*,*,k)*fH(*,*,k))#dVx))
            norm=max(abs([RxHCX,RxCI_CX]))
            for k=0,nx-1 do CI_CX_error(k)=abs(RxHCX(k)-RxCI_CX(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in CX collision operator: ',sval(max(CI_CX_Error))
         endif
         if H_P_EL then begin &; P -> H momentum transfer via full collision integral
            for k=0,nx-1 do RxCI_P_H(k)=-(1.0/2.0)*(mu*mH)*vth2*total(Vr2pidVr*((Alpha_H_P(*,*,k)*fH(*,*,k))#dVx))
            norm=max(abs([RxP_H,RxCI_P_H]))
            for k=0,nx-1 do CI_P_H_error(k)=abs(RxP_H(k)-RxCI_P_H(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in P -> H elastic BKG collision operator: ',sval(max(CI_P_H_Error))
         endif
         if H_H2_EL then begin &; H2 -> H momentum transfer via full collision integral
            for k=0,nx-1 do RxCI_H2_H(k)=-(2.0/3.0)*(mu*mH)*vth2*total(Vr2pidVr*((Alpha_H_H2(*,*,k)*fH(*,*,k))#dVx))
            norm=max(abs([RxH2_H,RxCI_H2_H]))
            for k=0,nx-1 do CI_H2_H_error(k)=abs(RxH2_H(k)-RxCI_H2_H(k))/norm
            print,prompt+'Maximum normalized momentum transfer error in H2 -> H elastic BKG collision operator: ',sval(max(CI_H2_H_Error))
         endif
         if H_H_EL then begin &; H -> H perp/parallel energy transfer via full collision integral
            for k=0,nx-1 do begin
               Work(*)=fH(*,*,k)
               Alpha_H_H(*)=SIG_H_H#Work
               Epara_Perp_CI(k)=0.5*(mu*mH)*vth3*total(Vr2pidVr*((Alpha_H_H*fH(*,*,k))#dVx))
            endfor
            norm=max(abs([Epara_PerpH_H,Epara_Perp_CI]))
            for k=0,nx-1 do CI_H_H_error(k)=abs(Epara_PerpH_H(k)-Epara_Perp_CI(k))/norm
            print,prompt+'Maximum normalized perp/parallel energy transfer error in H -> H elastic BKG collision operator: ',sval(max(CI_H_H_Error))
         endif
      endif
;
;  Mesh Point Error based on fH satisfying Boltzmann equation
;
      T1=dblarr(nvr,nvx,nx) & T2=T1 & T3=T1 & T4=T1& T5=T1
      for k=0,nx-2 do begin
         for j=0,nvx-1 do T1(*,j,k)=2*vx(j)*(fH(*,j,k+1)-fH(*,j,k))/(x(k+1)-x(k))
         T2(*,*,k)=Sn(*,*,k+1)+Sn(*,*,k)
         T3(*,*,k)=Beta_CX_sum(*,*,k+1)+Beta_CX_sum(*,*,k)
         T4(*,*,k)=alpha_c(*,*,k+1)*fH(*,*,k+1)+alpha_c(*,*,k)*fH(*,*,k)
         T5(*,*,k)=Omega_H_P(k+1)*MH_P_sum(*,*,k+1)+Omega_H_H2(k+1)*MH_H2_sum(*,*,k+1)+Omega_H_H(k+1)*MH_H_sum(*,*,k+1)+$
                   Omega_H_P(k)*  MH_P_sum(*,*,k)+  Omega_H_H2(k)*  MH_H2_sum(*,*,k)+  Omega_H_H(k)*  MH_H_sum(*,*,k)
         Mesh_Error(*,*,k)=abs(T1(*,*,k)-T2(*,*,k)-T3(*,*,k)+T4(*,*,k)-T5(*,*,k))/$
                              max(abs([T1(*,*,k),T2(*,*,k),T3(*,*,k),T4(*,*,k),T5(*,*,k)]))
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
            moment_error(k,m)=abs(MT1-MT2-MT3+MT4-MT5)/max(abs([MT1,MT2,MT3,MT4,MT5]))
         endfor
         max_moment_error(m)=max(moment_error(*,m))
      endfor
;
; Compute error in qxH_total
;
;    qxH_total2 total neutral heat flux profile (watts m^-2)
;               This is the total heat flux transported by the neutrals
;               computed in a different way from:
;
;               qxH_total2(k)=vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#(Vx*dVx)))*0.5*(mu*mH)
;
;               This should agree with qxH_total if the definitions of nH, pH, piH_xx,
;               TH, VxH, and qxH are coded correctly.
;
      qxH_total2=dblarr(nx)
      for k=0,nx-1 do qxH_total2(k)=0.5*(mu*mH)*vth3*total(Vr2pidVr*((vr2vx2(*,*,k)*fH(*,*,k))#(Vx*dVx)))
      qxH_total_error=abs(qxH_total-qxH_total2)/max(abs([qxH_total,qxH_total2]))
;
; Compute error in QH_total
;
      Q1=dblarr(nx)
      Q2=dblarr(nx)
      QH_total_error=fltarr(nx)
      for k=0,nx-2 do begin
         Q1(k)=(qxH_total(k+1)-qxH_total(k))/(x(k+1)-x(k))
         Q2(k)=0.5*(QH_total(k+1)+QH_total(k))
      endfor
      QH_total_error=abs(Q1-Q2)/max(abs([Q1,Q2]))
;
      if debrief gt 0 then begin
         print,prompt+'Maximum particle convervation error of total collision operator: ',sval(max(C_Error))
         print,prompt+'Maximum H_P_CX  particle convervation error: ',sval(max(CX_Error))
         print,prompt+'Maximum H_H_EL  particle conservation error: '+sval(max_H_H_error(0))
         print,prompt+'Maximum H_H_EL  x-momentum conservation error: '+sval(max_H_H_error(1))
         print,prompt+'Maximum H_H_EL  total energy conservation error: '+sval(max_H_H_error(2))
         print,prompt+'Maximum H_H2_EL particle conservation error: '+sval(max_H_H2_error(0))
         print,prompt+'Maximum H_P_EL  particle conservation error: '+sval(max_H_P_error(0))
         print,prompt+'Average mesh_error =',ave_mesh_error
         print,prompt+'Maximum mesh_error =',max_mesh_error
         for m=0,4 do print,prompt+'Maximum fH vx^'+sval(m)+' moment error: '+sval(max_moment_error(m))
         print,prompt+'Maximum qxH_total error =',max(qxH_total_Error)
         print,prompt+'Maximum QH_total error =',max(QH_total_Error)
         if debug gt 0 then press_return
      endif
   endif

   if plot gt 1 then begin
      fH1d=fltarr(nvx,nx)
      for k=0,nx-1 do begin
         fH1d(*,k)=Vr2pidVr#fH(*,*,k)
      endfor
      plot,vx,fH1d(*,0),yrange=[min(fH1d),max(fH1d)],/nodata,title=_H+' Velocity Distribution Function: fH(Vx)',xtitle='Vx/Vth'
      for i=0,nx-1 do oplot,vx,fH1d(*,i),color=(i mod 6)+2
      if pause then press_return
   endif

   mid1=locate(x,0.7*(max(x)+min(x))/2)
   mid2=locate(x,0.85*(max(x)+min(x))/2)
   mid3=locate(x,(max(x)+min(x))/2)
   mid4=locate(x,1.15*(max(x)+min(x))/2)
   mid5=locate(x,1.3*(max(x)+min(x))/2)
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
      data=[TH,Te,Ti,THP,TH2]
      jp=where(data gt 0)
      yrange=[min(data(jp)),max(data(jp))]
      plot,x,TH,/nodata,/ylog,yrange=yrange,title='Temperature Profiles',xtitle='x (meters)',ytitle='eV'
      oplot,x,TH,color=2
      xyouts,x(mid1),1.2*TH(mid1),_H,color=2
      oplot,x,Te,color=3
      xyouts,1.1*x(mid2),1.2*Ti(mid2),_P,color=5
      oplot,x,Te,color=5
      xyouts,x(mid3),1.2*Te(mid3),'e!U-!N',color=3
      oplot,x,TH2,color=4
      xyouts,x(mid4),1.2*TH2(mid4),_HH,color=4
      oplot,x,THP,color=6
      xyouts,x(mid5),1.2*THP(mid5),_Hp,color=6
      if pause then press_return
   endif
;
   if plot gt 0 then begin
      data=[SourceH,SRecomb,Sion]
      jp=where(data gt 0)
      yrange=[min(data(jp)),max(data(jp))]
      plot,x,SourceH,/nodata,/ylog,yrange=yrange,title=_H+' Source and Sink Profiles',xtitle='x (meters)',ytitle='m!U-3!N s!U-1!N'
      oplot,x,SourceH,color=2
      xyouts,x(mid1),1.2*SourceH(mid1),_HH+' Dissociation',color=2
      oplot,x,SRecomb,color=3
      xyouts,x(mid2),1.2*SRecomb(mid2),_P+' Recombination',color=3
      oplot,x,Sion,color=4
      xyouts,x(mid3),1.2*Sion(mid3),_H+' Ionization',color=4
      oplot,x,WallH,color=6
      xyouts,x(mid4),1.2*WallH(mid4),_H+' Side-Wall Loss',color=6
      if pause then press_return
   endif
;
   if plot gt 0 then begin
      gammaxH_plus=dblarr(nx)
      gammaxH_minus=dblarr(nx)
      for k=0,nx-1 do begin
         gammaxH_plus(k)=vth*total(Vr2pidVr*(fH(*,ip,k)#(Vx(ip)*dVx(ip))))
         gammaxH_minus(k)=vth*total(Vr2pidVr*(fH(*,in,k)#(Vx(in)*dVx(in))))
      endfor  
      data=[gammaxH_plus,gammaxH_minus,gammaxH]
      jp=where(data lt 1.0e32)
      yrange=[min(data(jp)),max(data(jp))]
      plot,x,gammaxH,/nodata,yrange=yrange,title=_H+' Fluxes',xtitle='x (meters)',ytitle='m!U-2!N s!U-1!N'
      oplot,x,gammaxH,color=2
      xyouts,x(mid1),GammaxH(mid1),'!7C!5',color=2
      oplot,x,gammaxH_plus,color=3
      xyouts,x(mid2),GammaxH_plus(mid2),'!7C!5!U+!N',color=3
      oplot,x,GammaxH_minus,color=4
      xyouts,x(mid3),GammaxH_minus(mid3),'!7C!5!U-!N',color=4
      if pause then press_return
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
   fHBC_s=fHBC
   GammaxHBC_s=GammaxHBC
   PipeDia_s=PipeDia
   fH2_s=fH2
   fSH_s=fSH
   nHP_s=nHP
   THP_s=THP
   fH_s=fH
   Simple_CX_s=Simple_CX
   JH_s=JH
   Recomb_s=Recomb
   H_H_EL_s=H_H_EL
   H_P_EL_s=H_P_EL
   H_H2_EL_s=H_H2_EL
   H_P_CX_s=H_P_CX
;
; Set output parameters to single precision
; 
;   fH=float(fH)
;   nH=float(nH)
;   GammaxH=float(GammaxH)
;   VxH=float(VxH)
;   pH=float(pH)
;   TH=float(TH)
;   qxH=float(qxH)
;   qxH_total=float(qxH_total)
;   NetHSource=float(NetHSource)
;   Sion=float(Sion)
;   QH=float(QH)
;   RxH=float(RxH)
;   QH_total=float(QH_total)
;   AlbedoH=float(AlbedoH)
;   piH_xx=float(piH_xx)
;   piH_yy=float(piH_yy)
;   piH_zz=float(piH_zz)
;   RxH2_H=float(RxH2_H)
;   RxP_H=float(RxP_H)
;   RxHCX=float(RxHCX)
;   EH2_H=float(EH2_H)
;   EP_H=float(EP_H)
;   EHCX=float(EHCX)
;   Epara_PerpH_H=float(Epara_PerpH_H)
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
