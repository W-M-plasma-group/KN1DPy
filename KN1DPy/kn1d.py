import numpy as np 
import os
from scipy import interpolate

from .read import sav_read, nc_read
from .edit_keys import edit_keys
from .create_kinetic_h2_mesh import create_kinetic_h2_mesh
from .create_shifted_maxwellian import create_shifted_maxwellian # fixed function name - nh
from .integ_bl import integ_bl
from .make_dvr_dvx import make_dvr_dvx
from .sval import sval
from .interp_fvrvxx import interp_fvrvxx
from .create_kinetic_h_mesh import create_kinetic_h_mesh
from .kinetic_h import kinetic_h 
from .kinetic_h2 import kinetic_h2 
from .interp_scalarx import interp_scalarx 
from .lyman_alpha import lyman_alpha
from .balmer_alpha import balmer_alpha

from .common import constants as CONST
from .common.KN1D import KN1D_Collisions, KN1D_Internal
from .global_vars import global_vars

#   Computes the molecular and atomic neutral profiles for inputted profiles
# of Ti(x), Te(x), n(x), and molecular neutral pressure, GaugeH2, at the boundary using
# IDL routines Kinetic_H and Kinetic_H2. Molecular densities, ionization profiles,
# atomic densities and moments of the atomic distribution function, such as
# T0(x), Qin(x), qx0_total(x),... are returned. 

#   It is assumed that molecular neutrals with temperature equal to the wall temperature
# (~ 1/40 eV) are attacking the plasma at x=x(0).
#
# History: First coding 5/1/2001  -  B. LaBombard
 

def kn1d(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
         error = 0, compute_errors = 0, plot = 0, debug = 0, debrief = 0, pause = 0, \
         Hplot = 0, Hdebug = 0, Hdebrief = 0, Hpause = 0, \
         H2plot = 0, H2debug = 0, H2debrief = 0, H2pause = 0):  # deleted what i added before dont know what I was thinking?? - GG 2/19, corrected typos

    # Input: 
    #	x	- fltarr(nx), cross-field coordinate (meters)
    #      xlimiter - float, cross-field coordinate of limiter edge (meters) (for graphic on plots)
    #	xsep	- float, cross-field coordinate separatrix (meters) (for graphic on plots)
    #	GaugeH2	- float, Molecular pressure (mtorr)
    #	mu	- float, 1=hydrogen, 2=deuterium
    #	Ti	- fltarr(nx), ion temperature profile (eV)
    #	Te	- fltarr(nx), electron temperature profile (eV)
    #	n	- fltarr(nx), density profile (m^-3)
    #	vxi	- fltarr(nx), plasma velocity profile [negative is towards 'wall' (m s^-1)]
    #	LC	- fltarr(nx), connection length (surface to surface) along field lines to nearest limiters (meters)
    #	          Zero values of LC are treated as LC=infinity.
    #      PipeDia	- fltarr(nx), effective pipe diameter (meters)
    #		  This variable allows collisions with the 'side-walls' to be simulated.
    #		  If this variable is undefined, then PipeDia set set to zero. Zero values
    #		  of PipeDia are ignored (i.e., treated as an infinite diameter).
    #
    #   Keyword Input:
    #      truncate	- float, this parameter is also passed to Kinetic_H and Kinetic_H2.
    #                 fH and fH2 are refined by iteration via routines Kinetic_H2 and Kinetic_H
    #		  until the maximum change in molecular neutral density (over its profile) normalized to 
    #		  the maximum value of molecular density is less than this 
    #	    	  value in a subsequent iteration. Default value is 1.0e-3
    #
    #       refine  - if set, then use previously computed atomic and molecular distribution functions
    #		  stored in internal common block (if any) or from FILE (see below) as the initial 
    #                 'seed' value'
    #
    #         file  - string, if not null, then read in 'file'.kn1d_mesh save set and compare contents
    #                 to the present input parameters and computational mesh. If these are the same
    #		  then read results from previous run in 'file'.kn1d_H2 and 'file'.kn1d_H.
    #
    #       Newfile - if set, then do not generate an error and exit if 'file'.KN1D_mesh or 'file'.KN1D_H2
    #                 or 'file'.KN1D_H do not exist or differ from the present input parameters. Instead, write 
    #                 new mesh and output files on exiting.
    #
    #     ReadInput - if set, then reset all input variables to that contained in 'file'.KN1D_input

    g=global_vars() # moved global_vars decleration up
    COLLISIONS = KN1D_Collisions()

    # Collision options inputted via common block KN1D_collisions (default parameter values is true for all collisions):

    # KN1D_Collisions common block - GG 2/15
    H2_H2_EL = COLLISIONS.H2_H2_EL # fixed typo - GG 2/19
    H2_P_EL = COLLISIONS.H2_P_EL
    H2_H_EL = COLLISIONS.H2_H_EL
    H2_HP_CX = COLLISIONS.H2_HP_CX
    H_H_EL = COLLISIONS.H_H_EL
    H_P_EL = COLLISIONS.H_P_EL
    H_P_CX = COLLISIONS.H_P_CX
    Simple_CX = COLLISIONS.Simple_CX

        # H2_H2_EL	- if set, then include H2 -> H2 elastic self collisions
        # H2_P_EL	- if set, then include H2 -> H(+) elastic collisions 
        # H2_H_EL	- if set, then include H2 <-> H elastic collisions 
        # H2_HP_CX	- if set, then include H2 -> H2(+) charge exchange collisions
        # H_H_EL	- if set, then include H -> H elastic self collisions
        # H_P_CX	- if set, then include H -> H(+) charge exchange collisions 
        # H_P_EL	- if set, then include H -> H(+) elastic collisions 
        # Simple_CX	- if set, then use CX source option (B): Neutrals are born
        #              in velocity with a distribution proportional to the local
        #              ion distribution function. Simple_CX=1 is default.

    #   Output:
    #   Molecular info
    #      xH2	- fltarr(nxH2), cross-field coordinate for molecular quantities (meters)
    #      nH2	- fltarr(nxH2), neutral moleular density profile (m^-3)
    #      GammaxH2 - fltarr(nxH2), neutral flux profile (# m^-2 s^-1)
    #      TH2	- fltarr(nxH2), molecular neutral temperature profile (m^-3)
    #    qxH2_total	- fltarr(nxH2), molecular neutral heat flux profile (watts m^-2)
    #      nHP	- fltarr(nxH2), molecular ion density profile (m^-3)
    #      THP	- fltarr(nxH2), molecular ion temperature profile (eV)
    #      SH	- fltarr(nxH2), atomic source profile (m^-3 s^-1)
    #      SP	- fltarr(nxH2), ion source profile (m^-3 s^-1)
    #
    #   Atomic info
    #      xH	- fltarr(nxH), cross-field coordinate for atomic quantities (meters)
    #      nH	- fltarr(nxH), neutral atomic density profile (m^-3)
    #      GammaxH 	- fltarr(nxH), neutral flux profile (# m^-2 s^-1)
    #      TH	- fltarr(nxH), atomic neutral temperature profile (m^-3)
    #    qxH_total	- fltarr(nxH), atomic neutral heat flux profile (watts m^-2)
    #   NetHSource	- fltarr(nxH), net source of atomic neutrals from molecular dissociation and recomb minus ionization (# m^-3) 
    #	Sion	- fltarr(nxH), atomic ionization rate (# m^-3) 
    #	QH_total- fltarr(nxH), net rate of total energy transfer to atomic neutral species (watts m^-3)
    #     SideWallH	- fltarr(nxH), atomic neutral sink rate arising from hitting the 'side walls' (m^-3 s^-1)
    #		  Unlike the molecules in Kinetic_H2, wall collisions result in the destruction of atoms.
    #                 This parameter is used to specify a resulting source of molecular
    #                 neutrals in Kinetic_H2. (molecular source = 2 times SideWallH)
    #	Lyman   - fltarr(nxH), Lyman-alpha emissivity (watts m^-3) using rate coefficients of L.C.Johnson and E. Hinnov
    #	Balmer  - fltarr(nxH), Balmer-alpha emissivity (watts m^-3) using rate coefficients of L.C.Johnson and E. Hinnov


    prompt = 'KN1D => '
    #global xH2,TiM,TeM,nM,PipeDiaM,vxM,vrM,TnormM,xH,TiA,TeA,nA,PipeDiaA,vxA,vrA,TnormA # necessary for setting variables from input_dict (lines 132-133) - nh
    #   Note that these are global to this file only, they are not used / cannot be called in other files

    

    #   KN1D_internal common block - GG 2/15
    
    kn1d_internal = KN1D_Internal()
    fH_s = kn1d_internal.fH_s
    fH2_s = kn1d_internal.fH2_s
    nH2_s = kn1d_internal.nH2_s
    SpH2_s = kn1d_internal.SpH2_s
    nHP_s = kn1d_internal.nHP_s
    THP_s = kn1d_internal.THP_s

    # Set defaults - GG - 2/15
    interp_debug = 0 
    max_gen = 100
    error = 1 # This i think should be moved to the keyword in the call line

        
    # Option: Read input parameters stored in file from previous run
    if ReadInput: 
        input = File
        fp = os.path.exists(input)
        if fp:
            if debrief:
                print(prompt, ' Reading input variables stored in ', input)
                if File[ len(File) - 4 : len(File) ] == '.sav':
                    input_dict = sav_read(File, '//Users/Gwen/Desktop/test.nc')
                elif File[ len(File) - 3 : len(File)]=='.nc': # fixed typo
                    input_dict = nc_read(File)
                else: 
                    if debug:
                        print(prompt, ' Error reading the file') # may be missing some statements 
                        print(prompt, ' finished')
                        # Press_return
                        return 
        else:
            try:
                input_dict=sav_read(File+'.KN1D_input',File+'_asnetcdf.nc') # allows for old-style inputs, might be removed later - nh
            except:
                print(prompt, 'Error reading the file')
                if debug:
                    print(prompt, ' Finished ')
                    return 
        # edit keys for C-Mod Files 
        edit_keys(input_dict) # Not sure if we need this in the final code - not sure the C-Mod files are really intended to be used as inputs - nh
        
        #for i in ['xH2','TiM','TeM','nM','PipeDiaM','vxM','vrM','TnormM','xH','TiA','TeA','nA','PipeDiaA','vxA','vrA','TnormA']:
        #    globals()[i]=input_dict[i.lower()] # Takes entries from input_dict and defines them as variables 
            #   slightly messy, but should solve comment below - nh
        # i think in the end we will have to use the dictionary to manually define the variables
    else:
        # determine optimized vr, vx, grid for kinetc_h2 (molecules, M)
        #nv = 6 NOTE Replaced with constants
        Eneut = np.array([0.003,0.01,0.03,0.1,0.3,1.0,3.0])
        fctr = 0.3
        if GaugeH2 > 15.0:
            fctr = fctr * 15 / GaugeH2

        # replaced function inputs, split output list into variables - nh // fixed keyword inputs - GG
        kh2_mesh = create_kinetic_h2_mesh(CONST.KH2_NV, mu, x, Ti, Te, n, PipeDia, E0 = Eneut, fctr = fctr) 
        xH2,TiM,TeM,nM,PipeDiaM,vxM,vrM,TnormM = kh2_mesh
        
        # determine optimized vr, vx grid for kinetic_h (atoms, A)
        #nv = 10 NOTE Replaced with Constatnts
        fctr = 0.3
        if GaugeH2 > 30.0 :
            fctr = fctr * 30 / GaugeH2

        # finished line since create_kinetic_h_mesh has been programmed - nh // fixed capitalization - GG // fixed keyword inputs - GG
        kh_mesh = create_kinetic_h_mesh(CONST.KH_NV, mu, x, Ti, Te, n, PipeDia, fctr = fctr) 
        xH,TiA,TeA,nA,PipeDiaA,vxA,vrA,TnormA = kh_mesh


    v0_bar=np.sqrt(8.0*CONST.TWALL*CONST.Q/(np.pi*2*mu*CONST.H_MASS))

    #   Set up molecular flux BC from inputted neutral pressure

    ipM=(vxM>0).nonzero()[0]
    inM=(vxM<0).nonzero()[0]
    nvrM=vrM.size
    nvxM=vxM.size
    nxH2=xH2.size
    
    #   Code below uses variables defined in create_kinetic_h_mesh
    
    ipA=(vxA>0).nonzero()[0]
    inA=(vxA<0).nonzero()[0]
    nvrA=vrA.size
    nvxA=vxA.size
    nxH=xH.size
    
    #  Code blocks below use common block variables, should be filled in later
    #  Initialize fH and fH2 (these may be over-written by data from and old run below)
    
    # NOTE Rewrite Deleted, possibly re-add later
    fH = np.zeros((nvrA,nvxA,nxH)).T
    fH2 = np.zeros((nvrM,nvxM,nxH2)).T
    nH2 = np.zeros(nxH2)
    nHP = np.zeros(nxH2)
    THP = np.zeros(nxH2)
        
    #   Convert pressure (mtorr) to molecular density and flux

    fh2BC=np.zeros((nvxM,nvrM)) # fixed mistake in defining the array - GG
    DensM=3.537e19*GaugeH2
    GammaxH2BC=0.25*DensM*v0_bar
    Tmaxwell=np.full(nvxM, CONST.TWALL) # changed list to numpy array 
    vx_shift=np.zeros(nvxM) # fixed size of arrays, its unclear in the original code if this is whats supposed to be done but the for loop in create shifted_maxwellian_include wont owrk otherwise - G
    mol=2
    Maxwell=create_shifted_maxwellian(vrM,vxM,Tmaxwell,vx_shift,mu,mol,TnormM)
    fh2BC[ipM]=Maxwell[0,ipM] # fixed indexing - GG

    # Compute NuLoss:
        # NuLoss = Cs/LC
    Cs_LC=np.zeros(LC.size)
    for ii in range(LC.size):
        if LC[ii]>0:
            Cs_LC[ii]=np.sqrt(CONST.Q*(Ti[ii]+Te[ii])/(mu*CONST.H_MASS))/LC[ii] # fixed notation of indexing arrays - GG
    interpfunc = interpolate.interp1d(x,Cs_LC) # fixed the way interpolation was called - GG
    NuLoss = interpfunc(xH2)
    
    #  Compute first guess SpH2

    #   If plasma recycling accounts for molecular source, then SpH2 = 1/2 n Cs/LC (1/2 accounts for H2 versus H)
    #   But, allow for SpH2 to be proportional to this function:
    #      SpH2 = beta n Cs/LC 
    #   with beta being an adjustable parameter, set by achieving a net H flux of zero onto the wall.
    #   For first guess of beta, set the total molecular source according to the formula
    
    # (See notes "Procedure to adjust the normalization of the molecular source at the 
    #   limiters (SpH2) to attain a net zero atom/molecule flux from wall")
    
    #	Integral{SpH2}dx =  (2/3) GammaxH2BC = beta Integral{n Cs/LC}dx
    
    nCs_LC=n*Cs_LC

    interpfunc = interpolate.interp1d(x,nCs_LC,fill_value="extrapolate") # fixed the way interpolation was called - GG
    SpH2_hat=interpfunc(xH2)

    SpH2_hat=SpH2_hat/integ_bl(xH2,SpH2_hat,value_only=1) # not sure if this line is correct
    beta=2/3*GammaxH2BC
    if refine: # readded this section now that we have the internal common block called - GG 2/15
        if SpH2_s!=None:
            SpH2=SpH2_s # from kn1d_internal common block
        else:
            SpH2=beta*SpH2_hat
    else:
        SpH2=beta*SpH2_hat
    SH2=SpH2

    #   Interpolate for vxiM and vxiA

    interpfunc = interpolate.interp1d(x,vxi,fill_value="extrapolate")
    vxiM = interpfunc(xH2)

    #interpfunc = interpolate.interp1d(x,vxi,fill_value="extrapolate") NOTE Is this Needed?
    vxiA = interpfunc(xH)

    iter=0
    EH_hist=np.array([0.0])
    SI_hist=np.array([0.0])
    oldrun=0

    #   Option: Read results from previous run

        #   Will do later after discussing save/restore

    #   Starting back at line 378 from IDL code
    #   Test for v0_bar consistency in the numerics by computing it from a half maxwellian at the wall temperature

    vthM=np.sqrt(2*CONST.Q*TnormM/(mu*CONST.H_MASS))
    Vr2pidVrM,VrVr4pidVrM,dVxM=make_dvr_dvx(vrM,vxM)[0:3]
    vthA=np.sqrt(2*CONST.Q*TnormA/(mu*CONST.H_MASS))
    Vr2pidVrA,VrVr4pidVrA,dVxA=make_dvr_dvx(vrA,vxA)[0:3]

    nbarHMax=np.sum(Vr2pidVrM*np.matmul(dVxM,fh2BC))
    vbarM=2*vthM*np.sum(Vr2pidVrM*np.matmul(vxM*dVxM,fh2BC))/nbarHMax
    vbarM_error=abs(vbarM-v0_bar)/max(vbarM,v0_bar)

    nvrM=vrM.size
    nvxM=vxM.size
    vr2vx2_ran2=np.zeros((nvrM,nvxM)).T # fixed indexing - GG

    mwell=Maxwell[0,:,:] #  variable named 'Max' in original code; changed here to avoid sharing name with built in function

    nbarMax=np.sum(Vr2pidVrM*np.matmul(dVxM,mwell))
    UxMax=vthM*np.sum(Vr2pidVrM*np.matmul(vxM*dVxM,mwell))/nbarMax
    for i in range(nvrM):
        vr2vx2_ran2[:,i]=vrM[i]**2+(vxM-UxMax/vthM)**2
    TMax=2*mu*CONST.H_MASS*vthM**2*np.sum(Vr2pidVrM*np.matmul(dVxM,vr2vx2_ran2*mwell))/(3*CONST.Q*nbarMax)

    UxHMax=vthM*np.sum(Vr2pidVrM*np.matmul(vxM*dVxM,fh2BC))/nbarHMax
    for i in range(nvrM):
        vr2vx2_ran2[:,i]=vrM[i]**2+(vxM-UxHMax/vthM)**2
    THMax=(2*mu*CONST.H_MASS)*vthM**2*np.sum(Vr2pidVrM*np.matmul(dVxM,vr2vx2_ran2*fh2BC))/(3*CONST.Q*nbarHMax)

    if compute_errors and debrief:
        print(prompt+'VbarM_error: '+sval(vbarM_error))
        print(prompt+'TWall Maxwellian: '+sval(TMax))
        print(prompt+'TWall Half Maxwellian: '+sval(THMax))

    #   Option to view inputted profiles

        #   Plotting - will maybe add later

    #   Starting back at line 429 from IDL code

    print("Satisfaction condition: ", truncate)
        
    if oldrun:
        # checks if the previous run satisfies the required conditions 
        if debrief: 
            print(prompt, 'Maximum Normalized change in nH2: ', sval(nDelta_nH2))
        if debrief and pause: 
            # press_return 
            return
        if nDelta_nH2 > truncate: 
            # goto fH_fH2_iterate I think we will have to make fH_fH2_iterate a function 
            # since we wont be reading old runs right now I am going to leave this as is 
            pass
    else:
        #   Entry point for fH_fH2 iteration : iterates through solving fh and fh2 until they satisfy boltzmans equation
        nDelta_nH2 = truncate + 1
        while nDelta_nH2 > truncate: # Used goto statements in IDL; changed here to while loop
            if debrief: 
                print(prompt, 'Maximum Normalized change in nH2: ', sval(nDelta_nH2))
            if debrief and pause: 
                # press_return
                pass

            # iter+=1 I dont think this line is necessary 
            if debrief:
                print(prompt+'fH/fH2 Iteration: '+sval(iter))
            nH2s = nH2

            # interpolate fH data onto H2 mesh: fH -> fHM
            do_warn=5e-3
            fHM=interp_fvrvxx(fH,vrA,vxA,xH,TnormA,vrM,vxM,xH2,TnormM,do_warn=do_warn, debug=interp_debug, g=g) 

            # Compute fH2 using Kinetic_H2
            ni_correct=1
            Compute_H_Source=1
            H2compute_errors=compute_errors and H2debrief # is this accurate, how can it be equal to both? - GG 2/15
            
            
            kh2_results = kinetic_h2(\
                    vxM, vrM, xH2, TnormM, mu, TiM, TeM, nM, vxiM, fh2BC, GammaxH2BC, NuLoss, PipeDiaM, fHM, SH2, fH2, nH2, THP, \
                    truncate=truncate, Simple_CX=Simple_CX, Max_Gen=max_gen, Compute_H_Source=Compute_H_Source,\
                    H2_H2_EL=H2_H2_EL,H2_P_EL=H2_P_EL,H2_H_EL=H2_H_EL,H2_HP_CX=H2_HP_CX, ni_correct=ni_correct,\
                    Compute_Errors=H2compute_errors, plot=H2plot,debug=H2debug,debrief=H2debrief,pause=H2pause, g=g) # fixed inputs - GG 2/26

            fH2, nHP, THP, nH2, GammaxH2, VxH2, pH2, TH2, qxH2, qxH2_total, Sloss, \
                QH2, RxH2, QH2_total, AlbedoH2, WallH2, fSH, SH, SP, SHP, NuE, NuDis, ESH, Eaxis, error = kh2_results
            

            # Interpolate H2 data onto H mesh: fH2 -> fH2A, fSH -> fSHA, nHP -> nHPA, THP -> THPA
            do_warn = 5.0E-3
            fH2A = interp_fvrvxx(fH2,vrM,vxM,xH2,TnormM,vrA,vxA,xH,TnormA, do_warn=do_warn, debug=interp_debug, g=g) 
            fSHA = interp_fvrvxx(fSH,vrM,vxM,xH2,TnormM,vrA,vxA,xH,TnormA, do_warn=do_warn, debug=interp_debug, g=g) 
            nHPA = interp_scalarx(nHP,xH2,xH, do_warn=do_warn, debug=interp_debug) 
            THPA = interp_scalarx(THP,xH2,xH, do_warn=do_warn, debug=interp_debug)     

            # Compute fH using Kinetic_H
            GammaxHBC = 0
            fHBC = np.zeros((nvxA,nvrA)) # original used (nxH, nvxA, nvrA) but kinetic_h requires shape (nvxA,nvrA)
            H2_H2_EL= H2_H_EL # fixed typo - GG 2/19
            ni_correct = 1
            Hcompute_errors = compute_errors and Hdebrief


            kh_results = kinetic_h(
                vxA,vrA,xH,TnormA,mu,TiA,TeA,nA,vxiA,fHBC,GammaxHBC,PipeDiaA,fH2A,fSHA,nHPA,THPA, fH=fH,\
                    truncate=truncate, Simple_CX=Simple_CX, Max_Gen=max_gen, \
                    H_H_EL=H_H_EL, H_P_EL=H2_P_EL, _H_H2_EL= H2_H2_EL, H_P_CX=H_P_CX, ni_correct=ni_correct, \
                    Compute_Errors=Hcompute_errors, plot=Hplot, debug=Hdebug, debrief=Hdebrief, pause=Hpause, g=g) # Not sure where some of the keywords are defined
            
            fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,SideWallH,error = kh_results


            # Interpolate SideWallH data onto H2 mesh: SideWallH -> SideWallHM
            SideWallHM= interp_scalarx(SideWallH, xH, xH2, do_warn=do_warn, debug=interp_debug)
            # Adjust SpH2 to achieve net zero hydrogen atom/molecule flux from wall
            # (See notes "Procedure to adjust the normalization of the molecular source at the 
            # limiters (SpH2) to attain a net zero atom/molecule flux from wall")

            # Compute SI, GammaH2Wall_minus, and GammaHWall_minus
            SI = integ_bl(xH2, SpH2, value_only=True)
            SwallI = integ_bl(xH2,0.5*SideWallHM, value_only = True)
            GammaH2Wall_minus = AlbedoH2*GammaxH2BC
            GammaHWall_minus=-GammaxH[0]

            # Compute Epsilon and alphaplus1RH0Dis
            Epsilon = 2*GammaH2Wall_minus/(SI+SwallI)
            alphaplus1RH0Dis = GammaHWall_minus/( (1-0.5*Epsilon)*(SI+SwallI)+GammaxH2BC)

            # Compute flux error, EH, and dEHdSI
            EH=2*GammaxH2[0]-GammaHWall_minus
            dEHdSI=-Epsilon-alphaplus1RH0Dis*(1-0.5*Epsilon)

            # Option: print normalized flux error
            nEH=np.abs(EH)/np.max(np.abs( np.array([2*GammaxH2[0],GammaHWall_minus] )))
            if debrief and compute_errors:
                print(prompt, 'Normalized Hydrogen Flux Error: ', sval(nEH))
            
            # Compute Adjustment 
            Delta_SI=-EH/dEHdSI
            SI=SI+Delta_SI

            # Rescale SpH2 to have new integral value, SI
            SpH2=SI*SpH2_hat
            EH_hist=np.append(EH_hist,EH)
            SI_hist=np.append(SI_hist,SI)

            # Set total H2 source
            SH2=SpH2+0.5*SideWallHM

            if compute_errors:
                _RxH_H2 = interp_scalarx(g.Kinetic_H2_Output_RxH_H2,xH2, xH, do_warn=do_warn, debug=interp_debug)
                DRx=_RxH_H2+g.Kinetic_H_Output_RxH2_H
                nDRx=np.max(np.abs(DRx))/np.max(np.abs(np.array([_RxH_H2,g.Kinetic_H_Output_RxH2_H])))
                if debrief:
                    print(prompt, 'Normalized H2 <-> H Momentum Transfer Error: ', sval(nDRx))
            Delta_nH2 = np.abs(nH2-nH2s)
            nDelta_nH2=np.max(Delta_nH2/np.max(nH2))

            print("nDelta_nH2: ", nDelta_nH2)
    
    # fH_fH2_done code section  

    # NOTE May be relevant for plotting, but not for output

    # error = 0
    # # Compute total H flux through crossing limiter radius
    # _GammaxH2 = interp_scalarx(GammaxH2, xH2, xH, do_warn=do_warn, debug=interp_debug)
    # Gam=2*_GammaxH2+GammaxH
    # interpfunc = interpolate.interp1d(xH, Gam, fill_value="extrapolate")
    # GammaHLim = interpfunc(xlimiter)

    # # Compute positive and negative particle flux contributions
    # gammaxH_plus = np.zeros(nxH)
    # gammaxH_minus = np.zeros(nxH)
    # i_p = np.argwhere(vxA > 0).T[0]
    # i_n = np.argwhere(vxA < 0).T[0] # changed formatting to avoid confusion with python function in
    # for k in range(0, nxH):
    #     gammaxH_plus[k] = vthA * np.sum(Vr2pidVrA* np.dot(fH[k][i_p][:], vxA[i_p] * dVxA[i_p]))
    #     gammaxH_minus[k] = vthA * np.sum(Vr2pidVrA* np.dot(fH[k][i_n][:], vxA[i_n] * dVxA[i_n]))
    
    # gammaxH2_plus = np.zeros(nxH2)
    # gammaxH2_minus = np.zeros(nxH2)
    # i_p = np.argwhere(vxM > 0).T[0]
    # i_n = np.argwhere(vxM < 0).T[0]
    # for k in range(0, nxH2):
    #     gammaxH2_plus[k] = vthM * np.sum(Vr2pidVrM* np.dot(fH2[k][i_p][:], vxM[i_p] * dVxM[i_p]))
    #     gammaxH2_minus[k] = vthM * np.sum(Vr2pidVrM* np.dot(fH2[k][i_n][:], vxM[i_n] * dVxM[i_n]))

    
    # Compute Lyman and Balmer
    Lyman = lyman_alpha(nA, TeA, nH, no_null = 1, g=g)
    Balmer = balmer_alpha(nA, TeA, nH, no_null = 1, g=g)

    #   Update KN1D_internal common block - GG 2/15
    #   NOTE What is the purpose of this, can it be removed?
    kn1d_internal.fH_s = fH 
    kn1d_internal.fH2_s = fH2 
    kn1d_internal.nH2_s = nH2
    kn1d_internal.SpH2_s = SpH2
    kn1d_internal.nHP_s = nHP
    kn1d_internal.THP_s = THP

    # The rest of the code for KN1D is for saving files and plotting which we can implement at a later date 
    return xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
        xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer
