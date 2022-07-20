import numpy as np 
import os
from read import sav_read, nc_read
from edit_keys import edit_keys
from create_kinetic_h2_mesh import create_kinetic_h2_mesh
from create_shifted_maxwellian import create_shifted_maxwellain
from scipy import interpolate

#   Computes the molecular and atomic neutral profiles for inputted profiles
# of Ti(x), Te(x), n(x), and molecular neutral pressure, GaugeH2, at the boundary using
# IDL routines Kinetic_H and Kinetic_H2. Molecular densities, ionization profiles,
# atomic densities and moments of the atomic distribution function, such as
# T0(x), Qin(x), qx0_total(x),... are returned. 

#   It is assumed that molecular neutrals with temperature equal to the wall temperature
# (~ 1/40 eV) are attacking the plasma at x=x(0).
#
# History: First coding 5/1/2001  -  B. LaBombard
 

def KN1D(x, xlimiter, xsep, GaugeH2, mu, Ti, Te, n, vxi, LC, PipeDia, \
         truncate = 1.0e-3, refine = 0, File = '', NewFile = 0, ReadInput = 0, \
         error = 0, compute_errors = 0, plot = 0, debug = 0, debreif = 0, pause = 0, \
         Hplot = 0, Hdebug = 0, Hdebreif = 0, Hpause = 0, \
         H2plot = 0, H2debug = 0, H2debreif = 0, H2pause = 0):

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
        # resets variables to be values from input file 
        # Option: Read input parameters stored in file from previous run
        if ReadInput: 
            input = File
            fp = os.path.exists(input)
            if fp:
                if debreif:
                    print(prompt, ' Reading input variables stored in ', input)
                    if File[ len(File) - 4 : len(File) ] == '.sav':
                        input_dict = sav_read(File, '//Users/Gwen/Desktop/test.nc')
                    elif File[ len(File) - 3 : len(File)]:
                        input_dict = nc_read(File)
                    else: 
                        if debug:
                            print(prompt, ' Error reading the file')
                            return 
            else:
                    print(prompt, 'Error reading the file')
                    if debug:
                        print(prompt, ' Finished ')
                        return 
            # edit keys for C-Mod Files 
            edit_keys(input_dict) 
            # i think in the end we will have to use the dictionary to manually define the variables
        else:
            # determine optimized vr, vx, grid for kinetc_h2 (molecules, M)
            nv = 6
            Eneut = np.array([0.003,0.01,0.03,0.1,0.3,1.0,3.0])
            fctr = 0.3
            if GaugeH2 > 15.0:
                fctr = fctr * 15 / GaugeH2
            xH2,TiM,TeM,nM,PipeDiaM,vxM,vrM,TnormM = create_kinetic_h2_mesh(n, mu, x, Ti, Te, n, PipeDia) # replaced function inputs, split output list into variables - nh
            
            # determine optimized vr, vx grid for kinetic_h (atoms, A)
            nv = 10
            fctr = 0.3
            if GaugeH2 > 30.0 :
                fctr = fctr * 30 / GaugeH2
            # create_kinetic_H
            
        if mu==1:
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
        else:
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

        plus=' + '
        arrow=' -> '
        elastic=' (elastic)'
        _e='e!U-!N'
        _hv='hv'

        #   constants can probably be referenced from elsewhere

        mH=1.6726231e-27
        q=1.602177e-19
        k_boltz=1.380658e-23
        Twall=293.0*k_boltz/q
        v0_bar=np.sqrt(8.0*Twall*q/(np.pi*2*mu*mH))

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
        
        if refine:
            pass
        else:
            pass
         
         #   Convert pressure (mtorr) to molecular density and flux
         #   (Currently untested)

        fh2BC=np.zeros(nvxM,nvrM)
        DensM=3.537e19*GaugeH2
        GammaxH2BC=0.25*DensM*v0_bar
        Tmaxwell=[Twall]
        vx_shift=[0.0]
        mol=2
        Maxwell=create_shifted_maxwellian(vrM,vxM,Tmaxwell,vx_shift,mu,mol,TnormM)
        fh2BC[ipM,:]=Maxwell[0,ipM,:]

        Cs_LC=np.zeros(LC.size)
        for ii in range(LC.size):
            if LC[ii]>0:
                Cs_LC[ii]=np.sqrt(q*(Ti(ii)+Te(ii))/(mu*mH))/LC(ii)
        NuLoss=interpolate.interp1d(x,Cs_LC)(xH2)
        
        #  Compute first guess SpH2
