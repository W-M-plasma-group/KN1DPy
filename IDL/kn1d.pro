;+
; KN1D.pro
;
;    Computes the molecular and atomic neutral profiles for inputted profiles
; of Ti(x), Te(x), n(x), and molecular neutral pressure, GaugeH2, at the boundary using
; IDL routines Kinetic_H and Kinetic_H2. Molecular densities, ionization profiles,
; atomic densities and moments of the atomic distribution function, such as
; T0(x), Qin(x), qx0_total(x),... are returned. 
;
;     It is assumed that molecular neutrals with temperature equal to the wall temperature
; (~ 1/40 eV) are attacking the plasma at x=x(0).
;
; History: First coding 5/1/2001  -  B. LaBombard
;
;________________________________________________________________________________
;
pro KN1D,x,xlimiter,xsep,GaugeH2,mu,Ti,Te,n,vxi,LC,PipeDia,$
       xH2,nH2,GammaxH2,TH2,qxH2_total,nHP,THP,SH,SP,$
       xH,nH,GammaxH,TH,qxH_total,NetHSource,Sion,QH_total,SideWallH,Lyman,Balmer,$
       GammaHLim,$
       truncate=truncate,refine=refine,File=File,NewFile=NewFile,ReadInput=ReadInput,$
       error=error,compute_errors=compute_errors,$
       plot=plot,debug=debug,debrief=debrief,pause=pause,$
       Hplot=Hplot,Hdebug=Hdebug,Hdebrief=Hdebrief,Hpause=Hpause,$
       H2plot=H2plot,H2debug=H2debug,H2debrief=H2debrief,H2pause=H2pause,$
       H2Gridfctr=H2GridfctrH2,HGridfctr=HGridfctr
;________________________________________________________________________________
;
; Input: 
;	x	- fltarr(nx), cross-field coordinate (meters)
;      xlimiter - float, cross-field coordinate of limiter edge (meters) (for graphic on plots)
;	xsep	- float, cross-field coordinate separatrix (meters) (for graphic on plots)
;	GaugeH2	- float, Molecular pressure (mtorr)
;	mu	- float, 1=hydrogen, 2=deuterium
;	Ti	- fltarr(nx), ion temperature profile (eV)
;	Te	- fltarr(nx), electron temperature profile (eV)
;	n	- fltarr(nx), density profile (m^-3)
;	vxi	- fltarr(nx), plasma velocity profile [negative is towards 'wall' (m s^-1)]
;	LC	- fltarr(nx), connection length (surface to surface) along field lines to nearest limiters (meters)
;	          Zero values of LC are treated as LC=infinity.
;      PipeDia	- fltarr(nx), effective pipe diameter (meters)
;		  This variable allows collisions with the 'side-walls' to be simulated.
;		  If this variable is undefined, then PipeDia set set to zero. Zero values
;		  of PipeDia are ignored (i.e., treated as an infinite diameter).
;
;   Keyword Input:
;      truncate	- float, this parameter is also passed to Kinetic_H and Kinetic_H2.
;                 fH and fH2 are refined by iteration via routines Kinetic_H2 and Kinetic_H
;		  until the maximum change in molecular neutral density (over its profile) normalized to 
;		  the maximum value of molecular density is less than this 
;	    	  value in a subsequent iteration. Default value is 1.0e-3
;
;       refine  - if set, then use previously computed atomic and molecular distribution functions
;		  stored in internal common block (if any) or from FILE (see below) as the initial 
;                 'seed' value'
;
;         file  - string, if not null, then read in 'file'.kn1d_mesh save set and compare contents
;                 to the present input parameters and computational mesh. If these are the same
;		  then read results from previous run in 'file'.kn1d_H2 and 'file'.kn1d_H.
;
;       Newfile - if set, then do not generate an error and exit if 'file'.KN1D_mesh or 'file'.KN1D_H2
;                 or 'file'.KN1D_H do not exist or differ from the present input parameters. Instead, write 
;                 new mesh and output files on exiting.
;
;     ReadInput - if set, then reset all input variables to that contained in 'file'.KN1D_input
;
;  Collision options inputted via common block KN1D_collisions (default parameter values is true for all collisions):
;
    common KN1D_collisions,H2_H2_EL,H2_P_EL,H2_H_EL,H2_HP_CX,H_H_EL,H_P_EL,H_P_CX,Simple_CX
;
;	       H2_H2_EL	- if set, then include H2 -> H2 elastic self collisions
;	        H2_P_EL	- if set, then include H2 -> H(+) elastic collisions 
;	        H2_H_EL	- if set, then include H2 <-> H elastic collisions 
;	       H2_HP_CX	- if set, then include H2 -> H2(+) charge exchange collisions
;	 	 H_H_EL	- if set, then include H -> H elastic self collisions
;	 	 H_P_CX	- if set, then include H -> H(+) charge exchange collisions 
;	         H_P_EL	- if set, then include H -> H(+) elastic collisions 
;             Simple_CX	- if set, then use CX source option (B): Neutrals are born
;                         in velocity with a distribution proportional to the local
;                         ion distribution function. Simple_CX=1 is default.
;             H2Gridfctr - This parameter sets the x grid spacing for the molecular mesh.
;                          Grid spacing is internally optimized according to neutral mean-free-paths.
;                          The default value of H2Gridfctr is 0.3, which results in a adequately fine mesh
;                          spacing in most cases. If a smaller mesh spacing is desired, reduce the value of GridfctrH2
;             HGridfctr -  This parameter performs the same function as H2Gridfctr, but for the atomic species.
;                          Default value is 0.3
; Output:
;   Molecular info
;      xH2	- fltarr(nxH2), cross-field coordinate for molecular quantities (meters)
;      nH2	- fltarr(nxH2), neutral moleular density profile (m^-3)
;      GammaxH2 - fltarr(nxH2), neutral flux profile (# m^-2 s^-1)
;      TH2	- fltarr(nxH2), molecular neutral temperature profile (m^-3)
;    qxH2_total	- fltarr(nxH2), molecular neutral heat flux profile (watts m^-2)
;      nHP	- fltarr(nxH2), molecular ion density profile (m^-3)
;      THP	- fltarr(nxH2), molecular ion temperature profile (eV)
;      SH	- fltarr(nxH2), atomic source profile (m^-3 s^-1)
;      SP	- fltarr(nxH2), ion source profile (m^-3 s^-1)
;
;   Atomic info
;      xH	- fltarr(nxH), cross-field coordinate for atomic quantities (meters)
;      nH	- fltarr(nxH), neutral atomic density profile (m^-3)
;      GammaxH 	- fltarr(nxH), neutral flux profile (# m^-2 s^-1)
;      TH	- fltarr(nxH), atomic neutral temperature profile (m^-3)
;    qxH_total	- fltarr(nxH), atomic neutral heat flux profile (watts m^-2)
;   NetHSource	- fltarr(nxH), net source of atomic neutrals from molecular dissociation and recomb minus ionization (# m^-3) 
;	Sion	- fltarr(nxH), atomic ionization rate (# m^-3 s^-1) 
;	QH_total- fltarr(nxH), net rate of total energy transfer to atomic neutral species (watts m^-3)
;     SideWallH	- fltarr(nxH), atomic neutral sink rate arising from hitting the 'side walls' (m^-3 s^-1)
;		  Unlike the molecules in Kinetic_H2, wall collisions result in the destruction of atoms.
;                 This parameter is used to specify a resulting source of molecular
;                 neutrals in Kinetic_H2. (molecular source = 2 times SideWallH)
;	Lyman   - fltarr(nxH), Lyman-alpha emissivity (watts m^-3) using rate coefficients of L.C.Johnson and E. Hinnov
;	Balmer  - fltarr(nxH), Balmer-alpha emissivity (watts m^-3) using rate coefficients of L.C.Johnson and E. Hinnov
;
;   Combined
;     GammaHLim - float, 2 times GammaxH2 plus GammaxH at edge of limiter (# m^-2 s^-1)
;________________________________________________________________________________
;-
   prompt='KN1D => '
   common KN1D_internal,fH_s,fH2_s,nH2_s,SpH2_s,nHP_s,THP_s
   interp_debug=0
   max_gen=100
   error=1
   key_default,truncate,1.0e-3
   key_default,refine,0
   key_default,File,''
   key_default,NewFile,0
   key_default,ReadInput,0
   key_default,plot,0
   key_default,debug,0
   key_default,debrief,0
   key_default,pause,0
   key_default,Hplot,0
   key_default,Hdebug,0
   key_default,Hdebrief,0
   key_default,Hpause,0
   key_default,H2plot,0
   key_default,H2debug,0
   key_default,H2debrief,0
   key_default,H2pause,0
   key_default,compute_errors,0
;
; Default collisions
;
   key_default,H2_H2_EL,1
   key_default,H2_P_EL,1
   key_default,H2_H_EL,1
   key_default,H2_HP_CX,1
   key_default,H_H_EL,1
   key_default,H_P_EL,1
   key_default,H_P_CX,1
   key_default,Simple_CX,1
   key_default,H2Gridfctr,0.3
   key_default,HGridfctr,0.3

;________________________________________________________________________________
; Option: Read input parameters stored in file from previous run
;________________________________________________________________________________
;
   if ReadInput then begin
      input=file+'.KN1D_input'
      fp=File_Present(Input,errmsg=errmsg)
      if fp then begin
         if debrief then printncr,prompt+' Reading input variables stored in ' & printcr,/rev,input
         restore,input
      endif else begin
         printncr,prompt+' Error reading file ' & printcr,/rev,errmsg
         error=1
         goto,return
      endelse
   endif else begin
;
;________________________________________________________________________________
; Determine optimized vr, vx, x grid for Kinetic_H2 (molecules, M)
;________________________________________________________________________________
      nv=6
      Eneut=[0.003,0.01,0.03,0.1,0.3,1.0,3.0]
      fctr=H2Gridfctr
      if GaugeH2 gt 15.0 then fctr=fctr*15/GaugeH2
      Create_Kinetic_H2_Mesh,nv,mu,x,Ti,Te,n,PipeDia,xH2,TiM,TeM,nM,PipeDiaM,vxM,vrM,TnormM,E0=Eneut,fctr=fctr

;
;________________________________________________________________________________
; Determine optimized vr, vx, x grid for Kinetic_H (atoms, A)
;________________________________________________________________________________
      nv=10
      fctr=HGridfctr
      if GaugeH2 gt 30.0 then fctr=fctr*30/GaugeH2
      Create_Kinetic_H_Mesh,nv,mu,x,Ti,Te,n,PipeDia,xH,TiA,TeA,nA,PipeDiaA,vxA,vrA,TnormA,fctr=fctr


   endelse
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
   _e='e!U-!N'
   _hv='hv'
   mH=1.6726231e-27		
   q=1.602177e-19				
   k_boltz=1.380658e-23				&;Bolzmann's constant, J K^-1
   Twall=293.0*k_boltz/q			&;room temperature (eV)
   v0_bar=sqrt(8.0*Twall*q/(!pi*2*mu*mH))	&;directed random velocity of diatomic molecule
;________________________________________________________________________________
; Set up molecular flux BC from inputted neutral pressure
;________________________________________________________________________________
;
   ipM=where(vxM gt 0)
   inM=where(vxM lt 0)
   nvrM=n_elements(vrM)
   nvxM=n_elements(vxM)
   nxH2=n_elements(xH2)

   ipA=where(vxA gt 0)
   inA=where(vxA lt 0)
   nvrA=n_elements(vrA)
   nvxA=n_elements(vxA)
   nxH=n_elements(xH)
;
; Initialize fH and fH2 (these may be over-written by data from and old run below)
;
   if refine then begin
      if type_of(fH_s) eq 0 then fH=dblarr(nvrA,nvxA,nxH) else fH=fH_s
      if type_of(fH2_s) eq 0 then fH2=dblarr(nvrM,nvxM,nxH2) else fH2=fH2_s
      if type_of(nH2_s) eq 0 then nH2=dblarr(nxH2) else nH2=nH2_s
      if type_of(nHP_s) eq 0 then nHP=dblarr(nxH2) else nHP=nHP_s
      if type_of(THP_s) eq 0 then THP=dblarr(nxH2) else THP=THP_s
   endif else begin
      fH=dblarr(nvrA,nvxA,nxH)
      fH2=dblarr(nvrM,nvxM,nxH2)
      nH2=dblarr(nxH2)
      nHP=dblarr(nxH2)
      THP=dblarr(nxH2)
   endelse
;
; Convert pressure (mtorr) to molecular density and flux
;
   fH2BC=fltarr(nvrM,nvxM)
   DensM=3.537e19*GaugeH2
   GammaxH2BC=0.25*DensM*v0_bar
   Tmaxwell=[Twall]
   vx_shift=[0.0]
   mol=2
   create_shifted_Maxwellian,vrM,vxM,Tmaxwell,vx_shift,mu,mol,TnormM,Maxwell
   fH2BC(*,ipM)=Maxwell(*,ipM,0)
;
; Compute NuLoss 
;   NuLoss = Cs/LC
;
   Cs_LC=LC & Cs_LC(*)=0.0
   ii=where(LC gt 0,count)
   if count gt 0 then begin
      Cs_LC(ii)=sqrt(q*(Ti(ii)+Te(ii))/(mu*mH))/LC(ii)
   endif
   NuLoss=interpol(Cs_LC,x,xH2)
;
; Compute first guess SpH2
;
;   If plasma recycling accounts for molecular source, then SpH2 = 1/2 n Cs/LC (1/2 accounts for H2 versus H)
;   But, allow for SpH2 to be proportional to this function:
;      SpH2 = beta n Cs/LC 
;   with beta being an adjustable parameter, set by achieving a net H flux of zero onto the wall.
;   For first guess of beta, set the total molecular source according to 
;   the formula
;
; (See notes "Procedure to adjust the normalization of the molecular source at the 
;   limiters (SpH2) to attain a net zero atom/molecule flux from wall")
;
;	Integral{SpH2}dx =  (2/3) GammaxH2BC = beta Integral{n Cs/LC}dx
;
   nCs_LC=n*Cs_LC
   SpH2_hat=interpol(nCs_LC,x,xH2)
   SpH2_hat=SpH2_hat/integ_bl(/value,xH2,SpH2_hat)
   beta=(2.0/3.0)*GammaxH2BC
   if refine then begin
      if type_of(SpH2_s) ne 0 then SpH2=SpH2_s else SpH2=beta*SpH2_hat
   endif else begin
      SpH2=beta*SpH2_hat
   endelse
   SH2=SpH2
;
; Interpolate for vxiM and vxiA
;
   vxiM=interpol(vxi,x,xH2)
   vxiA=interpol(vxi,x,xH)
   iter=0
   EH_hist=[0.0]
   SI_hist=[0.0]
   
   oldrun=0
;________________________________________________________________________________
; Option: Read results from previous run
;________________________________________________________________________________
   if strlen(file) gt 0 then begin
;
; Check for old data present
;
      mesh=file+'.KN1D_mesh'
      mp=File_Present(mesh,errmsg=errmsg1)
      H2output=file+'.KN1D_H2'
      H2p=File_Present(H2output,errmsg=errmsg2)
      Houtput=file+'.KN1D_H'
      Hp=File_Present(Houtput,errmsg=errmsg3)
      old_data_present=mp and H2p and Hp
      if old_data_present then begin
         if debrief then printncr,prompt+' Reading mesh variables stored in ' & printcr,/rev,mesh
         restore,mesh
;
; Test if mesh data is the same
;

         test=0
         ii=where(x_s ne x,count) & test=test+count
         ii=where(GaugeH2_s ne GaugeH2,count) & test=test+count
         ii=where(mu_s ne mu,count) & test=test+count
         ii=where(Ti_s ne Ti,count) & test=test+count
         ii=where(Te_s ne Te,count) & test=test+count
         ii=where(n_s ne n,count) & test=test+count
         ii=where(vxi_s ne vxi,count) & test=test+count
         ii=where(PipeDia_s ne PipeDia,count) & test=test+count
         ii=where(LC_s ne LC,count) & test=test+count
         ii=where(xH2_s ne xH2,count) & test=test+count
         ii=where(vxM_s ne vxM,count) & test=test+count
         ii=where(vrM_s ne vrM,count) & test=test+count
         ii=where(TnormM_s ne TnormM,count) & test=test+count
         ii=where(xH_s ne xH,count) & test=test+count
         ii=where(vxA_s ne vxA,count) & test=test+count
         ii=where(vrA_s ne vrA,count) & test=test+count
         ii=where(TnormA_s ne TnormA,count) & test=test+count
         if test le 0 then oldrun=1
         if oldrun then begin
;
; Restore old data
;
            if debrief then printncr,prompt+' Reading output variables stored in ' & printcr,/rev,H2output
            restore,H2output
            if debrief then printncr,prompt+' Reading output variables stored in ' & printcr,/rev,Houtput
            restore,Houtput
            fH_s=fH
            fH2_s=fH2
            nH2_s=nH2
            SpH2_s=SpH2
            nHP_s=nHP
            THP_s=THP
         endif else begin
; Or not
            if not NewFile then begin
               if debrief then printncr,prompt+' Mesh variables from previous run are different! Computing new output...'
            endif
         endelse
      endif else begin
         if not NewFile then begin
            printncr,prompt+' Can not read old data: '& printcr,/rev,errmsg1+' '+errmsg2+' '+errmsg3
            error=1
            goto,return
         endif
      endelse
   endif else begin
      if NewFile then begin
         print,prompt+' ERROR: No file name specified'
         error=1
         goto,return
      endif
   endelse
;
; Test for v0_bar consistency in the numerics by computing it from a half maxwellian at the wall temperature
;
   vthM=sqrt(2*q*TnormM/(mu*mH))
   Make_dVr_dVx,vrM,vxM,Vr2pidVrM,VrVr4pidVrM,dVxM
   vthA=sqrt(2*q*TnormA/(mu*mH))
   Make_dVr_dVx,vrA,vxA,Vr2pidVrA,VrVr4pidVrA,dVxA
   nbarHMax=total(Vr2pidVrM*(fH2BC#dVxM))
   vbarM=2*vthM*total(Vr2pidVrM*(fH2BC#(VxM*dVxM)))/nbarHMax
   vbarM_error=abs(vbarM-v0_bar)/max([vbarM,v0_bar])
   nvrm=n_elements(vrm)
   nvxm=n_elements(vxm)
   vr2vx2_ran2=dblarr(nvrm,nvxm)
   Max=Maxwell(*,*,0)
   nbarMax=total(Vr2pidVrM*(Max#dVxM))
   UxMax=vthM*total(Vr2pidVrM*(Max#(VxM*dVxM)))/nbarMax

   for i=0,nvrm-1 do vr2vx2_ran2(i,*)=vrm(i)^2+(vxm-UxMax/vthM)^2
   TMax=(2*mu*mH)*vthM*vthM*total(Vr2pidVrM*((vr2vx2_ran2*Max)#dVxM))/(3*q*nbarMax)

   UxHMax=vthM*total(Vr2pidVrM*(fH2BC#(VxM*dVxM)))/nbarHMax
   for i=0,nvrm-1 do vr2vx2_ran2(i,*)=vrm(i)^2+(vxm-UxHMax/vthM)^2
   THMax=(2*mu*mH)*vthM*vthM*total(Vr2pidVrM*((vr2vx2_ran2*fH2BC)#dVxM))/(3*q*nbarHMax)


   if compute_errors and debrief then print,prompt+'VbarM_error: '+sval(vbarM_error)
   if compute_errors and debrief then print,prompt+'TWall Maxwellian: '+sval(TMax)
   if compute_errors and debrief then print,prompt+'TWall Half Maxwellian: '+sval(THMax)

;________________________________________________________________________________
; Option to view inputted profiles
;________________________________________________________________________________
   if plot then begin
      multiplot,x,n,x,te,x,ti,x,LC,x,1.001*PipeDia,title='Inputted profiles with '+_HH+' Gauge Presure: '+sval(GaugeH2,l=5)+' mtorr',$
             xtitle='x (meters)',ytop=0.9,ybot=0.1,$
             ytitle=['m!U-2!N','eV','eV','m','m'],color=[2,3,4,6,8],$
             varlabel=['Density','Te','Ti','Connection Length','Pipe Diameter'],logy=[1,1,1,0,0]
      xl=data2normal(!x.crange(0))
      xR=data2normal(!x.crange(1))
      yl=0.1
      yh=0.9
      xlim=data2normal(xlimiter)
      xs=data2normal(xsep)
      h=.05
      plots,[xl,xlim,xlim,xl,xl],[yl,yl,yl+h,yl+h,yl],/norm
      xyouts,0.5*(xl+xlim),0.5*(yl+yl+h),/norm,'LIMITER',charsize=0.8,align=0.5
      xyouts,0.5*(xs+xlim),0.5*(yl+yl+h),/norm,'SOL',charsize=0.8,align=0.5
      xyouts,0.5*(xs+xR),0.5*(yl+yl+h),/norm,'CORE',charsize=0.8,align=0.5
      xhash=xl+(xlim-xl)*findgen(11)/10
      for ii=n_elements(xhash)-1,1,-1 do plots,[xhash(ii),xhash(ii)-h],[yl+h,yl],/norm
      plots,[xlim,xlim],[yl,yh],linestyle=5,/norm
      plots,[xs,xs],[yl,yh],linestyle=5,/norm
      if pause then press_return
   endif
;________________________________________________________________________________
; Entry point for fH/fH2 iteration
;________________________________________________________________________________
;


   if oldrun then goto,test_nDelta_nH2

fH_fH2_iterate:
   iter=iter+1
   if debrief then print,prompt+'fH/fH2 Iteration: '+sval(iter)
   nH2s=nH2

   print, 'breakpoint at the start of the next KN1D iteration'
;
; Interpolate fH data onto H2 mesh: fH -> fHM
;
   warn=5.0e-3
   Interp_fVrVxX,fH,vrA,vxA,xH,TnormA,fHM,vrM,vxM,xH2,TnormM,warn=warn,debug=interp_debug, correct = True ; zzzzz should be true but true breaks for the new version of this
;
; Compute fH2 using Kinetic_H2
;
   ni_correct=1
   Compute_H_Source=1
   H2compute_errors=compute_errors and H2debrief

   print, 'just before going into kinetic_h2'

   Kinetic_H2,vxM,vrM,xH2,TnormM,mu,TiM,TeM,nM,vxiM,fH2BC,GammaxH2BC,NuLoss,PipeDiaM,fHM,SH2,$
       fH2,nH2,GammaxH2,VxH2,pH2,TH2,qxH2,qxH2_total,Sloss,QH2,RxH2,QH2_total,AlbedoH2,WallH2,truncate=truncate,$
       nHP,THP,fSH,SH,SP,SHP,NuE,NuDis,$
       Simple_CX=Simple_CX,Max_Gen=Max_Gen,Compute_H_Source=Compute_H_Source,$
       No_Sawada=No_Sawada,H2_H2_EL=H2_H2_EL,H2_P_EL=H2_P_EL,H2_H_EL=H2_H_EL,H2_HP_CX=H2_HP_CX,ni_correct=ni_correct,$
       ESH=ESH,Eaxis=Eaxis,error=H2error,compute_errors=H2compute_errors,$
       plot=H2plot,debug=H2debug,debrief=H2debrief,pause=H2pause


   common Kinetic_H2_Output,piH2_xx,piH2_yy,piH2_zz,RxH2CX,RxH_H2,RxP_H2,RxW_H2,EH2CX,EH_H2,EP_H2,EW_H2,Epara_PerpH2_H2
   common Kinetic_H2_H_moments,nHM,VxHM,THM

;
; Interpolate H2 data onto H mesh: fH2 -> fH2A, fSH -> fSHA, nHP -> nHPA, THP -> THPA
;
   warn=5.0e-3
   Interp_fVrVxX,fH2,vrM,vxM,xH2,TnormM,fH2A,vrA,vxA,xH,TnormA,warn=warn,debug=interp_debug, correct = True
   Interp_fVrVxX,fSH,vrM,vxM,xH2,TnormM,fSHA,vrA,vxA,xH,TnormA,warn=warn,debug=interp_debug, correct = True
   Interp_ScalarX,nHP,xH2,nHPA,xH,warn=warn,debug=interp_debug
   Interp_ScalarX,THP,xH2,THPA,xH,warn=warn,debug=interp_debug
;
; Compute fH using Kinetic_H
;
   GammaxHBC=0
   fHBC=dblarr(nvrA,nvxA,nxH)
   H_H2_EL=H2_H_EL
   ni_correct=1
   Hcompute_errors=compute_errors and Hdebrief

   

   Kinetic_H,vxA,vrA,xH,TnormA,mu,TiA,TeA,nA,vxiA,fHBC,GammaxHBC,PipeDiaA,fH2A,fSHA,nHPA,THPA,$
       fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,QH,RxH,QH_total,AlbedoH,SideWallH,$
       truncate=truncate,Simple_CX=Simple_CX,Max_Gen=Max_Gen,$
       No_Johnson_Hinnov=No_Johnson_Hinnov,No_Recomb=No_Recomb,$
       H_H_EL=H_H_EL,H_P_EL=H_P_EL,H_H2_EL=H_H2_EL,H_P_CX=H_P_CX,ni_correct=ni_correct,$
       error=Herror,compute_errors=Hcompute_errors,$
       plot=Hplot,debug=Hdebug,debrief=Hdebrief,pause=Hpause

   common Kinetic_H_Output,piH_xx,piH_yy,piH_zz,RxHCX,RxH2_H,RxP_H,RxW_H,EHCX,EH2_H,EP_H,EW_H,Epara_PerpH_H,SourceH,SRecomb
   common Kinetic_H_H2_Moments,nH2A,VxH2A,TH2A
;
; Interpolate SideWallH data onto H2 mesh: SideWallH -> SideWallHM
;
   Interp_ScalarX,SideWallH,xH,SideWallHM,xH2,warn=warn,debug=interp_debug
;
; Adjust SpH2 to achieve net zero hydrogen atom/molecule flux from wall
; (See notes "Procedure to adjust the normalization of the molecular source at the 
;   limiters (SpH2) to attain a net zero atom/molecule flux from wall")
; 
; Compute SI, GammaH2Wall_minus, and GammaHWall_minus
;
   SI=integ_bl(/value,xH2,SpH2)
   SwallI=integ_bl(/value,xH2,0.5*SideWallHM)
   GammaH2Wall_minus=AlbedoH2*GammaxH2BC
   GammaHWall_minus=-GammaxH(0)
;
; Compute Epsilon and alphaplus1RH0Dis
;
   Epsilon=2*GammaH2Wall_minus/(SI+SwallI)
   alphaplus1RH0Dis=GammaHWall_minus/( (1-0.5*epsilon)*(SI+SwallI)+GammaxH2BC)
;
; Compute flux error, EH, and dEHdSI
;
   EH=2*GammaxH2(0)-GammaHWall_minus
   dEHdSI=-Epsilon-alphaplus1RH0Dis*(1-0.5*Epsilon)
;
; Option: print normalized flux error
;
   nEH=abs(EH)/max(abs([2*GammaxH2(0),GammaHWall_minus]))
   if debrief and compute_errors then print,prompt+'Normalized Hydrogen Flux Error: '+sval(nEH)
;
; Compute adjustment in SI, DeltaSI
;
   Delta_SI=-EH/dEHdSI
   SI=SI+Delta_SI
;
; Rescale SpH2 to have new integral value, SI
;
   SpH2=SI*SpH2_hat
   EH_hist=[EH_hist,EH]
   SI_hist=[SI_hist,SI]
;
; Set total H2 source
;
   SH2=SpH2+0.5*SideWallHM

   if compute_errors then begin
      Interp_ScalarX,RxH_H2,xH2,_RxH_H2,xH,warn=warn,debug=interp_debug
      DRx=_RxH_H2+RxH2_H
      nDRx=max(abs(DRx))/max(abs([_RxH_H2,RxH2_H]))
      if debrief then print,prompt+'Normalized H2 <-> H Momentum Transfer Error: '+sval(nDRx)
   endif

   Delta_nH2=abs(nH2-nH2s)
   nDelta_nH2=max(Delta_nH2/max(nH2))

Test_nDelta_nH2:
;
; Test for convergence/iterate
;
   if debrief then print,prompt+'Maximum Normalized change in nH2: ',sval(nDelta_nH2)
   if debrief and pause then press_return
   if nDelta_nH2 lt truncate then goto,fH_fH2_done
   goto,fH_fH2_iterate

fH_fH2_done:
   error=0
;
; Compute total H flux through crossing limiter radius
;

   Interp_ScalarX,GammaxH2,xH2,_GammaxH2,xH,warn=warn,debug=interp_debug
   Gam=2*_GammaxH2+GammaxH
   GammaHLim=interpol(Gam,xH,xlimiter)
;
; Compute positive and negative particle flux contributions
;
   gammaxH_plus=dblarr(nxH)
   gammaxH_minus=dblarr(nxH)
   ip=where(vxa gt 0)
   in=where(vxa lt 0)
   for k=0,nxH-1 do begin
      gammaxH_plus(k)=vthA*total(Vr2pidVrA*(fH(*,ip,k)#(VxA(ip)*dVxA(ip))))
      gammaxH_minus(k)=vthA*total(Vr2pidVrA*(fH(*,in,k)#(VxA(in)*dVxA(in))))
   endfor  
   gammaxH2_plus=dblarr(nxH2)
   gammaxH2_minus=dblarr(nxH2)
   ip=where(vxM gt 0)
   in=where(vxM lt 0)
   for k=0,nxH2-1 do begin
      gammaxH2_plus(k)=vthM*total(Vr2pidVrM*(fH2(*,ip,k)#(VxM(ip)*dVxM(ip))))
      gammaxH2_minus(k)=vthM*total(Vr2pidVrM*(fH2(*,in,k)#(VxM(in)*dVxM(in))))
   endfor  
;
; Compute Lyman and Balmer
;
   Lyman=Lyman_Alpha(nA,TeA,nH,/no_null)
   Balmer=Balmer_Alpha(nA,TeA,nH,/no_null)
;
   fH_s=fH
   fH2_s=fH2
   nH2_s=nH2
   SpH2_s=SpH2
   nHP_s=nHP
   THP_s=THP

   x_s=x
   GaugeH2_s=GaugeH2
   mu_s=mu
   Ti_s=Ti
   Te_s=Te
   n_s=n
   vxi_s=vxi
   PipeDia_s=PipeDia
   LC_s=LC
   xH2_s=xH2
   vxM_s=vxM
   vrM_s=vrM
   TnormM_s=TnormM
   xH_s=xH
   vxA_s=vxA
   vrA_s=vrA
   TnormA_s=TnormA
   EH_hist=EH_hist(1:*)
   SI_hist=SI_hist(1:*)

   if (strlen(file) gt 0) and (iter gt 0) then begin
      input=file+'.KN1D_input'
      if debrief then printncr,prompt+' Saving input variables in ' & printcr,/rev,input
      save,file=input,x,xlimiter,xsep,GaugeH2,mu,Ti,Te,n,vxi,LC,PipeDia,truncate,$
                    xH2,TiM,TeM,nM,PipeDiaM,vxM,vrM,TnormM,$
                    xH,TiA,TeA,nA,PipeDiaA,vxA,vrA,TnormA

      mesh=file+'.KN1D_mesh'
      if debrief then printncr,prompt+' Saving copy of input and mesh variables in ' & printcr,/rev,mesh
      save,file=mesh,x_s,GaugeH2_s,mu_s,Ti_s,Te_s,n_s,vxi_s,LC_s,PipeDia_s,$
                    xH2_s,vxM_s,vrM_s,TnormM_s,xH_s,vxA_s,vrA_s,TnormA_s

      H2output=file+'.KN1D_H2'
      if debrief then printncr,prompt+' Saving H2 output variables in ' & printcr,/rev,H2output
      save,file=H2output,xH2,fH2,nH2,GammaxH2,VxH2,pH2,TH2,qxH2,qxH2_total,Sloss,QH2,RxH2,QH2_total,AlbedoH2,$
              nHP,THP,fSH,SH,SP,SHP,NuE,NuDis,Nuloss,SpH2,$
              piH2_xx,piH2_yy,piH2_zz,RxH2CX,RxH_H2,RxP_H2,RxW_H2,EH2CX,EH_H2,EP_H2,EW_H2,Epara_PerpH2_H2,$

 Gam,GammaxH2_plus,GammaxH2_minus

      Houtput=file+'.KN1D_H'
      if debrief then printncr,prompt+' Saving H output variables in ' & printcr,/rev,Houtput
      save,file=Houtput,xH,fH,nH,GammaxH,VxH,pH,TH,qxH,qxH_total,NetHSource,Sion,SideWallH,QH,RxH,QH_total,AlbedoH,$
              GammaHLim,nDelta_nH2,$
              piH_xx,piH_yy,piH_zz,RxHCX,RxH2_H,RxP_H,RxW_H,EHCX,EH2_H,EP_H,EW_H,Epara_PerpH_H,SourceH,SRecomb,EH_hist,SI_hist,$
              GammaxH_plus,GammaxH_minus,Lyman,Balmer
   endif
   
;
plots:
;
;________________________________________________________________________________
; Optional plots
;________________________________________________________________________________
   xloc=[.15,.3,.45,.6,.75,.9]
   mid=lonarr(n_elements(xloc))
   midH=lonarr(n_elements(xloc))
   midH2=lonarr(n_elements(xloc))
;
; Density
;
   if plot gt 0 then begin
      ydata=[n,nH2,nHP,nH]
      jp=where(ydata gt 0)
      yrange=[(min(ydata(jp)) > 1.0e15),max(ydata(jp))]
      plot,x,n,/nodata,/ylog,yrange=yrange,title='KN1D: Density Profiles',xtitle='x (meters)',ytitle='m!U-3!N'
@kn1d.include
      oplot,x,n,color=2
      xyouts,x(mid(4)),1.1*n(mid(4)),_e,color=2
      oplot,xH2,nH2,color=4
      xyouts,xH2(midH2(1)),1.1*nH2(midH2(1)),_HH,color=4
      oplot,xH,nH,color=3
      xyouts,xH(midH(2)),1.1*nH(midH(2)),_H,color=3
      oplot,xH2,nHP,color=6
      xyouts,xH2(midH2(3)),1.1*nHP(midH2(3)),_Hp,color=6
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Temperature
;
   if plot gt 0 then begin
      ydata=[Ti,Te,TH2,THP,TH]
      jp=where(ydata gt 0)
      yrange=[.02,200.]
      plot,x,Ti,/nodata,/ylog,yrange=yrange,ystyle=1,title='KN1D: Temperature Profiles',xtitle='x (meters)',ytitle='eV'
@kn1d.include
      oplot,x,Ti,color=1
      xyouts,x(mid(5)),1.1*Ti(mid(5)),_p,color=1
      oplot,x,Te,color=2
      xyouts,x(mid(4)),1.1*Te(mid(4)),_e,color=2
      oplot,xH,TH,color=3
      xyouts,xH(midH(2)),1.1*TH(midH(2)),_H,color=3
      oplot,xH2,THP,color=6
      xyouts,xH2(midH2(3)),1.1*THP(midH2(3)),_Hp,color=6
      oplot,xH2,TH2,color=4
      xyouts,xH2(midH2(1)),1.1*TH2(midH2(1)),_HH,color=4
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Particle Fluxes
;
   if plot gt 0 then begin
      GammaxP=n*vxi
      ydata=[2*GammaxH2,GammaxH,GammaxP]
      f=1/1.0e21
      yrange=[min(ydata),1.2*max(ydata)]*f
      plot,x,n*f,/nodata,yrange=yrange,title='KN1D: Particle Fluxes',xtitle='x (meters)',ytitle='10 !U21!N m!U-2!N s!U-1!N'
@kn1d.include
      oplot,xH2,2*GammaxH2*f,color=4
      xyouts,xH2(midH2(0)),1.1*2*GammaxH2(midH2(0))*f,'2x'+_HH,color=4
      oplot,x,GammaxP*f,color=2
      xyouts,x(mid(1)),1.1*GammaxP(mid(1))*f,_e,color=2
      oplot,xH,GammaxH*f,color=3
      xyouts,xH(midH(2)),1.1*GammaxH(midH(2))*f,_H,color=3
      oplot,xH,Gam*f,color=1
      xyouts,xH(midH(3)),1.1*Gam(midH(3))*f,'2x'+_HH+'+'+_H,color=1
@kn1d_limiter.include
      xyouts,.6,.9,'Total '+_H+' Flux at limiter edge: '+string(GammaHLim,format='(E9.2)'),/norm,charsize=.8,color=2
      if pause then press_return
   endif
;
; Positive and Negative Particle Flux Components
;
   if plot gt 0 then begin
      if type_of(gammaxH_plus) eq 0 then begin
         gammaxH_plus=dblarr(nxH)
         gammaxH_minus=dblarr(nxH)
         ip=where(vxa gt 0)
         in=where(vxa lt 0)
         for k=0,nxH-1 do begin
            gammaxH_plus(k)=vthA*total(Vr2pidVrA*(fH(*,ip,k)#(VxA(ip)*dVxA(ip))))
            gammaxH_minus(k)=vthA*total(Vr2pidVrA*(fH(*,in,k)#(VxA(in)*dVxA(in))))
         endfor
         gammaxH2_plus=dblarr(nxH2)
         gammaxH2_minus=dblarr(nxH2)
         ip=where(vxM gt 0)
         in=where(vxM lt 0)
         for k=0,nxH2-1 do begin
            gammaxH2_plus(k)=vthM*total(Vr2pidVrM*(fH2(*,ip,k)#(VxM(ip)*dVxM(ip))))
            gammaxH2_minus(k)=vthM*total(Vr2pidVrM*(fH2(*,in,k)#(VxM(in)*dVxM(in))))
         endfor
      endif
      ydata=[gammaxH_plus,gammaxH_minus,gammaxH,2*gammaxH2_plus,2*gammaxH2_minus,2*gammaxH2]
      f=1/1.0e21
      yrange=[min(ydata),1.2*max(ydata)]*f
      plot,x,n*f,/nodata,yrange=yrange,title='KN1D: Particle Flux Components',xtitle='x (meters)',ytitle='10 !U21!N m!U-2!N s!U-1!N'
@kn1d.include
      oplot,xH2,2*GammaxH2_plus*f,color=3
      xyouts,xH2(midH2(0)),2*GammaxH2_plus(midH2(0))*f,'2x'+_HH+'(+)',color=3
      oplot,xH2,2*GammaxH2_minus*f,color=3
      xyouts,xH2(midH2(0)),2*GammaxH2_minus(midH2(0))*f,'2x'+_HH+'(-)',color=3
      oplot,xH2,2*GammaxH2*f,color=2
      xyouts,xH2(midH2(0)),2*GammaxH2(midH2(0))*f,'2x'+_HH,color=2
      oplot,xH,GammaxH_plus*f,color=4
      xyouts,xH(midH(2)),GammaxH_plus(midH(2))*f,_H+'(+)',color=4
      oplot,xH,GammaxH_minus*f,color=4
      xyouts,xH(midH(2)),GammaxH_minus(midH(2))*f,_H+'(-)',color=4
      oplot,xH,GammaxH*f,color=1
      xyouts,xH(midH(2)),GammaxH(midH(2))*f,_H,color=1
@kn1d_limiter.include
      xyouts,.6,.9,'Total '+_H+' Flux at limiter edge: '+string(GammaHLim,format='(E9.2)'),/norm,charsize=.8,color=2
      if pause then press_return
   endif
;
; Sources/Sinks
;
   if plot gt 0 then begin
      ydata=[SH,SP,SHP,SpH2,SideWallH,NuLoss*nHP,NuDis*nHP,Sion,Srecomb]
      jp=where(ydata gt 0)
      yrange=[min(ydata(jp)),max(ydata(jp))]
      plot,x,n,/nodata,/ylog,yrange=yrange,title='Source(+) and Sink(-) Profiles',$
              ytitle='m!U-3!N s!U-1!N',xtitle='x (meters)'
@kn1d.include
      oplot,xH2,SH,color=2
      xyouts,xH2(midH2(4)),2*SH(midH2(4)),'+'+_H+'('+_HH+')',color=2
      oplot,xH2,SHP,color=8
      xyouts,xH2(midH2(1)),1.2*SHP(midH2(1)),'+'+_Hp+'('+_HH+')',color=8
      oplot,xH2,SP,color=13
      xyouts,xH2(midH2(2)),1.2*SP(midH2(2)),'+'+_p+'('+_HH+')',color=13
      oplot,xH2,NuLoss*nHp,color=5
      xyouts,xH2(midH2(4)),1.2*NuLoss(midH2(4))*nHp(midH2(4)),'-'+_Hp+'(LIM)',color=5
      oplot,xH2,NuDis*nHp,color=6
      xyouts,xH2(midH2(0)),1.2*NuDis(midH2(0))*nHp(midH2(0)),'-'+_Hp+'(Dis)',color=6
      oplot,xH2,SpH2,color=1
      xyouts,xH2(midH2(2)),0.8*SpH2(midH2(2)),'+'+_HH+'(LIM)',color=1
      oplot,xH,SideWallH,color=3
      xyouts,xH(midH(0)),0.8*SideWallH(midH(0)),'+'+_HH+'(Side Wall)',color=3
      oplot,xH,Sion,color=4
      xyouts,xH(midH(5)),1.1*Sion(midH(5)),'-'+_H+'(Ion)',color=4
      oplot,xH,Srecomb,color=3
      xyouts,xH(midH(1)),1.1*Srecomb(midH(1)),'+'+_H+'(Rec)',color=3
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Momentum Transfer
;
   if plot gt 0 then begin
      ydata=[RxH2CX,RxH_H2,RxP_H2,RxW_H2,RxHCX,RxH2_H,RxP_H,RxW_H]
      yrange=[min(ydata),max(ydata)]
      plot,x,n,/nodata,yrange=yrange,title='KN1D: x-Momentum Transfer Rates',xtitle='x (meters)',ytitle='nt m!U-3!N'
@kn1d.include
      oplot,xH2,RxH2CX,color=8
      xyouts,xH2(midH2(3)),RxH2CX(midH2(3)),_HP+'->'+_HH,color=8
      oplot,xH2,RxH_H2,color=6
      xyouts,xH2(midH2(1)),1.2*RxH_H2(midH2(1)),_H+'->'+_HH,color=6
      oplot,xH2,RxP_H2,color=5
      xyouts,xH2(midH2(2)),RxP_H2(midH2(2)),_P+'->'+_HH,color=5
      oplot,xH2,-RxW_H2,color=6
      xyouts,xH2(midH2(0)),-1.2*RxW_H2(midH2(0)),_HH+'->Side Wall',color=6

      oplot,xH,RxHCX,color=4
      xyouts,xH(midH(4)),RxHCX(midH(4)),_P+'->'+_H+'(CX)',color=4
      oplot,xH,RxH2_H,color=3
      xyouts,xH(midH(1)),1.2*RxH2_H(midH(1)),_HH+'->'+_H,color=3
      oplot,xH,RxP_H,color=2
      xyouts,xH(midH(5)),RxP_H(midH(5)),_P+'->'+_H+'(EL)',color=2
      oplot,xH,-RxW_H,color=3
      xyouts,xH(midH(0)),-1.2*RxW_H(midH(0)),_H+'->Side Wall',color=3
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Energy Tranfer
;
   if plot gt 0 then begin
      ydata=[EH2CX,EH_H2,EP_H2,EW_H2,EHCX,EH2_H,EP_H,EW_H]
      f=1/1.0e6
      yrange=[min(ydata),max(ydata)]*f
      plot,x,n*f,/nodata,yrange=yrange,title='KN1D: Energy Transfer Rates',xtitle='x (meters)',ytitle='MW m!U-3!N'
@kn1d.include
      oplot,xH2,EH2CX*f,color=8
      xyouts,xH2(midH2(2)),EH2CX(midH2(2))*f,_HP+'->'+_HH,color=8
      oplot,xH2,EH_H2*f,color=6
      xyouts,xH2(midH2(1)),EH_H2(midH2(1))*f,_H+'->'+_HH,color=6
      oplot,xH2,EP_H2*f,color=5
      xyouts,xH2(midH2(3)),EP_H2(midH2(3))*f,_P+'->'+_HH,color=5
      oplot,xH2,-EW_H2*f,color=7
      xyouts,xH2(midH2(0)),-EW_H2(midH2(0))*f,_HH+'->Side Wall',color=7

      oplot,xH,EHCX*f,color=4
      xyouts,xH(midH(5)),EHCX(midH(5))*f,_P+'->'+_H+'(CX)',color=4
      oplot,xH,EH2_H*f,color=3
      xyouts,xH(midH(4)),EH2_H(midH(4))*f,_HH+'->'+_H,color=3
      oplot,xH,EP_H*f,color=2
      xyouts,xH(midH(4)),EP_H(midH(4))*f,_P+'->'+_H+'(EL)',color=2
      oplot,xH,-EW_H*f,color=9
      xyouts,xH(midH(0)),-EW_H(midH(0))*f,_H+'Side Wall',color=9
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Temperature Isotropization
;
   if plot gt 0 then begin
      ydata=[Epara_PerpH2_H2,Epara_PerpH_H]
      yrange=[min(ydata),max(ydata)]
      plot,x,n,/nodata,yrange=yrange,title='KN1D: T!D//!N -> T!D!9x!3!N Isotropization Rates',$
          xtitle='x (meters)',ytitle='watts m!U-3!N'
@kn1d.include
      oplot,xH2,Epara_PerpH2_H2,color=2
      xyouts,xH2(midH2(2)),Epara_PerpH2_H2(midH2(2)),_HH+'<->'+_HH+'(EL)',color=2
      oplot,xH,Epara_PerpH_H,color=4
      xyouts,xH(midH(3)),Epara_PerpH_H(midH(3)),_H+'<->'+_H+'(EL)',color=4
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Heat Fluxes
;
   if plot gt 0 then begin
      ydata=[qxH_total,qxH2_total]
      f=1/1.0e3
      yrange=[min(ydata),max(ydata)]*f
      plot,x,n*f,/nodata,yrange=yrange,title='KN1D: Heat Fluxes',$
          xtitle='x (meters)',ytitle='KW m!U-2!N'
@kn1d.include
      oplot,xH2,qxH2_total*f,color=2
      xyouts,xH2(midH2(2)),qxH2_total(midH2(2))*f,_HH,color=2
      oplot,xH,qxH_total*f,color=4
      xyouts,xH(midH(3)),qxH_total(midH(3))*f,_H,color=4
@kn1d_limiter.include
      if pause then press_return
   endif
;
; Lyman-alpha and Balmer-alpha emissivities
;
   if plot gt 0 then begin
      ydata=[100*Balmer,Lyman]
      f=1/1.0e3
      yrange=[min(ydata),max(ydata)]*f
      plot,x,n*f,/nodata,yrange=yrange,title='KN1D: H!D!7a!N!3 and L!D!7a!N!3 Emissivities',$
          xtitle='x (meters)',ytitle='KW m!U-3!N'
@kn1d.include
      oplot,xH,100*Balmer*f,color=2
      xyouts,xH(midH(4)),100*Balmer(midH(4))*f,'H!D!7a!N!3x100',color=2
      oplot,xH,Lyman*f,color=4
      xyouts,xH(midH(4)),Lyman(midH(4))*f,'L!D!7a!N!3',color=4
@kn1d_limiter.include
      if pause then press_return
   endif

return:
   if debug then begin
      print,prompt+' finished'
      press_return
   endif
   return
   end
