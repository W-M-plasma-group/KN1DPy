; Runs KN1D with inputs from sav file

    PRINT, 'Clearing IDL'
    .reset_session
    CLEAR
    
    PRINT, 'Compiling KN1D'
    .r kn1d

    sav_name = '1090904029_950to1050_towall.sav'
    PRINT, 'LOADING SAV FILE: ' + sav_name
    savfile = FILEPATH(sav_name, ROOT_DIR='/home/jduns/KN1D_python/IDL_version/', SUBDIR=['savfiles'])
    sObj = OBJ_NEW('IDL_Savefile', savfile)
    sNames = sObj->Names()
    PRINT, 'SAV FILE VARIABLE NAMES:'
    PRINT, '     ', sNames
    PRINT, ''

    PRINT, 'RUNNING KN1D'


;   Restoring variables from savfile
    RESTORE, 'savfiles/' + sav_name

;   Settings

    ;if !d.name ne 'PS' then tek
    newfile=1
    file='test_kn1d'
    plot=0
    pause=0
    Hpause=0
    H2pause=0
    debug=0
    Hdebug=0
    H2debug=0
    refine=0
    key_default,pause,0
    compute_errors=1
    debrief=0
    Hdebrief=1
    H2debrief=1

    common KN1D_collisions,H2_H2_EL,H2_P_EL,H2_H_EL,H2_HP_CX,H_H_EL,H_P_EL,H_P_CX,Simple_CX
    ; H_H_EL=0
    ; H2_H2_EL=0  
    ; Simple_CX = 0

    t0 = SYSTIME(/SECONDS)
    KN1D,X,X_LIM,X_SEP,P_WALL,MU,T_I,T_E,N_E,VX,LC,D_PIPE,$
        xH2,nH2,GammaxH2,TH2,qxH2_total,nHP,THP,SH,SP,$
        xH,nH,GammaxH,TH,qxH_total,NetHSource,Sion,QH_total,SideWallH,Lyman,Balmer,$
        GammaHLim,$
        truncate=truncate,refine=refine,File=File,NewFile=NewFile,ReadInput=ReadInput,$
        error=error,compute_errors=compute_errors,$
        plot=plot,debug=debug,debrief=debrief,pause=pause,$
        Hplot=Hplot,Hdebug=Hdebug,Hdebrief=Hdebrief,Hpause=Hpause,$
        H2plot=H2plot,H2debug=H2debug,H2debrief=H2debrief,H2pause=H2pause
    t1 = SYSTIME(/SECONDS)
    PRINT, 'Elapsed time (s): ', t1 - t0

    ;exit