;+
; jamie_test_make_dvr_dvx.pro
;
; Quick check of Make_dVr_dVx for a 3×3 (vr,vx) grid.
;+
PRO jamie_test_make_dvr_dvx

  ; 1 define a simple vr, vx
  vr = [0.0, 1.0, 2.0]
  vx = [-1.0, 0.0, 1.0]

  ; 2 call the procedure
  Make_dVr_dVx, vr, vx,                $
        Vr2pidVr, VrVr4pidVr, dVx,     $
        vrL=vrL, vrR=vrR,              $
        vxL=vxL, vxR=vxR,              $
        Vol=Vol,                       $
        Vth_DeltaVx=Vth_DVx,           $
        Vx_DeltaVx=Vx_DVx,             $
        Vr_DeltaVr=Vr_DVr,             $
        vr2vx2=vr2vx2,                 $
        jpa=jpa, jpb=jpb, jna=jna, jnb=jnb

  ; 3 print everything
  PRINT, 'Vr2pidVr =', Vr2pidVr
  PRINT, 'VrVr4pidVr =', VrVr4pidVr
  PRINT, 'dVx =', dVx
  PRINT, 'vrL =', vrL, ' vrR =', vrR
  PRINT, 'vxL =', vxL, ' vxR =', vxR
  PRINT, 'Vol ='
  PRINT, Vol
  PRINT, 'Vth_DeltaVx (center slice) =', Vth_DVx[1, *]
  PRINT, 'Vx_DeltaVx (center slice) =', Vx_DVx[1, *]
  PRINT, 'Vr_DeltaVr (center slice) =', Vr_DVr[*, 1]
  PRINT, 'vr2vx2 ='
  PRINT, vr2vx2
  PRINT, 'jpa, jpb, jna, jnb =', jpa, jpb, jna, jnb
END
