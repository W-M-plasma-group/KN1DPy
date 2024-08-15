; Create_Shifted_Maxwellian.include
;
; Este archivo INCLUDE es utilizado por Kinetic_H2.pro y Kinetic_H.pro
;
; Entrada:
;   Vx_shift  - arreglo de doble precisión (dblarr(nx)), (m s^-1)
;   Tmaxwell  - arreglo de doble precisión (dblarr(nx)), (eV)
;   Shifted_Maxwellian_Debug - si está activado, se imprime información de depuración
;   mol       - 1=átomo, 2=molécula diatómica
;
; Salida:
;   Maxwell   - arreglo de doble precisión (dblarr(nvr,nvx,nx)) que representa una función de distribución Maxwelliana desplazada
;               con momento numéricamente evaluado vx cercano a Vx_shift y temperatura cercana a Tmaxwell

; Inicialización de arreglos y variables
AN = dblarr(nvr, nvx, 2)             ; Arreglo para AN
BN = dblarr(nvr, nvx, 2)             ; Arreglo para BN
sgn = [1, -1]                       ; Signos para comparación
Maxwell = 0.0                       ; Inicializar Maxwell

; Iterar sobre todos los valores de k
for k = 0, nx-1 do begin
    ; Si Tmaxwell(k) es mayor que 0, calcular Maxwell
    if Tmaxwell(k) gt 0.0 then begin
        ; Calcular Maxwell
        for i = 0, nvr-1 do begin
            arg = -(vr(i)^2 + (vx - Vx_shift(k)/vth)^2) * mol * Tnorm / Tmaxwell(k)
            arg = arg < 0.0                    ; Esta línea indica que no se recogen los elementos de arg menores a 0
            Maxwell(i, *, k) = exp(arg > (-80))
        endfor

        ; Normalizar Maxwell
        Maxwell(*, *, k) = Maxwell(*, *, k) / total(Vr2pidVr * (Maxwell(*, *, k) # dVx))

        ; Si Shifted_Maxwellian_Debug está activado, calcular errores
        if shifted_Maxwellian_debug then begin
            vx_out1 = vth * total(Vr2pidVr * (Maxwell(*, *, k) # (Vx * dVx)))
            for i = 0, nvr-1 do vr2vx2_ran2(i, *) = vr(i)^2 + (vx - vx_out1 / vth)^2
            T_out1 = (mol * mu * mH) * vth2 * total(Vr2pidVr * ((vr2vx2_ran2 * Maxwell(*, *, k)) # dVx)) / (3 * q)
            vth_local = 0.1 * sqrt(2 * Tmaxwell(k) * q / (mol * mu * mH))
            Terror = abs(Tmaxwell(k) - T_out1) / Tmaxwell(k)
            Verror = abs(vx_out1 - Vx_shift(k)) / vth_local
        endif

        ; Calcular momentos deseados
        WxD = Vx_shift(k)
        ED = WxD^2 + 3 * q * Tmaxwell(k) / (mol * mu * mH)

        ; Calcular momentos actuales de Maxwell
        WxMax = vth * total(Vr2pidVr * (Maxwell(*, *, k) # (Vx * dVx)))
        EMax = vth2 * total(Vr2pidVr * ((vr2vx2_2D * Maxwell(*, *, k)) # dVx))

        ; Calcular Nij y agregar ceros alrededor
        Nij = dblarr(nvr+2, nvx+2)
        Nij(1:nvr, 1:nvx) = Maxwell(*, *, k) * vol

        ; Calcular variaciones de Nij
        Nijp1_vx_Dvx = shift(Nij * vx_Dvx, 0, -1)
        Nij_vx_Dvx = Nij * vx_Dvx
        Nijm1_vx_Dvx = shift(Nij * vx_Dvx, 0, 1)
        Nip1j_vr_Dvr = shift(Nij * vr_Dvr, -1, 0)
        Nij_vr_Dvr = Nij * vr_Dvr
        Nim1j_vr_Dvr = shift(Nij * vr_Dvr, 1, 0)

        ; Calcular AN y BN
        _AN = shift(Nij * vth_Dvx, 0, 1) - Nij * vth_Dvx
        AN(*, *, 0) = _AN(1:nvr, 1:nvx)
        _AN = -shift(Nij * vth_Dvx, 0, -1) + Nij * vth_Dvx
        AN(*, *, 1) = _AN(1:nvr, 1:nvx)

        BN(*,   jpa+1:jpb   ,0)= Nijm1_vx_Dvx(  1:nvr,  jpa+2:jpb+1 )-Nij_vx_Dvx(1:nvr  ,jpa+2:jpb+1)
        BN(*,   jpa         ,0)=-Nij_vx_Dvx(    1:nvr,  jpa+1       )
        BN(*,   jnb         ,0)= Nij_vx_Dvx(    1:nvr,  jnb+1       )
        BN(*,   jna:jnb-1   ,0)=-Nijp1_vx_Dvx(  1:nvr,  jna+1:jnb   )+Nij_vx_Dvx(1:nvr  ,jna+1:jnb  )
        BN(*,   *           ,0)= Nim1j_vr_Dvr(  1:nvr,      1:nvx   )-Nij_vr_Dvr(1:nvr  ,    1:nvx  ) + BN(*,*,0)

        BN(*,       jpa+1:jpb   ,1)=-Nijp1_vx_Dvx(1:nvr,jpa+2:jpb+1 )+Nij_vx_Dvx(1:nvr, jpa+2:jpb+1)
        BN(*,       jpa         ,1)=-Nijp1_vx_Dvx(1:nvr,jpa+1       )
        BN(*,       jnb         ,1)= Nijm1_vx_Dvx(1:nvr,jnb+1       )
        BN(*,       jna  :jnb-1 ,1)= Nijm1_vx_Dvx(1:nvr,jna+1:jnb   )-Nij_vx_Dvx(1:nvr, jna+1:jnb  )
        BN(1:nvr-1, *           ,1)=-Nip1j_vr_Dvr(2:nvr,    1:nvx   )+Nij_vr_Dvr(2:nvr,     1:nvx  ) +  BN(1:nvr-1, *,1)
        BN(0,       *           ,1)=-Nip1j_vr_Dvr(1    ,    1:nvx   )+                                  BN(0,       *,1) 

        ; Eliminar ceros agregados en Nij
        Nij = Nij(1:nvr, 1:nvx)

        ; Probar todas las combinaciones de signos para a_Max y b_Max
        TB1 = fltarr(2)
        TB2 = fltarr(2)
        ia = 0
        while ia lt 2 do begin
            TA1 = vth * total(AN(*, *, ia) # Vx)
            TA2 = vth2 * total(vr2vx2_2D * AN(*, *, ia))
            ib = 0
            while ib lt 2 do begin
                if TB1(ib) eq 0 then TB1(ib) = vth * total(BN(*, *, ib) # Vx)
                if TB2(ib) eq 0 then TB2(ib) = vth2 * total(vr2vx2_2D * BN(*, *, ib))
                denom = TA2 * TB1(ib) - TA1 * TB2(ib)
                b_Max = 0.0
                a_Max = 0.0
                if denom ne 0.0 and TA1 ne 0.0 then begin
                    b_Max = (TA2 * (WxD - WxMax) - TA1 * (ED - EMax)) / denom
                    a_Max = (WxD - WxMax - TB1(ib) * b_Max) / TA1
                endif
                if a_Max * sgn(ia) gt 0.0 and b_Max * sgn(ib) gt 0.0 then begin
                    Maxwell(*, *, k) = (Nij + AN(*, *, ia) * a_Max + BN(*, *, ib) * b_Max) / vol
                    ia = 2
                    ib = 2
                endif
                ib = ib + 1
            endwhile
            ia = ia + 1
        endwhile

        ; Normalizar Maxwell después de la corrección
        Maxwell(*, *, k) = Maxwell(*, *, k) / total(Vr2pidVr * (Maxwell(*, *, k) # dVx))

        ; Si Shifted_Maxwellian_Debug está activado, imprimir errores
        if shifted_Maxwellian_debug then begin
            vx_out2 = vth * total(Vr2pidVr * (Maxwell(*, *, k) # (Vx * dVx)))
            for i = 0, nvr-1 do vr2vx2_ran2(i, *) = vr(i)^2 + (vx - vx_out2 / vth)^2
            T_out2 = (mol * mu * mH) * vth2 * total(Vr2pidVr * ((vr2vx2_ran2 * Maxwell(*, *, k)) # dVx)) / (3 * q)
            Terror2 = abs(Tmaxwell(k) - T_out2) / Tmaxwell(k)
            Verror2 = abs(Vx_shift(k) - vx_out2) / vth_local
            print, 'CREATE_SHIFTED_MAXWELLIAN=> Terror: ' + sval(Terror) + ' -> ' + sval(Terror2) + '  Verror: ' + sval(Verror) + ' -> ' + sval(Verror2)
        endif
    endif
endfor
