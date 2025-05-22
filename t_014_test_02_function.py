# kinetic_h.py

from t_014_test_01_class import shared_data  # Importamos la instancia global

def kinetic_h(vx, vr, x, Tnorm, mu):
    # Inicializar New_Grid
    New_Grid = 1
    if shared_data.vx_s is None:
        print('shared_data.vx_s No está definida')
    # Verificar si vx_s está definida (no es None)
    if shared_data.vx_s is not None:
        print('shared_data.vx_s Sí está definida')
        test = 0

        # Comparar vx_s con vx
        test += sum(1 for a, b in zip(shared_data.vx_s, vx) if a != b)
        # Comparar vr_s con vr
        test += sum(1 for a, b in zip(shared_data.vr_s, vr) if a != b)
        # Comparar x_s con x
        test += sum(1 for a, b in zip(shared_data.x_s, x) if a != b)
        # Comparar Tnorm_s con Tnorm
        test += sum(1 for a, b in zip(shared_data.Tnorm_s, Tnorm) if a != b)
        # Comparar mu_s con mu
        test += sum(1 for a, b in zip(shared_data.mu_s, mu) if a != b)

        # Si no hay diferencias, establecer New_Grid = 0
        if test <= 0:
            New_Grid = 0

    # Actualizar las variables compartidas con los valores actuales
    shared_data.vx_s = vx
    shared_data.vr_s = vr
    shared_data.x_s = x
    shared_data.Tnorm_s = Tnorm
    shared_data.mu_s = mu

    return New_Grid
