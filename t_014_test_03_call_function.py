# main.py
from t_014_test_01_class import shared_data
from t_014_test_02_function import kinetic_h

# Simulación de un bucle
for i in range(5):
    vx =    [i, i + 1, i + 2]
    vr =    [i, i + 1, i + 2]
    x =     [i, i + 1, i + 2]
    Tnorm = [i, i + 1, i + 2]
    mu =    [i, i + 1, i + 2]

for i in range(5):
    print(66*'-')
    print('Antes')
    print(f'shared_data.vx_s: {shared_data.vx_s}')
    print(f'vx: {vx}')
    result = kinetic_h(vx, vr, x, Tnorm, mu)
    print(f"Iteración {i}: New_Grid = {result}")
    print('Después')
    print(f'shared_data.vx_s: {shared_data.vx_s}')
    print(f'vx: {vx}')
    if i==2:
        print('Reasignando None a shared_data.vx_s')
        shared_data.vx_s = None
        print(f'shared_data.vx_s: {shared_data.vx_s}')
    print(33*'$&')