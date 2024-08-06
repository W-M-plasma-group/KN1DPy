from KN1D import KN1D
import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time
from loguru import logger

# Configurar el logger (opcional, dependiendo de tus necesidades)
logger.add("file.log", format="{time} {level} {message}", level="INFO")
# Registrar diferentes niveles de mensajes
# logger.debug("Este es un mensaje de depuración, generalmente para desarrollo")
# logger.info("Este es un mensaje informativo")
# logger.success("Este es un mensaje de éxito")
# logger.warning("Este es un mensaje de advertencia")
# logger.error("Este es un mensaje de error")
# logger.critical("Este es un mensaje crítico")

name = "Carlo"
logger.info(f"Lectura del archivo por: {name}")
# read sav file
fname = '1090904024_950to1050.sav'
data_file = readsav(fname)


# params = {
#     'x': data_file['x'],
#     'xlimiter': data_file['x_lim'],
#     'xsep': data_file['x_sep'],
#     'GaugeH2': data_file['p_wall'],
#     'mu': data_file['mu'],
#     'Ti': data_file['t_i'],
#     'Te': data_file['t_e'],
#     'n': data_file['n_e'],
#     'vxi': data_file['vx'],
#     'LC': data_file['lc'],
#     'PipeDia': data_file['d_pipe'],
#     'nv_h2': 20,
#     'nv_h': 20
# }

# xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
# xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, \
# SideWallH, Lyman, Balmer = KN1D(**params)



# logger.info('xH2',len(xH2))
# logger.info('nH2',len(nH2))
# logger.info('data_file[xH2]',len(data_file['xH2']))
# logger.info('data_file[nH2]',len(data_file['nH2']))

# plt.plot(xH2,nH2)
# plt.plot(data_file['xH2'],data_file['nH2'])
# plt.title(r'n$_{H_2}$ Comparison: Shot '+fname[:10])
# plt.yscale('log')
# plt.legend(['Python','IDL'])
# plt.xlabel('x (m)')
# plt.ylabel(r'Density (m$^{-3}$)')
# #plt.savefig(fileloc + fname+'_nH2_testplotb')
# plt.show()
# plt.clf()
# logger.info('xH',len(xH))
# logger.info('nH',len(nH))
# logger.info('data_file[xH]',len(data_file['xH']))
# logger.info('data_file[nH]',len(data_file['nH']))
# plt.plot(xH,nH)
# plt.plot(data_file['xH'],data_file['nH'])
# plt.yscale('log')
# plt.title(r'n$_{H}$ Comparison: Shot '+fname[:10])
# plt.legend(['Python','IDL'])
# plt.xlabel('x (m)')
# plt.ylabel(r'Density (m$^{-3}$)')
# #plt.savefig(fileloc + fname+'_nH_testplotb')
# plt.show()

import time
from scipy.io import readsav

filename = '1090904024_950to1050.sav'
data_file = readsav(filename)

# Valores a iterar para nv_h2 y nv_h
nv_values = [10, 15, 20, 25, 30, 35, 40, 45, 50]

# Diccionario de parámetros base
params = {
    'x': data_file['x'],
    'xlimiter': data_file['x_lim'],
    'xsep': data_file['x_sep'],
    'GaugeH2': data_file['p_wall'],
    'mu': data_file['mu'],
    'Ti': data_file['t_i'],
    'Te': data_file['t_e'],
    'n': data_file['n_e'],
    'vxi': data_file['vx'],
    'LC': data_file['lc'],
    'PipeDia': data_file['d_pipe'],
}

# Función para guardar resultados en un archivo de texto
def save_results(filename, results):
    with open(filename, 'w') as file:
        for name, value in results.items():
            file.write(f"{name}: {value}\n")

# Iterar sobre los valores de nv_h2 y nv_h
for nv_h2 in nv_values:
    for nv_h in nv_values:
        # Actualizar el diccionario de parámetros
        params.update({'nv_h2': nv_h2, 'nv_h': nv_h})

        # Medir el tiempo de inicio
        start_time = time.time()

        # Llamar a la función KN1D
        xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
        xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D(**params)

        # Medir el tiempo de fin
        end_time = time.time()
        
        # Calcular el tiempo de ejecución
        execution_time = end_time - start_time

        # Guardar los resultados en un archivo de texto
        results = {
            'xH2': xH2, 'nH2': nH2, 'GammaxH2': GammaxH2, 'TH2': TH2, 'qxH2_total': qxH2_total,
            'nHP': nHP, 'THP': THP, 'SH': SH, 'SP': SP,
            'xH': xH, 'nH': nH, 'GammaxH': GammaxH, 'TH': TH, 'qxH_total': qxH_total,
            'NetHSource': NetHSource, 'Sion': Sion, 'QH_total': QH_total, 'SideWallH': SideWallH,
            'Lyman': Lyman, 'Balmer': Balmer,
            'tiempo_ejecución': execution_time
        }

        # Nombre del archivo de salida
        output_filename = f"txt_file_nvh_{nv_h}_nvh2_{nv_h2}.txt"

        # Guardar los resultados en el archivo
        save_results(output_filename, results)
