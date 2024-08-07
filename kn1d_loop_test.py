from KN1D import KN1D
import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time
from loguru import logger

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
import time
from scipy.io import readsav

filename = '1090904024_950to1050.sav'
data_file = readsav(filename)
nv_values = [10, 15, 20, 25, 30, 35, 40, 45, 50]
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
def save_results(filename, results):
    with open(filename, 'w') as file:
        for name, value in results.items():
            file.write(f"{name}: {value}\n")
for nv_h2 in nv_values:
    for nv_h in nv_values:
        params.update({'nv_h2': nv_h2, 'nv_h': nv_h})
        start_time = time.time()

        xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
        xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer = KN1D(**params)

        end_time = time.time()
        execution_time = end_time - start_time
        results = {
            'xH2': xH2, 'nH2': nH2, 'GammaxH2': GammaxH2, 'TH2': TH2, 'qxH2_total': qxH2_total,
            'nHP': nHP, 'THP': THP, 'SH': SH, 'SP': SP,
            'xH': xH, 'nH': nH, 'GammaxH': GammaxH, 'TH': TH, 'qxH_total': qxH_total,
            'NetHSource': NetHSource, 'Sion': Sion, 'QH_total': QH_total, 'SideWallH': SideWallH,
            'Lyman': Lyman, 'Balmer': Balmer,
            'tiempo_ejecución': execution_time
        }
        output_filename = f"txt_file_nvh_{nv_h}_nvh2_{nv_h2}.txt"
        save_results(output_filename, results)
