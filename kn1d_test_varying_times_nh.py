from KN1D import KN1D
import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time
import copy
from loguru import logger
logger.add("file.log", format="{time} {level} {message}", level="INFO")

fname='1090904024_950to1050.sav'
data_file = readsav(fname)
def save_results(filename, results):
    with open(filename, 'w') as file:
        for name, value in results.items():
            file.write(f"{name}: {value}\n")
for times_nh in range(2,13,1):
    logger.info(f"time_nh = {times_nh}")
    x_rd = np.linspace(data_file['x'][0],data_file['x'][-1],times_nh*len(data_file['x']))
    ####
    func_ti = interpolate.interp1d(data_file['x'],data_file['t_i'],fill_value='extrapolate')
    func_te = interpolate.interp1d(data_file['x'],data_file['t_e'],fill_value='extrapolate')
    func_ne = interpolate.interp1d(data_file['x'],data_file['n_e'],fill_value='extrapolate')
    func_d_pipe = interpolate.interp1d(data_file['x'],data_file['d_pipe'],fill_value='extrapolate')
    func_lc = interpolate.interp1d(data_file['x'],data_file['lc'],fill_value='extrapolate')
    func_vx = interpolate.interp1d(data_file['x'],data_file['vx'],fill_value='extrapolate')
    logger.info('Las interpolaciones fueron creadas')
    ####
    ti_rd = copy.copy(func_ti(x_rd))
    te_rd = copy.copy(func_te(x_rd))
    ne_rd = copy.copy(func_ne(x_rd))
    d_pipe_rd = copy.copy(func_d_pipe(x_rd))
    lc_rd = copy.copy(func_lc(x_rd))
    vx_rd = copy.copy(func_vx(x_rd))
    logger.info('Variables refinadas')
    logger.info(f'len(xh2)= {len(x_rd)}')
    ####
    ####
    t0=time.time()
    params = {
        'x':        x_rd,
        'xlimiter': data_file['x_lim'],
        'xsep':     data_file['x_sep'],
        'GaugeH2':  data_file['p_wall'],
        'mu':       data_file['mu'],
        'Ti':       ti_rd,
        'Te':       te_rd,
        'n':        ne_rd,
        'vxi':      vx_rd,
        'LC':       lc_rd,
        'PipeDia':  d_pipe_rd,
        'nv_h':     40,
        'nv_h2':    30,
            }
    xH2, nH2, GammaxH2, TH2, qxH2_total, nHP, THP, SH, SP, \
    xH, nH, GammaxH, TH, qxH_total, NetHSource, Sion, QH_total, SideWallH, Lyman, Balmer=\
                KN1D(**params)
    
    logger.info(f'len(xH2)_KN1D= {len(xH2)}')
    tf=time.time()
    t=tf-t0

    results = {
            'xH2': xH2, 'nH2': nH2, 'GammaxH2': GammaxH2, 'TH2': TH2, 'qxH2_total': qxH2_total,
            'nHP': nHP, 'THP': THP, 'SH': SH, 'SP': SP,
            'xH': xH, 'nH': nH, 'GammaxH': GammaxH, 'TH': TH, 'qxH_total': qxH_total,
            'NetHSource': NetHSource, 'Sion': Sion, 'QH_total': QH_total, 'SideWallH': SideWallH,
            'Lyman': Lyman, 'Balmer': Balmer,
            'tiempo_ejecución': t
        }
    output_filename = f"txt_file_times_nh_{times_nh}.txt"
    save_results(output_filename, results)
    logger.info(f'Completed in {str(t)} seconds')

    print('xH2',len(xH2))
    print('nH2',len(nH2))
    print('data_file[xH2]',len(data_file['xH2']))
    print('data_file[nH2]',len(data_file['nH2']))


    plt.plot(xH2,nH2)
    plt.plot(data_file['xH2'],data_file['nH2'])
    plt.title(r'n$_{H_2}$ Comparison: Shot '+fname[:10])
    plt.yscale('log')
    plt.legend(['Python','IDL'])
    plt.xlabel('x (m)')
    plt.ylabel(r'Density (m$^{-3}$)')
    #plt.savefig(fileloc + fname+'_nH2_testplotb')
    plt.show()
    plt.clf()
    print('xH',len(xH))
    print('nH',len(nH))
    print('data_file[xH]',len(data_file['xH']))
    print('data_file[nH]',len(data_file['nH']))
    plt.plot(xH,nH)
    plt.plot(data_file['xH'],data_file['nH'])
    plt.yscale('log')
    plt.title(r'n$_{H}$ Comparison: Shot '+fname[:10])
    plt.legend(['Python','IDL'])
    plt.xlabel('x (m)')
    plt.ylabel(r'Density (m$^{-3}$)')
    plt.show()