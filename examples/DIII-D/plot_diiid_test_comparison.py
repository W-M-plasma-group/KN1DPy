"""
DIII-D KN1D comparison plots — Python (KN1DPy) vs IDL
======================================================
Run from the KN1DPy root directory:

    python examples/DIII-D/plot_diiid_test_comparison.py

Loads outputs from:
  - examples/DIII-D/python_output/    (Python run)
  - examples/DIII-D/IDL_output/       (IDL run)

and plots the key profiles side-by-side for comparison.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import readsav

# ------------------------------------------------------------------ #
#  Load Python outputs
# ------------------------------------------------------------------ #
py_h  = np.load('examples/DIII-D/python_output/KN1D_H.npz')
py_h2 = np.load('examples/DIII-D/python_output/KN1D_H2.npz')
inp   = np.load('examples/DIII-D/python_output/KN1D_input.npz')

xA  = inp['xH']    # H mesh x coordinates
Ti  = inp['TiA']   # ion temperature on H mesh
Te  = inp['TeA']   # electron temperature on H mesh
n   = inp['nA']    # electron density on H mesh

xH     = py_h['xH'];      nH     = py_h['nH'];      TH     = py_h['TH']
GamxH  = py_h['GammaxH']; Sion   = py_h['Sion'];    Lyman  = py_h['Lyman']
Balmer = py_h['Balmer'];  qxH    = py_h['qxH_total']
SideWallH = py_h['SideWallH'];  SRecomb = py_h['SRecomb']

xH2    = py_h2['xH2'];    nH2    = py_h2['nH2'];    TH2    = py_h2['TH2']
GamxH2 = py_h2['GammaxH2'];      nHP = py_h2['nHP']; THP = py_h2['THP']
SH     = py_h2['SH'];     SP     = py_h2['SP'];      qxH2   = py_h2['qxH2_total']
NuLoss = py_h2['NuE'];    NuDis  = py_h2['NuDis']

# ------------------------------------------------------------------ #
#  Load IDL outputs  (IDL sav files use lowercase keys)
# ------------------------------------------------------------------ #
idl_h  = readsav('examples/DIII-D/IDL_output/diiid_test.KN1D_H')
idl_h2 = readsav('examples/DIII-D/IDL_output/diiid_test.KN1D_H2')

i_xH  = idl_h['xh'];   i_nH = idl_h['nh'];   i_TH = idl_h['th']
i_GamxH  = idl_h['gammaxh'];   i_Sion  = idl_h['sion']
i_Lyman  = idl_h['lyman'];     i_Balmer = idl_h['balmer']
i_qxH    = idl_h['qxh_total']; i_SRecomb = idl_h['srecomb']
i_SideWallH = idl_h['sidewallh']

i_xH2 = idl_h2['xh2'];  i_nH2 = idl_h2['nh2'];   i_TH2 = idl_h2['th2']
i_GamxH2 = idl_h2['gammaxh2']; i_nHP = idl_h2['nhp'];  i_THP = idl_h2['thp']
i_SH  = idl_h2['sh'];   i_SP  = idl_h2['sp'];     i_qxH2 = idl_h2['qxh2_total']
i_NuLoss = idl_h2['nue'];      i_NuDis = idl_h2['nudis']

# xlimiter for reference line
xlimiter = float(inp['xlimiter'])

# ------------------------------------------------------------------ #
#  Plot helpers
# ------------------------------------------------------------------ #
C_PY  = '#1f77b4'   # blue  — Python
C_IDL = '#d62728'   # red   — IDL
LW = 1.5

def add_limiter(ax):
    ax.axvline(xlimiter, color='k', lw=0.8, ls='--', alpha=0.5, label='limiter')

def legend(ax):
    ax.legend(fontsize=7, loc='best')

fig, axes = plt.subplots(3, 3, figsize=(14, 11))
fig.suptitle('KN1D DIII-D: Python vs IDL comparison', fontsize=13, y=1.01)

# ------------------------------------------------------------------ #
#  (0,0) Density profiles
# ------------------------------------------------------------------ #
ax = axes[0, 0]
ax.semilogy(xA,  n,    color='gray',  lw=LW, ls=':',  label='$n_e$ (input)')
ax.semilogy(xH,  nH,   color=C_PY,   lw=LW,           label='$n_H$ Python')
ax.semilogy(i_xH, i_nH, color=C_IDL, lw=LW, ls='--',  label='$n_H$ IDL')
ax.semilogy(xH2, nH2,  color=C_PY,   lw=LW, ls='-.',   label='$n_{H_2}$ Python')
ax.semilogy(i_xH2, i_nH2, color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='$n_{H_2}$ IDL')
ax.semilogy(xH2, nHP,  color='purple', lw=LW,           label='$n_{H^+}$ Python')
ax.semilogy(i_xH2, i_nHP, color='orange', lw=LW, ls='--', label='$n_{H^+}$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('m$^{-3}$')
ax.set_title('Density profiles')
legend(ax)

# ------------------------------------------------------------------ #
#  (0,1) Temperature profiles
# ------------------------------------------------------------------ #
ax = axes[0, 1]
ax.semilogy(xA,  Ti,   color='gray',  lw=LW, ls=':',   label='$T_i$ (input)')
ax.semilogy(xA,  Te,   color='gray',  lw=LW, ls='--',  label='$T_e$ (input)')
ax.semilogy(xH,  TH,   color=C_PY,   lw=LW,            label='$T_H$ Python')
ax.semilogy(i_xH, i_TH, color=C_IDL, lw=LW, ls='--',   label='$T_H$ IDL')
ax.semilogy(xH2, TH2,  color=C_PY,   lw=LW, ls='-.',    label='$T_{H_2}$ Python')
ax.semilogy(i_xH2, i_TH2, color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='$T_{H_2}$ IDL')
ax.semilogy(xH2, THP,  color='purple', lw=LW,            label='$T_{H^+}$ Python')
ax.semilogy(i_xH2, i_THP, color='orange', lw=LW, ls='--', label='$T_{H^+}$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('eV')
ax.set_title('Temperature profiles')
legend(ax)

# ------------------------------------------------------------------ #
#  (0,2) Particle fluxes
# ------------------------------------------------------------------ #
ax = axes[0, 2]
f = 1e-21
ax.plot(xH2, 2*GamxH2*f,   color=C_PY,  lw=LW,           label='$2\\Gamma_{H_2}$ Python')
ax.plot(i_xH2, 2*i_GamxH2*f, color=C_IDL, lw=LW, ls='--', label='$2\\Gamma_{H_2}$ IDL')
ax.plot(xH,  GamxH*f,      color=C_PY,  lw=LW, ls='-.',   label='$\\Gamma_H$ Python')
ax.plot(i_xH, i_GamxH*f,   color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='$\\Gamma_H$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('$10^{21}$ m$^{-2}$ s$^{-1}$')
ax.set_title('Particle fluxes')
legend(ax)

# ------------------------------------------------------------------ #
#  (1,0) Ionisation source
# ------------------------------------------------------------------ #
ax = axes[1, 0]
ax.semilogy(xH,  Sion,   color=C_PY,  lw=LW,           label='$S_{ion}$ Python')
ax.semilogy(i_xH, i_Sion, color=C_IDL, lw=LW, ls='--',  label='$S_{ion}$ IDL')
ax.semilogy(xH,  SideWallH, color=C_PY, lw=LW, ls='-.',  label='SideWallH Python')
ax.semilogy(i_xH, i_SideWallH, color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='SideWallH IDL')
ax.semilogy(xH, np.abs(SRecomb)+1, color=C_PY, lw=LW, ls=(0,(5,2)), label='$S_{recomb}$ Python')
ax.semilogy(i_xH, np.abs(i_SRecomb)+1, color=C_IDL, lw=LW, ls=':', label='$S_{recomb}$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('m$^{-3}$ s$^{-1}$')
ax.set_title('Sources/Sinks')
legend(ax)

# ------------------------------------------------------------------ #
#  (1,1) Molecular sources
# ------------------------------------------------------------------ #
ax = axes[1, 1]
pos = lambda a: np.where(a > 0, a, np.nan)
ax.semilogy(xH2, pos(SH), color=C_PY,  lw=LW,           label='$S_H$ (from H2) Python')
ax.semilogy(i_xH2, pos(i_SH), color=C_IDL, lw=LW, ls='--', label='$S_H$ (from H2) IDL')
ax.semilogy(xH2, pos(SP), color=C_PY,  lw=LW, ls='-.',   label='$S_p$ (from H2) Python')
ax.semilogy(i_xH2, pos(i_SP), color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='$S_p$ (from H2) IDL')
ax.semilogy(xH2, pos(NuLoss*nHP), color=C_PY, lw=LW, ls=(0,(5,2)), label='NuLoss·$n_{H^+}$ Python')
ax.semilogy(i_xH2, pos(i_NuLoss*i_nHP), color=C_IDL, lw=LW, ls=':', label='NuLoss·$n_{H^+}$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('m$^{-3}$ s$^{-1}$')
ax.set_title('H2 sources/losses')
legend(ax)

# ------------------------------------------------------------------ #
#  (1,2) Heat fluxes
# ------------------------------------------------------------------ #
ax = axes[1, 2]
fkw = 1e-3
ax.plot(xH,  qxH*fkw,   color=C_PY,  lw=LW,           label='$q_{xH}$ Python')
ax.plot(i_xH, i_qxH*fkw, color=C_IDL, lw=LW, ls='--',  label='$q_{xH}$ IDL')
ax.plot(xH2, qxH2*fkw,  color=C_PY,  lw=LW, ls='-.',   label='$q_{xH_2}$ Python')
ax.plot(i_xH2, i_qxH2*fkw, color=C_IDL, lw=LW, ls=(0,(3,1,1,1)), label='$q_{xH_2}$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('kW m$^{-2}$')
ax.set_title('Heat fluxes')
legend(ax)

# ------------------------------------------------------------------ #
#  (2,0) Lyman-alpha emissivity
# ------------------------------------------------------------------ #
ax = axes[2, 0]
ax.plot(xH,  Lyman*1e-3,   color=C_PY,  lw=LW,           label='$L_\\alpha$ Python')
ax.plot(i_xH, i_Lyman*1e-3, color=C_IDL, lw=LW, ls='--',  label='$L_\\alpha$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('kW m$^{-3}$')
ax.set_title('Lyman-$\\alpha$ emissivity')
legend(ax)

# ------------------------------------------------------------------ #
#  (2,1) Balmer-alpha emissivity
# ------------------------------------------------------------------ #
ax = axes[2, 1]
ax.plot(xH,  Balmer*1e-3,   color=C_PY,  lw=LW,           label='$H_\\alpha$ Python')
ax.plot(i_xH, i_Balmer*1e-3, color=C_IDL, lw=LW, ls='--',  label='$H_\\alpha$ IDL')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('kW m$^{-3}$')
ax.set_title('Balmer-$\\alpha$ emissivity')
legend(ax)

# ------------------------------------------------------------------ #
#  (2,2) nH relative difference
# ------------------------------------------------------------------ #
ax = axes[2, 2]
nH_idl_interp = np.interp(xH, i_xH, i_nH)
reldiff = (nH - nH_idl_interp) / np.maximum(nH_idl_interp, 1e10) * 100
ax.plot(xH, reldiff, color='k', lw=LW)
ax.axhline(0, color='gray', lw=0.8, ls='--')
add_limiter(ax)
ax.set_xlabel('x (m)'); ax.set_ylabel('(%)')
ax.set_title('$n_H$ relative difference\n(Python $-$ IDL) / IDL')

fig.tight_layout()
out = 'examples/DIII-D/diiid_test_comparison.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
