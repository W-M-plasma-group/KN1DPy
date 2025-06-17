import numpy as np
from create_jh_bscoef import create_jh_bscoef
from jhr_coef import jhr_coef
from nh_saha import nh_saha
from jhs_coef import jhs_coef
from jhalpha_coef import jhalpha_coef
import os.path


def N0_from_balmer_alpha(B_Alpha, Density, Te , Source = 0, Ionization = 0,Recombination = 0, create = 0, g=None):
    #________________________________________________________________________________
# Input:
#  	B_alpha 	- fltarr, local Balmer-alpha emissivity (watts m^-3)
#  	Density		- fltarr, electron density (=hydrogen ion density) (m^-3)
#  	Te		- fltarr, electron temperature (eV)
#
# Keywords:
#	Source       - returns the local net ionization source rate m^-3 s^-1
#	Ionization   - returns the local ionization rate m^-3 s^-1
#	Recombination- returns the local recombination rate m^-3 s^-1
#	create	- if set, then create bi-cubic spline coefficients for
#		  interpolation of r0(p) r1(p) and save them in the
#		  default save set. 
#________________________________________________________________________________
# History:
#    Coding by B. LaBombard  6/29/99
#    Coefficients from J. Terry's idl code JH_RATES.PRO

# variables in JH_coef common block - this is only temporary bc we havent finished discussing common blocks
    Dknot = g.JH_Coef_DKnot
    Tknot = g.JH_Coef_TKnot
    LogR_BSCoef=g.JH_Coef_LogR_BSCoef
    LogS_BSCoef=g.JH_Coef_LogS_BSCoef
    LogAlpha_BSCoef=g.JH_Coef_LogAlpha_BSCoef
    A_Lyman=g.JH_Coef_A_Lyman
    A_Balmer=g.JH_Coef_A_Balmer
    if create or not os.path.exists('jh_bscoef.npz'):
        create_jh_bscoef()
    if LogR_BSCoef is None:
        # this is where old data is restored 
        s=np.load('jh_bscoef.npz')
        Dknot=s['DKnot']
        Tknot=s['TKnot']
        order=s['order']
        LogR_BSCoef=s['LogR_BSCoef']
        LogS_BSCoef=s['LogS_BSCoef']
        LogAlpha_BSCoef=s['LogAlpha_BSCoef']
        A_Lyman=s['A_Lyman']
        A_Balmer=s['A_Balmer']
        # update global vars JH_coef common block
        g.JH_Coef_DKnot=Dknot
        g.JH_Coef_TKnot=Tknot
        g.JH_Coef_order=order
        g.JH_Coef_LogR_BSCoef=LogR_BSCoef
        g.JH_Coef_LogS_BSCoef=LogS_BSCoef
        g.JH_Coef_LogAlpha_BSCoef=LogAlpha_BSCoef
        g.JH_Coef_A_Lyman=A_Lyman
        g.JH_Coef_A_Balmer=A_Balmer

    # From Johnson-Hinnov, eq(11):
    #   n(3) =  ( r0(3) + r1(3) * n(1) / NHsaha(1) ) * NHsaha(3)
    # Inverting: 
    #   n(1) =  ( n(3) / NHsaha(3) - r0(3) )* NHsaha(1) / r1(3)
    if np.size(Density) != np.size(Te):
        raise Exception('Number of elements of Density and Te are different!')
    if np.size(Density) != np.size(B_Alpha):
        raise Exception(' Number of elements of Density and B_Alpha are different! ')
    N0 = Density ; N0[:] = 1.0e32
    n3 = N0
    Source = N0
    Ionization = N0
    Recombination = N0
    r03 = jhr_coef(Density, Te, 0, 3, g=g)
    r13 = jhr_coef(Density, Te, 1, 3, g=g)
    NHSaha1 = nh_saha(Density, Te, 1)
    NHSaha3 = nh_saha(Density, Te, 3)
    S = jhs_coef(Density, Te, g=g)
    Alpha = jhalpha_coef(Density, Te, g=g)
    for i in range(0, np.size(Density)):
        if 0 < B_Alpha[i] < 1e32 and r03[i] < 1.0e32 and r13[i] < 1.0e32 and \
            NHSaha1 < 1.0e32 and NHSaha3 < 1.0e32 and 0 <  S < 1.0e32 and 0 < Alpha < 1.0e32:
            ok = np.append(ok, i)
    count = np.size(ok)
    if count > 0:
        n3[ok] = B_Alpha[ok] / ( A_Balmer[0] * 13.6057 * (0.25 - 1.0 / 9.0) * 1.6e-19)
        N0[ok] = ( n3[ok] / NHSaha3[ok] - r03[ok] ) * NHSaha1[ok] / r13[ok]
        # Evaluate ionization, recombination and net source from Johnson-Hinnov equation (12)
        # Ionization = n(1)*Density*S
        # Recombination = Density^2*Alpha
        # Source = Ionization - Recombination
        Ionization[ok] = N0[ok] * (Density[ok] * S[ok])
        Recombination[ok] = Density[ok] * (Density[ok] * Alpha[ok])
        Source[ok] = Ionization[ok] - Recombination[ok]
    return N0


