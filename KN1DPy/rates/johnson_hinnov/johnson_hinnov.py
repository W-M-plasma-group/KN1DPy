import os
import numpy as np
from numpy.typing import NDArray
from scipy import interpolate

from ...utils import get_local_directory, bs2dr

class Johnson_Hinnov():
    '''
    This class handles various operations involving Johnson-Hinnov rate equations and coefficients 
    [L.C.Johnson and E. Hinnov, Journal of Quantitative Spectroscopy & Radiative Transfer. vol. 13 pp.333-358]

    Attributes
    ----------
        dknot : ndarray
            B-spline knot locations in log electron density (m^-3).
        tknot : ndarray
            B-spline knot locations in log electron temperature (eV).
        order : int
            Order of the tensor-product B-spline basis (4 = cubic).
        logr_bscoef : ndarray
            B-spline coefficients for log Johnson-Hinnov population correction
            factors R(ne, Te) for neutral and ionized hydrogen, levels n = 2-6.
        logs_bscoef : ndarray
            B-spline coefficients for log effective ionization rate coefficient
            S(ne, Te) (m^3 s^-1).
        logalpha_bscoef : ndarray
            B-spline coefficients for log effective recombination coefficient
            alpha(ne, Te) (m^3 s^-1).
        a_lyman : ndarray
            Einstein A coefficients (s^-1) for H Lyman transitions (n = 2-1...16-1).
        a_balmer : ndarray
            Einstein A coefficients [s^-1] for H Balmer transitions (n = 3-2...17-2).

    History
    -------
        Coding by B. LaBombard  6/29/99
        Coefficients from J. Terry's idl code JH_RATES.PRO
    '''

    def __init__(self, create: bool=False):
        '''
        Initializes the Johnson-Hinov Coefficients

        Parameters
        ----------
            create : bool, default = False
                If set, regenerate coefficients
                NOTE This feature is currently not working, use file in github
        '''

        path = get_local_directory(__file__)

        if create or not os.path.exists(path+"/jh_bscoef.npz"):
            self._create_jh_bscoef()
        
        jh_data = np.load(path+"/jh_bscoef.npz")

        self.dknot = jh_data['dknot']
        self.tknot = jh_data['tknot']
        self.order = jh_data['order']
        self.logr_bscoef = jh_data['logr_bscoef']
        self.logs_bscoef = jh_data['logs_bscoef']
        self.logalpha_bscoef = jh_data['logalpha_bscoef']
        self.a_lyman = jh_data['a_lyman']
        self.a_balmer = jh_data['a_balmer']


    def jhr_coef(self, Density: NDArray, Te: NDArray, Ion: int, p: int, no_null: bool=False) -> NDArray:
        '''
        Evaluates the r0(p) and r1(p) coefficients from Johnson-Hinnov tables 1 or 2.

        Parameters
        ----------
            Density	: ndarray
                electron density (=hydrogen ion density) (m^-3)
            Te : ndarray
                electron temperature (eV
            Ion	: int
                - 0=return "recombination" coeffcient, r0(p)
                - 1=return "ionization" coeffcient, r1(p)
            p : int
                hydrogen energy level, p=1 is ground state 
            no_null : bool, default=False
                if true, then rather than generate a NULL value when Density and Te
                are outside the data range, compute the rate based on the min or max
                data range values.
        '''
        
        if np.size(Density) != np.size(Te):
            raise Exception('Number of elements of Density and Te are different!')
        if not np.isscalar(Ion):
            raise Exception('"Ion" must be a scalar')
        if Ion < 0 or Ion > 1:
            raise Exception('"Ion" must 0 or 1')
        if not np.isscalar(p):
            raise Exception('"p" must be a scalar')
        if not (1 < p < 7):
            raise Exception('"p" must in range 1 < p < 7')
        
        result = np.full_like(Density, 1.0e32)
        LDensity = np.log(Density)
        LTe = np.log(Te)

        if no_null:
            LDensity = np.clip(LDensity, self.dknot[0] + 0.001, self.dknot[-1] - 0.001)
            LTe = np.clip(LTe, self.tknot[0] + 0.001, self.tknot[-1] - 0.001)
            ok = np.arange(LDensity.size)
        else: # NOTE Not Tested, Might not be needed?
            for i in range(0, len(Density)):
                if min(self.dknot) <= LDensity[i] <= max(self.dknot) and min(self.tknot) <= LTe[i] <= max(self.tknot):
                    ok = np.append(ok, i)

        if ok.size > 0:
            result[ok] = np.exp(bs2dr(LTe[ok], LDensity[ok], self.order, self.order, self.tknot, self.dknot, self.logr_bscoef.T[:,Ion,p-2]))

        return result 
    

    def jhs_coef(self, Density: NDArray, Te: NDArray, no_null: bool=False) -> NDArray:
        '''
        Evaluates the ionization rate coefficients, S (m^-3 s^-1), from  Johnson-Hinnov table 2 (MKS units).
        
        Parameters
        ----------
            Density	: ndarray
                electron density (=hydrogen ion density) (m^-3)
            Te : ndarray
                electron temperature (eV)
            no_null : bool, default=False
                if true, then rather than generate a NULL value when Density and Te
                are outside the data range, compute the rate based on the min or max
                data range values.
        '''

        if np.size(Density) != np.size(Te):
            raise Exception('Number of elements of Density and Te are different!')
        
        result = np.full_like(Density, 1.0e32)
        LDensity = np.log(Density)
        LTe = np.log(Te)

        if no_null:
            LDensity = np.clip(LDensity, self.dknot[0] + 0.001, self.dknot[-1] - 0.001)
            LTe = np.clip(LTe, self.tknot[0] + 0.001, self.tknot[-1] - 0.001)
            ok = np.arange(LDensity.size)
        else: # NOTE Not Tested, Might not be needed?
            for i in range(0, len(Density)):
                if min(self.dknot) <= LDensity[i] <= max(self.dknot) and min(self.tknot) <= LTe[i] <= max(self.tknot):
                    ok = np.append(ok, i)

        if ok.size > 0: 
            result[ok] = np.exp(bs2dr(LTe[ok], LDensity[ok], self.order, self.order, self.tknot, self.dknot, self.logs_bscoef))

        return result
    

    def jhalpha_coef(self, Density: NDArray, Te: NDArray, no_null: bool=False) -> NDArray:
        '''
        Evaluates the recombination rate coefficient, S (m^-6 s^-1), from  Johnson-Hinnov table 1 (MKS units).
        
        Parameters
        ----------
            Density	: ndarray
                electron density (=hydrogen ion density) (m^-3)
            Te : ndarray
                electron temperature (eV)
            no_null : bool, default=False
                if true, then rather than generate a NULL value when Density and Te
                are outside the data range, compute the rate based on the min or max
                data range values.
        '''

        if np.size(Density) != np.size(Te):
            raise Exception('Number of elements of Density and Te are different!')
        
        result = np.full_like(Density, 1.0e32)
        LDensity = np.log(Density)
        LTe = np.log(Te)

        if no_null:
            LDensity = np.clip(LDensity, self.dknot[0] + 0.001, self.dknot[-1] - 0.001)
            LTe = np.clip(LTe, self.tknot[0] + 0.001, self.tknot[-1] - 0.001)
            ok = np.arange(LDensity.size)
        else: # NOTE Not Tested, Might not be needed?
            for i in range(0, len(Density)):
                if min(self.dknot) <= LDensity[i] <= max(self.dknot) and min(self.tknot) <= LTe[i] <= max(self.tknot):
                    ok = np.append(ok, i)

        if ok.size > 0: 
            result[ok] = np.exp(bs2dr(LTe[ok], LDensity[ok], self.order, self.order, self.tknot, self.dknot, self.logalpha_bscoef))

        return result
    

    def nh_saha(self, Density: NDArray, Te: NDArray, p: int) -> NDArray:
        '''
        Evaluates the Saha equilibrium population density (m^-3)
        for atomic hydrogen 
        
        Parameters
        ----------
            Density : ndarray
                electron density (=hydrogen ion density) (m^-3)
            Te : ndarray
                electron temperature (eV)
            p : int, hydrogen energy level, p=1 is ground state
        '''

        if np.size(Density) != np.size(Te):
            raise Exception('Number of Elements of Density and Te are different!')
        if not np.isscalar(p):
            raise Exception('"p" must be a scalar')
        if p < 0:
            raise Exception('“p” must greater than 0')
        
        result = np.full_like(Density, 1.0e32)
        ok = np.where((0.0 < Density) & (Density < 1.0e32) & (0.0 < Te) & (Te < 1.0e32))[0]

        result[ok] = 3.310E-28*((Density[ok]*p)**2)*np.exp(13.6057 / ((p**2)*Te[ok])) / (Te[ok]**1.5)
        return result


    def lyman_alpha(self, Density: NDArray, Te: NDArray, N0: NDArray, no_null: bool=False) -> NDArray:
        '''
        Computes Lyman-alpha emissivity (watts m^-3) given the local
        electron density, electron temperature, and ground-state
        neutral density.
        
        Method:
            (1) Compute the local n=2 population density using the Johnson-Hinnov
                rate equations and coefficients [L.C.Johnson and E. Hinnov, J. Quant. 
                Spectrosc. Radiat. Transfer. vol. 13 pp.333-358]
            (2) Multiply by the n=2->1 spontaneous emission coefficient
            (3) Convert to watts/m^3
        
        Parameters
        ----------
         	Density : ndarray
                electron density (=hydrogen ion density) (m^-3)
         	Te : ndarray
                electron temperature (eV)
         	N0 : ndarray
                ground state neutral density (m^-3)
        	no_null : bool, default=false
                if true, then rather than generate a NULL value when Density and Te
                are not null but still outside the data range, compute the rate based on the min or max
        		data range values.
        '''
        
        # From Johnson-Hinnov, eq (11):
        # n(2) =  ( r0(2) + r1(2) * n(1) / NHsaha(1) ) * NHsaha(2)
        if np.size(Density) != np.size(Te):
            raise Exception('Number of elements of Density and Te are different!')
        if np.size(Density) != np.size(N0):
            raise Exception(' Number of elements of Density and N0 are different! ')
        
        result = np.full(Density.shape,1.0e32)
        photons = np.full(Density.shape,1.0e32)
        r02 = self.jhr_coef(Density, Te, 0, 2, no_null=no_null)
        r12 = self.jhr_coef(Density, Te, 1, 2, no_null=no_null)
        NHSaha1 = self.nh_saha(Density, Te, 1)
        NHSaha2 = self.nh_saha(Density, Te, 2)

        ok = np.where((0 < N0) & (N0 < 1e32) & (r02 < 1.0e32) & (r12 < 1.0e32) & (NHSaha1 < 1.0e32) & (NHSaha2 < 1.0e32))[0]

        photons[ok] = self.a_lyman[0]*(r02[ok] + (r12[ok]*N0[ok] / NHSaha1[ok]))*NHSaha2[ok]
        result[ok] = 13.6057*0.75*photons[ok]*1.6e-19
        return result


    def balmer_alpha(self, Density : NDArray, Te : NDArray, N0 : NDArray, photons = 0, no_null = 0) -> NDArray:
        '''
        Computes Balmer-alpha emissivity (watts m^-3) given the local
        electron density, electron temperature, and ground-state
        neutral density.

        Method:
            (1) Compute the local n=3 population density using the Johnson-Hinnov
                rate equations and coefficients [L.C.Johnson and E. Hinnov, J. Quant. 
                Spectrosc. Radiat. Transfer. vol. 13 pp.333-358]
            (2) Multiply by the n=3->2 spontaneous emission coefficient
            (3) Convert to watts/m^3
        
        Parameters
        ----------
            Density : ndarray
                electron density (=hydrogen ion density) (m^-3)
            Te : ndarray
                electron temperature (eV
            N0 : ndarray
                ground state neutral density (m^-3)
            No_Null : bool, default=false
                if set, then rather than generate a NULL value when Density and Te
                are not null but still outside the data range, compute the rate based on the min or max
                data range values.
        '''

        # From Johnson-Hinnov, eq (11):
        # n(3) = ( r(0) + r1(3) * n(1) / NHsaha(1) ) * NHsaha(3)

        if np.size(Density) != np.size(Te):
            raise Exception('Number of elements of Density and Te are different!')
        if np.size(Density) != np.size(N0):
            raise Exception(' Number of elements of Density and N0 are different! ')
        
        result = np.full(Density.shape,1.0e32)
        photons = np.full(Density.shape,1.0e32)
        r03 = self.jhr_coef(Density, Te, 0, 3, no_null = no_null)
        r13 = self.jhr_coef(Density, Te, 1, 3, no_null = no_null)
        NHSaha1 = self.nh_saha(Density, Te, 1)
        NHSaha3 = self.nh_saha(Density, Te, 3)

        ok = np.where((0 < N0) & (N0 < 1e32) & (r03 < 1.0e32) & (r13 < 1.0e32) & (NHSaha1 < 1.0e32) & (NHSaha3 < 1.0e32))[0]
        
        photons[ok] = self.a_balmer[0]*(r03[ok] + (r13[ok]*N0[ok] / NHSaha1[ok]))*NHSaha3[ok]
        result[ok] = 13.6057*(0.25 - 1.0/9.0)*photons[ok]*1.6e-19
        return result
    



    def _create_jh_bscoef(self, path): 
        # NOTE This function does not currently produce correct results
        # Use the stored version of the file in github

        print("ERROR: JH_BSCoef creation is not working. Please make sure the 'jh_bscoef.npz' file from the github" \
                "is in the same directory as this file")
        
        return

        DENS = np.array([1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16])
        TEMP = np.array([0.345, 0.69, 1.38, 2.76, 5.52, 11.0, 22.1, 44.1, 88.0, 176.5, 706.0])

        r = np.zeros((7,11,2,5))
        
        r[:,:,0,0] = np.array([ [7.6e-6, 1.1e-5, 1.9e-5, 4.9e-5, 2.4e-4, 2.2e-3, 1.8e-2],
                                [1.5e-3, 1.8e-3, 2.5e-3, 4.5e-3, 1.3e-2, 7.1e-2, 3.7e-1],
                                [2.6e-2, 2.9e-2, 3.5e-2, 4.9e-2, 9.6e-2, 3.2e-1, 7.8e-1],
                                [1.3e-1, 1.4e-1, 1.5e-1, 1.9e-1, 2.8e-1, 6.1e-1, 9.2e-1],
                                [3.6e-1, 3.7e-1, 3.8e-1, 4.2e-1, 5.2e-1, 8.0e-1, 9.6e-1],
                                [6.9e-1, 6.9e-1, 7.0e-1, 7.3e-1, 7.9e-1, 9.2e-1, 9.8e-1],
                                [1.1,    1.1,    1.1,    1.1,    1.1,    1.0,    1.0   ],
                                [1.5,    1.5,    1.5,    1.5,    1.4,    1.1,    1.0   ],
                                [2.0,    2.0,    1.9,    1.9,    1.7,    1.3,    1.0   ],
                                [2.4,    2.4,    2.4,    2.3,    2.1,    1.4,    1.1   ],
                                [3.4,    3.4,    3.3,    3.2,    2.9,    2.0,    1.2   ] ]).T

        r[:,:,1,0] = np.array([ [2.5e-7, 2.5e-6, 2.5e-5, 2.5e-4, 2.5e-3, 2.4e-2, 2.0e-1],
                                [1.9e-7, 1.9e-6, 1.9e-5, 1.9e-4, 1.9e-3, 1.8e-2, 1.0e-1],
                                [1.6e-7, 1.6e-6, 1.6e-5, 1.6e-4, 1.5e-3, 1.1e-2, 3.2e-2],
                                [1.5e-7, 1.5e-6, 1.5e-5, 1.5e-4, 1.3e-3, 7.2e-3, 1.3e-2],
                                [1.6e-7, 1.6e-6, 1.6e-5, 1.5e-4, 1.3e-3, 5.4e-3, 8.0e-3],
                                [1.8e-7, 1.8e-6, 1.8e-5, 1.7e-4, 1.4e-3, 5.1e-3, 7.0e-3],
                                [2.1e-7, 2.1e-6, 2.1e-5, 2.0e-4, 1.6e-3, 5.6e-3, 7.5e-3],
                                [2.3e-7, 2.3e-6, 2.3e-5, 2.2e-4, 1.7e-3, 6.3e-3, 8.7e-3],
                                [2.3e-7, 2.3e-6, 2.3e-5, 2.2e-4, 1.8e-3, 7.0e-3, 1.0e-2],
                                [2.2e-7, 2.2e-6, 2.1e-5, 2.1e-4, 1.7e-3, 7.4e-3, 1.1e-2],
                                [1.6e-7, 1.6e-6, 1.6e-5, 1.6e-4, 1.4e-3, 7.2e-3, 1.3e-2] ]).T

        r[:,:,0,1] = np.array([ [2.2e-3, 3.1e-3, 6.0e-3, 2.2e-2, 1.3e-1, 3.5e-1, 4.2e-1],
                                [2.6e-2, 3.3e-2, 5.0e-2, 1.2e-1, 4.3e-1, 7.2e-1, 8.5e-1],
                                [1.1e-1, 1.3e-1, 1.6e-1, 3.0e-1, 6.8e-1, 8.9e-1, 9.7e-1],
                                [2.7e-1, 2.9e-1, 3.4e-1, 5.0e-1, 8.2e-1, 9.5e-1, 9.9e-1],
                                [4.8e-1, 5.0e-1, 5.4e-1, 6.8e-1, 9.0e-1, 9.8e-1, 1.0   ],
                                [7.3e-1, 7.4e-1, 7.7e-1, 8.5e-1, 9.5e-1, 9.9e-1, 1.0   ],
                                [1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0   ],
                                [1.3,    1.3,    1.3,    1.2,    1.1,    1.0,    1.0   ],
                                [1.6,    1.6,    1.5,    1.4,    1.1,    1.0,    1.0   ],
                                [1.9,    1.9,    1.8,    1.6,    1.2,    1.1,    1.0   ],
                                [2.5,    2.4,    2.4,    2.1,    1.5,    1.1,    1.0   ] ]).T

        r[:,:,1,1] = np.array([ [1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.1e-3, 1.3e-2, 1.1e-1],
                                [8.2e-8, 8.1e-7, 8.0e-6, 7.7e-5, 6.1e-4, 4.5e-3, 4.2e-2],
                                [7.1e-8, 7.0e-7, 6.8e-6, 5.9e-5, 3.3e-4, 1.6e-3, 4.1e-3],
                                [6.8e-8, 6.7e-7, 6.3e-6, 4.9e-5, 2.1e-4, 7.3e-4, 1.2e-3],
                                [7.2e-8, 7.0e-7, 6.5e-6, 4.7e-5, 1.8e-4, 4.9e-4, 6.8e-4],
                                [8.1e-8, 7.8e-7, 7.2e-6, 5.1e-5, 1.9e-4, 4.5e-4, 5.8e-4],
                                [9.1e-8, 8.9e-7, 8.2e-6, 5.8e-5, 2.1e-4, 5.0e-4, 6.4e-4],
                                [9.7e-8, 9.5e-7, 8.8e-6, 6.5e-5, 2.5e-4, 6.0e-4, 7.6e-4],
                                [9.7e-8, 9.4e-7, 8.8e-6, 6.7e-5, 2.7e-4, 6.9e-4, 9.1e-4],
                                [8.9e-8, 8.7e-7, 8.2e-6, 6.5e-5, 2.8e-4, 7.6e-4, 1.0e-3],
                                [6.5e-8, 6.4e-7, 6.1e-6, 5.2e-5, 2.7e-4, 8.0e-4, 1.2e-3] ]).T

        r[:,:,0,2] = np.array([ [1.8e-2, 2.8e-2, 7.3e-2, 3.1e-1, 6.0e-1, 7.4e-1, 7.7e-1],
                                [8.2e-2, 1.1e-1, 2.2e-1, 5.6e-1, 8.3e-1, 9.3e-1, 9.6e-1],
                                [2.0e-1, 2.4e-1, 3.9e-1, 7.4e-1, 9.2e-1, 0.98,   0.99  ],
                                [3.7e-1, 4.1e-1, 5.7e-1, 8.4e-1, 9.6e-1, 0.99,   1.0   ],
                                [5.6e-1, 6.0e-1, 7.2e-1, 9.0e-1, 9.8e-1, 1.0,    1.0   ],
                                [7.7e-1, 7.9e-1, 8.5e-1, 9.5e-1, 9.9e-1, 1.0,    1.0   ],
                                [9.9e-1, 9.9e-1, 9.9e-1, 1.0,    1.0,    1.0,    1.0   ],
                                [1.2,    1.2,    1.1,    1.1,    1.0,    1.0,    1.0   ],
                                [1.4,    1.4,    1.3,    1.1,    1.0,    1.0,    1.0   ],
                                [1.7,    1.6,    1.5,    1.2,    1.1,    1.0,    1.0   ],
                                [2.1,    2.1,    1.9,    1.5,    1.1,    1.0,    1.0   ] ]).T

        r[:,:,1,2] = np.array([ [7.2e-8, 7.1e-7, 6.9e-6, 5.7e-5, 4.8e-4, 5.3e-3, 4.5e-2],
                                [5.9e-8, 5.7e-7, 5.1e-6, 3.1e-5, 1.7e-4, 1.1e-3, 5.9e-3],
                                [5.1e-8, 4.9e-7, 4.0e-6, 1.9e-5, 7.1e-5, 3.0e-4, 7.8e-4],
                                [4.8e-8, 4.5e-7, 3.4e-6, 1.4e-5, 4.2e-5, 1.3e-4, 2.1e-4],
                                [5.0e-8, 4.7e-7, 3.4e-6, 1.3e-5, 3.5e-5, 8.6e-5, 1.2e-4],
                                [5.6e-8, 5.2e-7, 3.7e-6, 1.4e-5, 3.6e-5, 8.1e-5, 1.0e-4],
                                [6.3e-8, 5.9e-7, 4.3e-6, 1.6e-5, 4.3e-5, 9.3e-5, 1.2e-4],
                                [6.7e-8, 6.3e-7, 4.8e-6, 1.9e-5, 5.2e-5, 1.1e-4, 1.4e-4],
                                [6.6e-8, 6.3e-7, 4.9e-6, 2.1e-5, 5.9e-5, 1.3e-4, 1.7e-4],
                                [6.1e-8, 5.8e-7, 4.7e-6, 2.2e-5, 6.3e-5, 1.4e-4, 2.0e-4],
                                [4.4e-8, 4.2e-7, 3.7e-6, 2.0e-5, 6.4e-5, 1.6e-4, 2.5e-4] ]).T
        
        r[:,:,0,3] = np.array([ [5.5e-2, 1.0e-1, 3.3e-1, 6.8e-1, 8.5e-1, 0.9,    0.92  ],
                                [1.5e-1, 2.4e-1, 5.5e-1, 8.4e-1, 9.5e-1, 0.98,   0.99  ],
                                [2.9e-1, 4.0e-1, 7.0e-1, 9.1e-1, 9.8e-1, 0.99,   1.0   ],
                                [4.5e-1, 5.5e-1, 8.0e-1, 9.5e-1, 9.9e-1, 1.0,    1.0   ],
                                [6.2e-1, 7.0e-1, 8.7e-1, 9.7e-1, 9.9e-1, 1.0,    1.0   ],
                                [8.0e-1, 8.4e-1, 9.3e-1, 9.8e-1, 1.0,    1.0,    1.0   ],
                                [9.8e-1, 9.8e-1, 9.9e-1, 1.0,    1.0,    1.0,    1.0   ],
                                [1.2,    1.1,    1.1,    1.0,    1.0,    1.0,    1.0   ],
                                [1.4,    1.3,    1.2,    1.0,    1.0,    1.0,    1.0   ],
                                [1.5,    1.5,    1.3,    1.1,    1.0,    1.0,    1.0   ],
                                [1.9,    1.9,    1.6,    1.2,    1.0,    1.0,    1.0   ] ]).T

        r[:,:,1,3] = np.array([ [6.0e-8, 5.7e-7, 4.4e-6, 2.5e-5, 1.8e-4, 2.0e-3, 1.6e-2],
                                [4.8e-8, 4.4e-7, 2.7e-6, 1.1e-5, 5.0e-5, 3.2e-4, 1.7e-3],
                                [4.1e-8, 3.5e-7, 1.8e-6, 5.9e-6, 1.9e-5, 8.1e-5, 2.0e-4],
                                [3.8e-8, 3.2e-7, 1.5e-6, 4.2e-6, 1.1e-5, 3.4e-5, 5.5e-5],
                                [4.0e-8, 3.2e-7, 1.4e-6, 3.8e-6, 9.4e-6, 2.3e-5, 3.1e-5],
                                [4.4e-8, 3.6e-7, 1.6e-6, 4.3e-6, 1.0e-5, 2.2e-5, 2.8e-5],
                                [5.0e-8, 4.1e-7, 1.9e-6, 5.2e-6, 1.2e-5, 2.5e-5, 3.2e-5],
                                [5.3e-8, 4.5e-7, 2.2e-6, 6.2e-6, 1.5e-5, 3.1e-5, 3.9e-5],
                                [5.2e-8, 4.5e-7, 2.4e-6, 7.1e-6, 1.7e-5, 3.7e-5, 4.8e-5],
                                [4.8e-8, 4.3e-7, 2.4e-6, 7.6e-6, 1.9e-5, 4.3e-5, 5.7e-5],
                                [3.5e-8, 3.2e-7, 2.1e-6, 7.5e-6, 2.0e-5, 4.8e-5, 7.2e-5] ]).T

        r[:,:,0,4] = np.array([ [1.1e-1, 2.7e-1, 6.4e-1, 8.6e-1, 9.4e-1, 0.96,   0.97  ],
                                [2.4e-1, 4.5e-1, 7.9e-1, 9.4e-1, 9.8e-1, 0.99,   1.0   ],
                                [3.8e-1, 6.0e-1, 8.7e-1, 9.7e-1, 9.9e-1, 1.0,    1.0   ],
                                [5.3e-1, 7.2e-1, 9.1e-1, 9.8e-1, 1.0,    1.0,    1.0   ],
                                [6.8e-1, 8.1e-1, 9.4e-1, 9.9e-1, 1.0,    1.0,    1.0   ],
                                [8.2e-1, 9.0e-1, 9.7e-1, 9.9e-1, 1.0,    1.0,    1.0   ],
                                [9.7e-1, 9.9e-1, 1.0,    1.0,    1.0,    1.0,    1.0   ],
                                [1.1,    1.1,    1.0,    1.0,    1.0,    1.0,    1.0   ],
                                [1.3,    1.2,    1.1,    1.0,    1.0,    1.0,    1.0   ],
                                [1.5,    1.3,    1.1,    1.0,    1.0,    1.0,    1.0   ],
                                [1.8,    1.7,    1.3,    1.1,    1.0,    1.0,    1.0   ] ]).T

        r[:,:,1,4] = np.array([ [5.2e-8, 4.3e-7, 2.3e-6, 1.0e-5, 7.2e-5, 7.7e-4, 6.5e-3],
                                [4.1e-8, 3.0e-7, 1.2e-6, 4.0e-6, 1.7e-5, 1.1e-4, 5.9e-4],
                                [3.4e-8, 2.2e-7, 7.7e-7, 2.1e-6, 6.6e-6, 2.7e-5, 6.8e-5],
                                [3.1e-8, 1.9e-7, 6.0e-7, 1.5e-6, 3.8e-6, 1.1e-5, 1.8e-5],
                                [3.2e-8, 1.9e-7, 5.9e-7, 1.4e-6, 3.2e-6, 7.7e-6, 1.0e-5],
                                [3.6e-8, 2.1e-7, 6.7e-7, 1.6e-6, 3.5e-6, 7.5e-6, 9.5e-6],
                                [4.1e-8, 2.5e-7, 8.2e-7, 1.9e-6, 4.3e-6, 8.9e-6, 1.1e-5],
                                [4.4e-8, 2.9e-7, 9.8e-7, 2.4e-6, 5.3e-6, 1.1e-5, 1.4e-5],
                                [4.4e-8, 3.0e-7, 1.1e-6, 2.7e-6, 6.3e-6, 1.3e-5, 1.7e-5],
                                [4.0e-8, 3.0e-7, 1.2e-6, 3.0e-6, 7.0e-6, 1.5e-5, 2.1e-5],
                                [3.0e-8, 2.4e-7, 1.1e-6, 3.1e-6, 7.5e-6, 1.8e-5, 2.6e-5] ]).T
        
        # print("R", r.T)
        # input()


        s = np.array([  [2.1e-26, 3.2e-26, 6.5e-26, 2.1e-25, 1.3e-24, 1.4e-23, 1.2e-22],
                        [1.0e-17, 1.3e-17, 2.0e-17, 4.3e-17, 1.5e-16, 9.4e-16, 5.0e-15],
                        [3.0e-13, 3.4e-13, 4.4e-13, 7.1e-13, 1.7e-12, 6.1e-12, 1.5e-11],
                        [6.7e-11, 7.3e-11, 8.6e-11, 1.1e-10, 2.0e-10, 4.9e-10, 7.6e-10],
                        [1.3e-09, 1.4e-09, 1.5e-09, 1.9e-09, 2.7e-09, 5.0e-09, 6.4e-09],
                        [6.9e-09, 7.2e-09, 7.7e-09, 8.9e-09, 1.2e-08, 1.9e-08, 2.2e-08],
                        [1.8e-08, 1.8e-08, 1.9e-08, 2.1e-08, 2.7e-08, 4.0e-08, 4.5e-08],
                        [2.8e-08, 2.9e-08, 3.0e-08, 3.3e-08, 4.1e-08, 5.8e-08, 6.7e-08],
                        [3.4e-08, 3.5e-08, 3.6e-08, 3.9e-08, 4.8e-08, 6.7e-08, 7.7e-08],
                        [3.4e-08, 3.4e-08, 3.6e-08, 3.9e-08, 4.7e-08, 6.5e-08, 7.7e-08],
                        [2.5e-08, 2.6e-08, 2.6e-08, 2.8e-08, 3.3e-08, 4.6e-08, 5.8e-08] ]).T
        # convert s from cm^3 s^-1 to m^3 s^-1
        s = s*1.0e-6

        # print("S", s.T)

        alpha = np.array([  [1.2e-12, 1.7e-12, 2.9e-12, 7.1e-12, 2.7e-11, 1.6e-10, 1.4e-09],
                            [6.1e-13, 7.3e-13, 1.0e-12, 1.7e-12, 3.9e-12, 1.4e-11, 7.1e-11],
                            [3.3e-13, 3.6e-13, 4.3e-13, 5.7e-13, 9.2e-13, 2.0e-12, 4.8e-12],
                            [1.8e-13, 1.9e-13, 2.1e-13, 2.4e-13, 3.1e-13, 4.8e-13, 7.0e-13],
                            [1.0e-13, 1.0e-13, 1.1e-13, 1.2e-13, 1.3e-13, 1.6e-13, 1.9e-13],
                            [5.6e-14, 5.7e-14, 5.7e-14, 5.9e-14, 6.1e-14, 6.5e-14, 7.2e-14],
                            [3.0e-14, 3.0e-14, 3.0e-14, 3.0e-14, 3.0e-14, 3.0e-14, 3.2e-14],
                            [1.5e-14, 1.5e-14, 1.5e-14, 1.5e-14, 1.5e-14, 1.4e-14, 1.5e-14],
                            [7.3e-15, 7.3e-15, 7.2e-15, 7.1e-15, 6.9e-15, 6.6e-15, 6.7e-15],
                            [3.4e-15, 3.4e-15, 3.3e-15, 3.3e-15, 3.2e-15, 3.0e-15, 3.0e-15],
                            [6.5e-16, 6.5e-16, 6.4e-16, 6.4e-16, 6.2e-16, 5.8e-16, 5.7e-16] ]).T
        # convert alpha from cm^3 s^-1 to m^3 s^-1
        alpha = alpha*1.0e-6

        # print("Alpha", alpha.T)

        #the following are the spontaneous emission coeffs for n = 2 to 1
        #   3 to 1, ... , 16 to 1
        A_lyman = np.array([4.699e8, 5.575e7, 1.278e7, 4.125e6, 1.644e6, 7.568e5, 3.869e5,
                            2.143e5, 1.263e5, 7.834e4, 5.066e4, 3.393e4, 2.341e4, 1.657e4,1.200e4])

        #the following are the spontaneous emission coeffs for n = 3 to 2
        #   4 to 2, ... 17 to 2
        A_balmer = np.array([4.410e7, 8.420e6, 2.530e6, 9.732e5, 4.389e5, 2.215e5, 1.216e5,
                            7.122e4, 4.397e4, 2.830e4, 18288.8, 12249.1, 8451.26, 5981.95, 4332.13])

        # convert to MKS, take natural log
        LogDensity = np.log(DENS*1.0e6)
        LogTe = np.log(TEMP)
        LogR = np.log(r)
        LogS = np.log(s)
        LogAlpha = np.log(alpha)
        # print("LogDensity", LogDensity)
        # print("LogTe", LogTe)
        # print("LogR", LogR.T)
        # print("LogS", LogS.T)
        # print("LogAlpha", LogAlpha.T)
        # input()

        # Loop through ION = 0, 1 and p =2, 6 (i = 0, 4)
        # Fit BSCoef to each 
        order = 4

        #NOTE Some values are close, but they are not quite correct, revisit if necessary
        print('Computing B-Spline coefficients for r0 and r1 values')
        LogR_BSCoef = np.zeros((np.size(LogDensity)*np.size(LogTe), 2, 5))
        # print("LogR_BSCoef", LogR_BSCoef)
        # input()

        for nIon in range(0,2):
            for nP in range(2,7):
                # LogR_Interp = scipy.interpolate.RectBivariateSpline(LogDensity, LogTe, LogR[:,:,nIon,nP-2])
                # LogR_BSCoef[:,nIon,nP-2] = LogR_Interp.get_coeffs()

                # RegularGridInterpolater runs a 
                LogR_Interp = interpolate.RegularGridInterpolator((LogDensity, LogTe), LogR[:,:,nIon,nP-2], method='cubic')
                LogR_BSCoef[:,nIon,nP-2] = LogR_Interp._spline.c.reshape((77), order='F')
                # print("LogR_BSCoef", LogR_BSCoef)
                # input()

        # plt.plot(LogR_BSCoef.T)
        # plt.show()
        # print("LogR_BSCoef", LogR_BSCoef.T)
        # input()

        # Do S and Alpha 

        print('Computing B-Spline coefficients for S and alpha values')
        # LogS_Interp = scipy.interpolate.RectBivariateSpline(LogTe, LogDensity, LogS.T, kx=order, ky=order)
        # LogS_BSCoef = LogS_Interp.get_coeffs()

        LogS_Interp = interpolate.RegularGridInterpolator((LogDensity, LogTe), LogS, method='cubic')
        LogS_BSCoef = LogS_Interp._spline.c.reshape((77), order='F')
        print("LogR_BSCoef", LogS_BSCoef.T)
        input()
        
        LogAlpha_Interp = interpolate.RectBivariateSpline(LogTe, LogDensity, LogAlpha.T, kx=order, ky=order)
        LogAlpha_BSCoef = LogAlpha_Interp.get_coeffs()
        TKnot, DKnot = LogAlpha_Interp.get_knots() # get knot locations

        print('Saving results in file: jh_bscoef.npz')
        np.savez(path+"/jh_bscoef",
                DKnot = DKnot,
                TKnot = TKnot,
                order = order,
                LogR_BSCoef = LogR_BSCoef,
                LogS_BSCoef = LogS_BSCoef,
                LogAlpha_BSCoef = LogAlpha_BSCoef,
                A_Lyman = A_lyman,
                A_Balmer = A_balmer)

        return
