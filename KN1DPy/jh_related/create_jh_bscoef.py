import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

#    Creates a save set string Bi-cubic spline interpolation 
# coefficients for parameters in Johnson - Hinov rate equations. 
#Gwendolyn Galleher 

#NOTE This File is not currently working - 9/20/25

def create_jh_bscoef():
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
            LogR_Interp = scipy.interpolate.RegularGridInterpolator((LogDensity, LogTe), LogR[:,:,nIon,nP-2], method='cubic')
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

    LogS_Interp = scipy.interpolate.RegularGridInterpolator((LogDensity, LogTe), LogS, method='cubic')
    LogS_BSCoef = LogS_Interp._spline.c.reshape((77), order='F')
    print("LogR_BSCoef", LogS_BSCoef.T)
    input()
    
    LogAlpha_Interp = scipy.interpolate.RectBivariateSpline(LogTe, LogDensity, LogAlpha.T, kx=order, ky=order)
    LogAlpha_BSCoef = LogAlpha_Interp.get_coeffs()
    TKnot, DKnot = LogAlpha_Interp.get_knots() # get knot locations

    print('Saving results in file: jh_bscoef.npz')
    np.savez("jh_bscoef",
             DKnot = DKnot,
             TKnot = TKnot,
             order = order,
             LogR_BSCoef = LogR_BSCoef,
             LogS_BSCoef = LogS_BSCoef,
             LogAlpha_BSCoef = LogAlpha_BSCoef,
             A_Lyman = A_lyman,
             A_Balmer = A_balmer)

    return
    
