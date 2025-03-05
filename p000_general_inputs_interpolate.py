import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def original_data():
    n_data = np.array([1.E10, 1.E11, 1.E12, 1.E13, 1.E14, 1.E15, 1.E16]) * 1.0e6  # Data converted from cm^-3 to m^-3
    t_data = np.array([0.345, 0.69, 1.38, 2.76, 5.52, 11.0, 22.1, 44.1, 88.0, 176.5, 706.])

    s_data = np.zeros((len(n_data), len(t_data)), dtype=np.float64)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    s_data[:,0]=np.array([2.1e-26,3.2e-26,6.5e-26,2.1e-25,1.3e-24,1.4e-23,1.2e-22])
    s_data[:,1]=np.array([1.0e-17,1.3e-17,2.0e-17,4.3e-17,1.5e-16,9.4e-16,5.0e-15])
    s_data[:,2]=np.array([3.0e-13,3.4e-13,4.4e-13,7.1e-13,1.7e-12,6.1e-12,1.5e-11])
    s_data[:,3]=np.array([6.7e-11,7.3e-11,8.6e-11,1.1e-10,2.0e-10,4.9e-10,7.6e-10])
    s_data[:,4]=np.array([1.3e-9,1.4e-9,1.5e-9,1.9e-9,2.7e-9,5.0e-9,6.4e-9])
    s_data[:,5]=np.array([6.9e-9,7.2e-9,7.7e-9,8.9e-9,1.2e-8,1.9e-8,2.2e-8])
    s_data[:,6]=np.array([1.8e-8,1.8e-8,1.9e-8,2.1e-8,2.7e-8,4.0e-8,4.5e-8])
    s_data[:,7]=np.array([2.8e-8,2.9e-8,3.0e-8,3.3e-8,4.1e-8,5.8e-8,6.7e-8])
    s_data[:,8]=np.array([3.4e-8,3.5e-8,3.6e-8,3.9e-8,4.8e-8,6.7e-8,7.7e-8])
    s_data[:,9]=np.array([3.4e-8,3.4e-8,3.6e-8,3.9e-8,4.7e-8,6.5e-8,7.7e-8])
    s_data[:,10]=np.array([2.5e-8,2.6e-8,2.6e-8,2.8e-8,3.3e-8,4.6e-8,5.8e-8])
    # ------------------------------------------------------------------------------------------
    # Convert data from cm^{3} s^{-1} to m^{3} s^{-1}
    s_data = s_data * 1.0e-6   
    # ------------------------------------------------------------------------------------------
    alpha_data=np.zeros((len(n_data), len(t_data)), dtype=np.float64)
    alpha_data[:,0] = np.array([1.2e-12,1.7e-12,2.9e-12,7.1e-12,2.7e-11,1.6e-10,1.4e-9])
    alpha_data[:,1] = np.array([6.1e-13,7.3e-13,1.0e-12,1.7e-12,3.9e-12,1.4e-11,7.1e-11])
    alpha_data[:,2] = np.array([3.3e-13,3.6e-13,4.3e-13,5.7e-13,9.2e-13,2.0e-12,4.8e-12])
    alpha_data[:,3] = np.array([1.8e-13,1.9e-13,2.1e-13,2.4e-13,3.1e-13,4.8e-13,7.0e-13])
    alpha_data[:,4] = np.array([1.0e-13,1.0e-13,1.1e-13,1.2e-13,1.3e-13,1.6e-13,1.9e-13])
    alpha_data[:,5] = np.array([5.6e-14,5.7e-14,5.7e-14,5.9e-14,6.1e-14,6.5e-14,7.2e-14])
    alpha_data[:,6] = np.array([3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.2e-14])
    alpha_data[:,7] = np.array([1.5e-14,1.5e-14,1.5e-14,1.5e-14,1.5e-14,1.4e-14,1.5e-14])
    alpha_data[:,8] = np.array([7.3e-15,7.3e-15,7.2e-15,7.1e-15,6.9e-15,6.6e-15,6.7e-15])
    alpha_data[:,9] = np.array([3.4e-15,3.4e-15,3.3e-15,3.3e-15,3.2e-15,3.0e-15,3.0e-15])
    alpha_data[:,10] = np.array([6.5e-16,6.5e-16,6.4e-16,6.4e-16,6.2e-16,5.8e-16,5.7e-16])
    # ------------------------------------------------------------------------------------------
    # Convert data from cm^{3} s^{-1} to m^{3} s^{-1}
    alpha_data = alpha_data*1.0e-6

    log_n_data  = np.log(n_data)
    log_t_data  = np.log(t_data)
    log_s_data  = np.log(s_data)
    log_alpha_data = np.log(alpha_data)

    return log_n_data,log_t_data,log_s_data,log_alpha_data


if __name__ == "__main__":
    log_n_data, log_t_data, log_s_data, log_alpha_data = original_data()

    log_n_grid, log_t_grid = np.meshgrid(log_n_data, log_t_data, indexing="ij")

    # Plot
    fig = plt.figure(figsize=(12, 8))

    # 1st subplot (log(S))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot_surface(log_n_grid, log_t_grid, log_s_data, cmap='viridis', alpha=0.7, edgecolor='k')
    ax1.scatter(log_n_grid, log_t_grid, log_s_data, color='r')
    ax1.set_xlabel(r'log(n) [$m^{-3}$]')
    ax1.set_ylabel(r'log(T) [$eV$]')
    ax1.set_zlabel(r'log(S) [$m^3 .s^{-1}$]')
    ax1.set_title(r'log($T_e$) $v_s$ log($n_e$) $v_s$ log($S$)')

    # 2nd subplot (log(Alpha))
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_surface(log_n_grid, log_t_grid, log_alpha_data, cmap='viridis', alpha=0.7, edgecolor='k')
    ax2.scatter(log_n_grid, log_t_grid, log_alpha_data, color='r')
    ax2.set_xlabel(r'log(n) [$m^{-3}$]')
    ax2.set_ylabel(r'log(T) [$eV$]')
    ax2.set_zlabel(r'log($\alpha$) [$m^3 .s^{-1}$]')
    ax2.set_title(r'log($T_e$) $v_s$ log($n_e$) $v_s$ log($\alpha$)')

    # show
    plt.tight_layout()
    plt.show()


