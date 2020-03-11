import numpy as np
import matplotlib.pyplot as plt

from util import round_to_nearest

def plot_solution(U, V, G, mesh, cmap="jet"):
    "Plot the intracellular, extracellular and membrane potential"

    z_idx = round_to_nearest(G.Nz/2)

    # Plot the intracellular potential
    V_reshape = np.zeros(G.N, V.dtype)
    V_reshape[mesh.v] = V

    # Remove extracellular points
    U_tmp = U.copy()
    extracellular_index_sets = [mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse,
            mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn,
            mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw,
            mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h,
            mesh.e_l, mesh.e, mesh.gx_he, mesh.gx_le, mesh.gx_se,
            mesh.gx_ne, mesh.gx_lse, mesh.gx_hse, mesh.gx_lne, mesh.gx_hne,
            mesh.gy_ne, mesh.gy_nw, mesh.gy_ln, mesh.gy_hn, mesh.gy_lnw,
            mesh.gy_hnw, mesh.gy_lne, mesh.gy_hne]
    for index_set in extracellular_index_sets:
        U_tmp[index_set] = np.nan
    if hasattr(mesh, 'd'):
        U_tmp[mesh.d] = np.nan

    # Calculate Ui on the membrane
    U_tmp[mesh.v] = U[mesh.v] + V_reshape[mesh.v]

    fig, axs = plt.subplots(3)

    u_tmp = U_tmp.reshape(G.Nz, G.Nx*G.Ny)
    #z_idx = round_to_nearest(G.Nz/2)
    u = u_tmp[z_idx, :].reshape(G.Ny, G.Nx)
    x_1d = np.arange(G.Nx)*G.dx
    y_1d = np.arange(G.Ny)*G.dy
    x, y = np.meshgrid(1e4*x_1d, 1e4*y_1d)
    pc = axs[0].pcolor(x,y,(u), cmap=cmap)
    cbar = fig.colorbar(pc, ax=axs[0])
    cbar.ax.set_ylabel('mV')

    axs[0].set_xlim([0, G.Lx*1e4])
    axs[0].set_ylim([0, G.Ly*1e4])
    axs[0].set_ylabel('y ($\mu$m)')
    axs[0].set_title('Intracellular potential')


    # Plot the extracellular potential
    U_tmp = np.nan*U
    U_tmp[mesh.e_all] = U[mesh.e_all]
    if hasattr(mesh, 'd'):
        U_tmp[mesh.d] = 0

    u_tmp = U_tmp.reshape(G.Nz, G.Nx*G.Ny)
    #z_idx = round_to_nearest((G.Nz-1)/2)
    u = u_tmp[z_idx, :].reshape(G.Ny, G.Nx)
    x_1d = np.arange(G.Nx)*G.dx
    y_1d = np.arange(G.Ny)*G.dy
    x, y = np.meshgrid(1e4*x_1d, 1e4*y_1d)
    pc = axs[1].pcolor(x,y,(u), cmap=cmap)
    cbar = fig.colorbar(pc, ax=axs[1])
    cbar.ax.set_ylabel('mV')

    axs[1].set_xlim([0, G.Lx*1e4])
    axs[1].set_ylim([0, G.Ly*1e4])
    axs[1].set_ylabel('y ($\mu$m)')
    axs[1].set_title('Extracellular potential')


    # Plot the membrane potential
    U_tmp = np.zeros(G.N, dtype=U.dtype)
    U_tmp[:] = np.nan
    U_tmp[mesh.v] = V

    u_tmp = U_tmp.reshape(G.Nz, G.Nx*G.Ny)
    #z_idx = round_to_nearest((G.Nz)/2)
    u = u_tmp[z_idx, :].reshape(G.Ny, G.Nx)
    x_1d = np.arange(G.Nx)*G.dx
    y_1d = np.arange(G.Ny)*G.dy
    x, y = np.meshgrid(1e4*x_1d, 1e4*y_1d)

    pc = axs[2].pcolor(x,y,(u), cmap=cmap)
    cbar = fig.colorbar(pc, ax=axs[2])
    cbar.ax.set_ylabel('mV')

    axs[2].set_xlim([0, G.Lx*1e4])
    axs[2].set_ylim([0, G.Ly*1e4])
    axs[2].set_xlabel('x ($\mu$m)')
    axs[2].set_ylabel('y ($\mu$m)')
    axs[2].set_title('Membrane potential')

    plt.tight_layout()
    plt.show()
