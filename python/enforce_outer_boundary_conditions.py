import numpy as np

from util import find_local_numbering

def enforce_outer_boundary_conditions(G, mesh, cells):
    empty_array = np.zeros(0, dtype=np.int64)
    if G.bc == "Dirichlet":
        # Select indices for Dirichlet boundary condition
        mesh.d = np.unique(np.concatenate([
            mesh.e_w, mesh.e_e, mesh.e_n, mesh.e_s, mesh.e_h, mesh.e_l,
            mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le,
            mesh.e_ls, mesh.e_ln, mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se,
            mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne,
            mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne
        ]))

        # Remove indices for Neumann boundary conditions
        mesh.e_w = empty_array
        mesh.e_e = empty_array
        mesh.e_s = empty_array
        mesh.e_n = empty_array
        mesh.e_l = empty_array
        mesh.e_h = empty_array
    elif G.bc == "Neumann":
        mesh.d = empty_array
    elif G.bc == "Neumann_zero_z":
        # Select indices for Dirichlet boundary condition
        mesh.d = np.unique(np.concatenate([
            mesh.e_w, mesh.e_e, mesh.e_n, mesh.e_s,
            mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le,
            mesh.e_ls, mesh.e_ln, mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se,
            mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne,
            mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne
        ]))

        # Remove indices for Neumann boundary conditions
        mesh.e_w = empty_array
        mesh.e_e = empty_array
        mesh.e_s = empty_array
        mesh.e_n = empty_array
    elif G.bc == "Neumann_zero_yz":
        # Select indices for Dirichlet boundary condition
        mesh.d = np.unique(np.concatenate([
            mesh.e_w, mesh.e_e,
            mesh.e_hw, mesh.e_he, mesh.e_lw, mesh.e_le,
            mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se,
            mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne,
            mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne
        ]))

        # Remove indices for Neumann boundary conditions
        mesh.e_w = empty_array
        mesh.e_e = empty_array

    # Update full list of extracellular mesh points (without Dirichlet and corner points)
    mesh.e_all = np.unique(np.concatenate([
        mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, mesh.e_l,
        mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne,
        mesh.m_hsw, mesh.m_hse, mesh.m_hnw, mesh.m_hne,
        mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn,
        mesh.m_lw, mesh.m_le, mesh.m_ls, mesh.m_ln,
        mesh.m_ne, mesh.m_sw, mesh.m_se, mesh.m_nw,
        mesh.m_w, mesh.m_e, mesh.m_s, mesh.m_n, mesh.m_h, mesh.m_l,
        mesh.gx_he, mesh.gx_le, mesh.gx_se, mesh.gx_ne,
        mesh.gy_hn, mesh.gy_ln, mesh.gy_nw, mesh.gy_ne,
        mesh.gx_lse, mesh.gx_lne, mesh.gx_hse, mesh.gx_hne,
        mesh.gy_lnw, mesh.gy_lne, mesh.gy_hnw, mesh.gy_hne,
        mesh.e
    ]))

    N = G.N

    # Update mesh.v_of_e
    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.v] = 1
    vec_of_e = vec[mesh.e_all]
    mesh.v_of_e = np.flatnonzero(vec_of_e)

    # Update mesh.v_inner_of_e
    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.v_inner] = 1
    vec_of_e = vec[mesh.e_all]
    mesh.v_inner_of_e = np.flatnonzero(vec_of_e)

    # update cell.v_of_e
    cells_v_of_e = find_local_numbering(N, mesh.e_all, [cell.v for cell in cells])
    for i, cell in enumerate(cells):
        cell.v_of_e = cells_v_of_e[i]
    del cells_v_of_e

    return mesh, cells
