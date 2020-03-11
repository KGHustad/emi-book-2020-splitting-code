import numpy as np

from index_mesh import IndexMesh
from combine_cells import combine_cells
from util import find_local_numbering

class EMIMesh(IndexMesh):
    def __init__(self, Nx, Ny, Nz, N):
        super().__init__(Nx, Ny, Nz, N)



def set_up_mesh(G):
    "Set up vectors containing the indices of the different types of nodes"


    # Read geometry
    Nx = G.Nx
    Ny = G.Ny
    Nz = G.Nz
    N = G.N

    mesh = EMIMesh(Nx, Ny, Nz, N)

    # x-direction ranges from west (w) to east (e)
    # y-direction ranges from south (s) to north (n)
    # z-direction ranges from low (l) to high (h)

    # Set up indices for outer boundary
    # Set up indices for outer boundary (corners)
    mesh.e_lsw = np.array([mesh[   0,    0,    0]])
    mesh.e_lse = np.array([mesh[Nx-1,    0,    0]])
    mesh.e_lnw = np.array([mesh[   0, Ny-1,    0]])
    mesh.e_lne = np.array([mesh[Nx-1, Ny-1,    0]])
    mesh.e_hsw = np.array([mesh[   0,    0, Nz-1]])
    mesh.e_hse = np.array([mesh[Nx-1,    0, Nz-1]])
    mesh.e_hnw = np.array([mesh[   0, Ny-1, Nz-1]])
    mesh.e_hne = np.array([mesh[Nx-1, Ny-1, Nz-1]])


    # Set up indices for outer boundary (lines)
    # x direction
    mesh.e_ls = mesh[ 1:Nx-1,     0,      0]
    mesh.e_ln = mesh[ 1:Nx-1,  Ny-1,      0]
    mesh.e_hs = mesh[ 1:Nx-1,     0,   Nz-1]
    mesh.e_hn = mesh[ 1:Nx-1,  Ny-1,   Nz-1]

    # y direction
    mesh.e_lw = mesh[     0, 1:Ny-1,      0]
    mesh.e_le = mesh[  Nx-1, 1:Ny-1,      0]
    mesh.e_hw = mesh[     0, 1:Ny-1,   Nz-1]
    mesh.e_he = mesh[  Nx-1, 1:Ny-1,   Nz-1]

    # z direction
    mesh.e_sw = mesh[     0,      0, 1:Nz-1]
    mesh.e_se = mesh[  Nx-1,      0, 1:Nz-1]
    mesh.e_nw = mesh[     0,   Ny-1, 1:Nz-1]
    mesh.e_ne = mesh[  Nx-1,   Ny-1, 1:Nz-1]


    # Set up indices for outer boundary (sides)
    mesh.e_l = mesh[1:-1, 1:-1, 0]
    mesh.e_h = mesh[1:-1, 1:-1, Nz-1]
    mesh.e_s = mesh[1:-1, 0, 1:-1]
    mesh.e_n = mesh[1:-1, Ny-1, 1:-1]
    mesh.e_w = mesh[0, 1:-1, 1:-1]
    mesh.e_e = mesh[Nx-1, 1:-1, 1:-1]


    # Set up indices for the cell
    mesh, cells = combine_cells(G, mesh)
    print("Combined cells")

    # Set up indices for the remaining extracellular cells
    vec = np.ones(N)

    index_sets_to_exclude_from_e = (
        mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne,
        mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne,
        mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn,
        mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln,
        mesh.e_ne, mesh.e_sw, mesh.e_se, mesh.e_nw,
        mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, mesh.e_l,
        mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne,
        mesh.m_hsw, mesh.m_hse, mesh.m_hnw, mesh.m_hne,
        mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn,
        mesh.m_lw, mesh.m_le, mesh.m_ls, mesh.m_ln,
        mesh.m_ne, mesh.m_sw, mesh.m_se, mesh.m_nw,
        mesh.m_w, mesh.m_e, mesh.m_s, mesh.m_n, mesh.m_h, mesh.m_l,
        mesh.i_all
    )
    for index_set in index_sets_to_exclude_from_e:
        vec[index_set] = 0
    mesh.e = np.flatnonzero(vec)


    # Set up indices for the membrane potential
    mesh.v = np.concatenate([
        mesh.m_lsw,
        mesh.m_lse,
        mesh.m_lnw,
        mesh.m_lne,
        mesh.m_hsw,
        mesh.m_hse,
        mesh.m_hnw,
        mesh.m_hne,
        mesh.m_hw,
        mesh.m_he,
        mesh.m_hs,
        mesh.m_hn,
        mesh.m_lw,
        mesh.m_le,
        mesh.m_ls,
        mesh.m_ln,
        mesh.m_ne,
        mesh.m_sw,
        mesh.m_se,
        mesh.m_nw,
        mesh.m_w,
        mesh.m_e,
        mesh.m_s,
        mesh.m_n,
        mesh.m_h,
        mesh.m_l
    ])
    mesh.v.sort()

    # Set up indices for the gap junction
    mesh.w = np.unique(np.concatenate([mesh.gx_e, mesh.gy_n]))
    mesh.w.sort()

    # Set up indices for the nodes not in the intracellular domain
    mesh.not_i = np.concatenate([
        mesh.e_lsw,
        mesh.e_lse,
        mesh.e_lnw,
        mesh.e_lne,
        mesh.e_hsw,
        mesh.e_hse,
        mesh.e_hnw,
        mesh.e_hne,
        mesh.e_hw,
        mesh.e_he,
        mesh.e_hs,
        mesh.e_hn,
        mesh.e_lw,
        mesh.e_le,
        mesh.e_ls,
        mesh.e_ln,
        mesh.e_ne,
        mesh.e_sw,
        mesh.e_se,
        mesh.e_nw,
        mesh.e_w,
        mesh.e_e,
        mesh.e_s,
        mesh.e_n,
        mesh.e_h,
        mesh.e_l,
        mesh.e,
        mesh.gx_he,
        mesh.gx_le,
        mesh.gx_se,
        mesh.gx_ne,
        mesh.gx_lse,
        mesh.gx_hse,
        mesh.gx_lne,
        mesh.gx_hne,
        mesh.gy_ne,
        mesh.gy_nw,
        mesh.gy_ln,
        mesh.gy_hn,
        mesh.gy_lnw,
        mesh.gy_hnw,
        mesh.gy_lne,
        mesh.gy_hne
    ])
    mesh.not_i.sort()

    # Set up indices for the extracellular domain
    mesh.e_all = np.concatenate([
        mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne,
        mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne,
        mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn,
        mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln,
        mesh.e_ne, mesh.e_sw, mesh.e_se, mesh.e_nw,
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
    ])
    mesh.e_all.sort()

    mesh.i_all = np.unique(np.concatenate([
        mesh.v, mesh.w, mesh.i, mesh.i_lsw, mesh.i_lse, mesh.i_lnw, mesh.i_lne, mesh.i_hsw, mesh.i_hse,
        mesh.i_hnw, mesh.i_hne, mesh.i_hw, mesh.i_he, mesh.i_hs, mesh.i_hn,
        mesh.i_lw, mesh.i_le, mesh.i_ls, mesh.i_ln, mesh.i_ne, mesh.i_sw,
        mesh.i_se, mesh.i_nw, mesh.i_w, mesh.i_e, mesh.i_s, mesh.i_n, mesh.i_h,
        mesh.i_l, mesh.gx_e, mesh.gy_n, mesh.gxi_w, mesh.gxi_hw, mesh.gxi_lw,
        mesh.gxi_sw, mesh.gxi_nw, mesh.gyi_hs, mesh.gyi_ls, mesh.gyi_sw,
        mesh.gyi_se, mesh.gyi_s, mesh.gxi_lsw, mesh.gxi_lnw, mesh.gxi_hsw,
        mesh.gxi_hnw, mesh.gyi_lsw, mesh.gyi_lse, mesh.gyi_hsw, mesh.gyi_hse
    ]))
    mesh.i_all.sort()

    mesh.i_inner_all = np.unique(np.concatenate([
        mesh.i,
        mesh.i_lsw, mesh.i_lse, mesh.i_lnw, mesh.i_lne,
        mesh.i_hsw, mesh.i_hse, mesh.i_hnw, mesh.i_hne,
        mesh.i_hw, mesh.i_he, mesh.i_hs, mesh.i_hn,
        mesh.i_lw, mesh.i_le, mesh.i_ls, mesh.i_ln,
        mesh.i_ne, mesh.i_sw, mesh.i_se, mesh.i_nw,
        mesh.i_w, mesh.i_e, mesh.i_s, mesh.i_n, mesh.i_l, mesh.i_h
    ]))
    mesh.i_inner_all.sort()


    # local numbering in mesh.v
    # cells[n].v should not overlap with cells[k].v, where k != n
    cells_v_idx = find_local_numbering(N, mesh.v, [cell.v for cell in cells])
    for i, cell in enumerate(cells):
        cell.v_idx = cells_v_idx[i]
    del cells_v_idx


    # local numbering in e_all
    cells_v_of_e = find_local_numbering(N, mesh.e_all, [cell.v for cell in cells])
    for i, cell in enumerate(cells):
        cell.v_of_e = cells_v_of_e[i]
    del cells_v_of_e


    # local numbering in i_all
    cells_v_of_i = find_local_numbering(N, mesh.i_all, [cell.v for cell in cells])
    cells_gxw_of_i = find_local_numbering(N, mesh.i_all, [cell.gx_w for cell in cells])
    cells_gxe_of_i = find_local_numbering(N, mesh.i_all, [cell.gx_e for cell in cells])
    cells_gys_of_i = find_local_numbering(N, mesh.i_all, [cell.gy_s for cell in cells])
    cells_gyn_of_i = find_local_numbering(N, mesh.i_all, [cell.gy_n for cell in cells])
    for i, cell in enumerate(cells):
        cell.v_of_i = cells_v_of_i[i]
        cell.gxw_of_i = cells_gxw_of_i[i]
        cell.gxe_of_i = cells_gxe_of_i[i]
        cell.gys_of_i = cells_gys_of_i[i]
        cell.gyn_of_i = cells_gyn_of_i[i]
    del cells_v_of_i
    del cells_gxw_of_i
    del cells_gxe_of_i
    del cells_gys_of_i
    del cells_gyn_of_i

    # cells[n].c_all can overlap with cells[k].c_all, where k != n
    # solution: treat cell grid as checkered board so that direct neighbours are split in different sets.
    # Swipe over the white first, then the black.
    #cells[n].c_of_i

    cell_indices_white = []
    cell_indices_black = []
    for j in range(G.num_cells_y):
        for i in range(G.num_cells_x):
            cell_index = j*G.num_cells_x + i
            if (i + j) & 1 == 1:
                cell_indices_black.append(cell_index)
            else:
                cell_indices_white.append(cell_index)

    cells_c_of_i_white = find_local_numbering(N, mesh.i_all, [cells[i].c_all for i in cell_indices_white])
    cells_c_of_i_black = find_local_numbering(N, mesh.i_all, [cells[i].c_all for i in cell_indices_black])
    for i, cell_index in enumerate(cell_indices_white):
        cells[cell_index].c_of_i = cells_c_of_i_white[i]
    for i, cell_index in enumerate(cell_indices_black):
        cells[cell_index].c_of_i = cells_c_of_i_black[i]

    del cell_indices_white
    del cell_indices_black
    del cells_c_of_i_white
    del cells_c_of_i_black


    # local numbering in w
    cells_gxw_of_w = find_local_numbering(N, mesh.w, [cell.gx_w for cell in cells])
    cells_gxe_of_w = find_local_numbering(N, mesh.w, [cell.gx_e for cell in cells])
    cells_gys_of_w = find_local_numbering(N, mesh.w, [cell.gy_s for cell in cells])
    cells_gyn_of_w = find_local_numbering(N, mesh.w, [cell.gy_n for cell in cells])
    for i, cell in enumerate(cells):
        cell.gxw_of_w = cells_gxw_of_w[i]
        cell.gxe_of_w = cells_gxe_of_w[i]
        cell.gys_of_w = cells_gys_of_w[i]
        cell.gyn_of_w = cells_gyn_of_w[i]
    del cells_gxw_of_w
    del cells_gxe_of_w
    del cells_gys_of_w
    del cells_gyn_of_w


    # Inner membrane points
    mesh.v_inner = np.concatenate([mesh.m_l, mesh.m_h, mesh.m_s, mesh.m_n, mesh.m_e, mesh.m_w])
    mesh.v_inner.sort()

    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.v_inner] = 1
    vec = vec[mesh.v]
    mesh.v_inner_of_v = np.flatnonzero(vec)


    # Set up way to extract membrane points from the intracellular and extracellular domains
    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.v] = 1
    vec_of_i = vec[mesh.i_all]
    vec_of_e = vec[mesh.e_all]
    mesh.v_of_i = np.flatnonzero(vec_of_i)
    mesh.v_of_e = np.flatnonzero(vec_of_e)

    # v inner
    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.v_inner] = 1
    vec_of_i = vec[mesh.i_all]
    vec_of_e = vec[mesh.e_all]
    mesh.v_inner_of_i = np.flatnonzero(vec_of_i)
    mesh.v_inner_of_e = np.flatnonzero(vec_of_e)


    # Set up way to extract intercalated disc points from the intracellular domain
    vec = np.zeros(N, dtype=np.int64)
    vec[mesh.w] = 1
    vec_of_i = vec[mesh.i_all]
    mesh.w_of_i = np.flatnonzero(vec_of_i)


    return mesh, cells
