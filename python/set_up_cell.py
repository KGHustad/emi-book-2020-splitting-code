import numpy as np

from util import round_to_nearest

class Cell(object):
    pass

def set_up_cell(G, mesh):

    # Load parameters
    #idx_x = idx % G.num_cells_x
    #idx_y = idx // G.num_cells_x
    idx_x = 0
    idx_y = 0
    dx = G.dx
    dy = G.dy
    dz = G.dz
    Nx = G.Nx
    Ny = G.Ny
    Nz = G.Nz
    cell_start_x = G.cell_start_x[idx_x]
    cell_start_y = G.cell_start_y[idx_y]
    cell_start_z = G.cell_start_z
    cell_length_x = G.cell_length_x
    cell_length_y = G.cell_length_y
    cell_length_z = G.cell_length_z
    j_length_x = G.j_length_x
    j_length_y = G.j_length_y
    j_width_y_x = G.j_width_y_x # width (in x-direction) of gap junction in y-direction
    j_width_y_z = G.j_width_y_z # width (in z-direction) of gap junction in y-direction
    j_width_x_y = G.j_width_x_y # width (in y-direction) of gap junction in x-direction
    j_width_x_z = G.j_width_x_z # width (in z-direction) of gap junction in x-direction
    """
    idx_x = rem(idx-1, G.num_cells_x) + 1;
    idx_y = floor((idx-1)/G.num_cells_x) + 1;
    dx = G.dx;
    dy = G.dy;
    dz = G.dz;
    Nx = G.Nx;
    Ny = G.Ny;
    Nz = G.Nz;
    cell_start_x = G.cell_start_x(idx_x);
    cell_start_y = G.cell_start_y(idx_y);
    cell_start_z = G.cell_start_z;
    cell_length_x = G.cell_length_x;
    cell_length_y = G.cell_length_y;
    cell_length_z = G.cell_length_z;
    j_length_x = G.j_length_x;
    j_length_y = G.j_length_y;
    j_width_y_x = G.j_width_y_x; # width (in x-direction) of gap junction in y-direction
    j_width_y_z = G.j_width_y_z; # width (in z-direction) of gap junction in y-direction
    j_width_x_y = G.j_width_x_y; # width (in y-direction) of gap junction in x-direction
    j_width_x_z = G.j_width_x_z; # width (in z-direction) of gap junction in x-direction
    """

    cell = Cell()

    # Define cell indices
    cell_start_ix = int(round_to_nearest(cell_start_x/dx))
    cell_start_iy = int(round_to_nearest(cell_start_y/dy))
    cell_start_iz = int(round_to_nearest(cell_start_z/dz))
    cell_end_ix = cell_start_ix + int(round_to_nearest(cell_length_x/dx))
    cell_end_iy = cell_start_iy + int(round_to_nearest(cell_length_y/dy))
    cell_end_iz = cell_start_iz + int(round_to_nearest(cell_length_z/dz))
    #print("cell start (%d, %d, %d)  cell end (%d, %d, %d)" % (cell_start_ix, cell_start_iy, cell_start_iz, cell_end_ix, cell_end_iy, cell_end_iz))
    """
    cell_start_ix = round(cell_start_x/dx)+1;
    cell_start_iy = round(cell_start_y/dy)+1;
    cell_start_iz = round(cell_start_z/dz)+1;
    cell_end_ix = cell_start_ix + round(cell_length_x/dx);
    cell_end_iy = cell_start_iy + round(cell_length_y/dy);
    cell_end_iz = cell_start_iz + round(cell_length_z/dz);
    """

    # Number of intervals for the junction parts
    j_length_x_n = int(round_to_nearest(j_length_x/dx))
    j_length_y_n = int(round_to_nearest(j_length_y/dy))
    j_width_x_n = int(round_to_nearest(j_width_y_x/dx)) # width (in x-direction) of gap junction in y-direction
    j_width_y_z_n = int(round_to_nearest(j_width_y_z/dz)) # width (in z-direction) of gap junction in y-direction
    j_width_y_n = int(round_to_nearest(j_width_x_y/dy)) # width (in y-direction) of gap junction in x-direction
    j_width_x_z_n = int(round_to_nearest(j_width_x_z/dz)) # width (in z-direction) of gap junction in x-direction
    """
    j_length_x_n = round(j_length_x/dx);
    j_length_y_n = round(j_length_y/dy);
    j_width_x_n = round(j_width_y_x/dx); # width (in x-direction) of gap junction in y-direction
    j_width_y_z_n = round(j_width_y_z/dz); # width (in z-direction) of gap junction in y-direction
    j_width_y_n = round(j_width_x_y/dy); # width (in y-direction) of gap junction in x-direction
    j_width_x_z_n = round(j_width_x_z/dz); # width (in z-direction) of gap junction in x-direction
    """

    # Number of intervals for the non-junction part
    n_cell_x = int(round_to_nearest(cell_length_x/dx)) # width (in x-direction) of full cell (including gap junctions)
    n_cell_x_m = n_cell_x - j_width_x_n # width (in x-direction) of full cell minus the width of north/south (y-directional) gap junction)
    n_cell_x_half = int(round_to_nearest(n_cell_x_m/2))
    n_cell_y = int(round_to_nearest(cell_length_y/dy))
    n_cell_y_m = n_cell_y - j_width_y_n
    n_cell_y_half = int(round_to_nearest(n_cell_y_m/2))
    n_cell_z = int(round_to_nearest(cell_length_z/dz))
    n_cell_xz_m = n_cell_z - j_width_x_z_n
    n_cell_xz_half = int(round_to_nearest(n_cell_xz_m/2))
    n_cell_yz_m = n_cell_z - j_width_y_z_n
    n_cell_yz_half = int(round_to_nearest(n_cell_yz_m/2))
    """
    n_cell_x = round(cell_length_x/dx);
    n_cell_x_m = n_cell_x - j_width_x_n;
    n_cell_x_half = round(n_cell_x_m/2);
    n_cell_y = round(cell_length_y/dy);
    n_cell_y_m = n_cell_y - j_width_y_n;
    n_cell_y_half = round(n_cell_y_m/2);
    n_cell_z = round(cell_length_z/dz);
    n_cell_xz_m = n_cell_z - j_width_x_z_n;
    n_cell_xz_half = round(n_cell_xz_m/2);
    n_cell_yz_m = n_cell_z - j_width_y_z_n;
    n_cell_yz_half = round(n_cell_yz_m/2);
    """

    # Points for start and end of junction parts
    cut_start_ix = cell_start_ix + n_cell_x_half
    cut_end_ix = cut_start_ix + j_width_x_n
    cut_start_iy = cell_start_iy + n_cell_y_half
    cut_end_iy = cut_start_iy + j_width_y_n

    cut_start_ixz = cell_start_iz + n_cell_xz_half
    cut_end_ixz = cut_start_ixz + j_width_x_z_n
    cut_start_iyz = cell_start_iz + n_cell_yz_half
    cut_end_iyz = cut_start_iyz + j_width_y_z_n
    #print("cut_start_ixz: %d   cut_start_ixz: %d" % (cut_start_ixz, cut_start_iyz))

    j_start_ix = cell_start_ix - j_length_x_n
    j_end_ix = cell_end_ix + j_length_x_n
    j_start_iy = cell_start_iy - j_length_y_n
    j_end_iy = cell_end_iy + j_length_y_n
    #print("j start (%d, %d)  j end (%d, %d)\n" % (j_start_ix, j_start_iy, j_end_ix, j_end_iy))

    cell.j_start_ix = j_start_ix
    cell.j_start_iy = j_start_iy
    cell.j_end_ix = j_end_ix
    cell.j_end_iy = j_end_iy
    """
    cut_start_ix = cell_start_ix + n_cell_x_half;
    cut_end_ix = cut_start_ix + j_width_x_n;
    cut_start_iy = cell_start_iy + n_cell_y_half;
    cut_end_iy = cut_start_iy + j_width_y_n;

    cut_start_ixz = cell_start_iz + n_cell_xz_half;
    cut_end_ixz = cut_start_ixz + j_width_x_z_n;
    cut_start_iyz = cell_start_iz + n_cell_yz_half;
    cut_end_iyz = cut_start_iyz + j_width_y_z_n;

    j_start_ix = cell_start_ix - j_length_x_n;
    j_end_ix = cell_end_ix + j_length_x_n;
    j_start_iy = cell_start_iy - j_length_y_n;
    j_end_iy = cell_end_iy + j_length_y_n;

    cell.j_start_ix = j_start_ix;
    """

    if hasattr(G, 'simple_cell') and G.simple_cell:
        raise NotImplementedError("Simple cells are not supported")

    # TODO: Fix the rest of this function

    # Set up indices for the membrane points
    cell.m_lsw = ((cell_start_iz)*Ny + (cell_start_iy))*Nx + cell_start_ix
    cell.m_lse = ((cell_start_iz)*Ny +(cell_start_iy))*Nx + cell_end_ix
    cell.m_lnw = ((cell_start_iz)*Ny + (cell_end_iy))*Nx + cell_start_ix
    cell.m_lne = ((cell_start_iz)*Ny + (cell_end_iy))*Nx + cell_end_ix
    cell.m_hsw = ((cell_end_iz)*Ny + (cell_start_iy))*Nx + cell_start_ix
    cell.m_hse = ((cell_end_iz)*Ny +(cell_start_iy))*Nx + cell_end_ix
    cell.m_hnw = ((cell_end_iz)*Ny + (cell_end_iy))*Nx + cell_start_ix
    cell.m_hne = ((cell_end_iz)*Ny + (cell_end_iy))*Nx + cell_end_ix
    """
    # Set up indices for the membrane points
    cell.m_lsw = ((cell_start_iz-1)*Ny + (cell_start_iy-1))*Nx + cell_start_ix;
    cell.m_lse = ((cell_start_iz-1)*Ny +(cell_start_iy-1))*Nx + cell_end_ix;
    cell.m_lnw = ((cell_start_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_start_ix;
    cell.m_lne = ((cell_start_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_end_ix;
    cell.m_hsw = ((cell_end_iz-1)*Ny + (cell_start_iy-1))*Nx + cell_start_ix;
    cell.m_hse = ((cell_end_iz-1)*Ny +(cell_start_iy-1))*Nx + cell_end_ix;
    cell.m_hnw = ((cell_end_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_start_ix;
    cell.m_hne = ((cell_end_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_end_ix;
    """

    # Set up indices for the membrane lines
    # Lines in x direction
    m_ls_main = mesh[cell_start_ix+1:cell_end_ix, cell_start_iy, cell_start_iz] #np.arange(cell.m_lsw + 1, cell.m_lse)
    m_ls_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_start_iy, cut_start_ixz]
    m_ls_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_start_iy, cut_start_ixz]
    cell.m_ls = np.concatenate([m_ls_main, m_ls_west_apx, m_ls_east_apx])

    m_hs_main = mesh[cell_start_ix+1:cell_end_ix, cell_start_iy, cell_end_iz] #np.arange(cell.m_lsw + 1, cell.m_lse)
    m_hs_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_start_iy, cut_end_ixz]
    m_hs_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_start_iy, cut_end_ixz]
    cell.m_hs = np.concatenate([m_hs_main, m_hs_west_apx, m_hs_east_apx])

    m_ln_main = mesh[cell_start_ix+1:cell_end_ix, cell_end_iy, cell_start_iz] #np.arange(cell.m_lsw + 1, cell.m_lse)
    m_ln_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_end_iy, cut_start_ixz]
    m_ln_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_end_iy, cut_start_ixz]
    cell.m_ln = np.concatenate([m_ln_main, m_ln_west_apx, m_ln_east_apx])

    m_hn_main = mesh[cell_start_ix+1:cell_end_ix, cell_end_iy, cell_end_iz] #np.arange(cell.m_lsw + 1, cell.m_lse)
    m_hn_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_end_iy, cut_end_ixz]
    m_hn_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_end_iy, cut_end_ixz]
    cell.m_hn = np.concatenate([m_hn_main, m_hn_west_apx, m_hn_east_apx])

    # Lines in y direction
    m_lw_main = mesh[cell_start_ix, cell_start_iy+1:cell_end_iy, cell_start_iz]
    m_lw_south_apx = mesh[cut_start_ix, j_start_iy+1:cell_start_iy, cut_start_iyz]
    m_lw_north_apx = mesh[cut_start_ix, cell_end_iy+1:j_end_iy:, cut_start_iyz]
    cell.m_lw = np.concatenate([m_lw_main, m_lw_south_apx, m_lw_north_apx])

    m_hw_main = mesh[cell_start_ix, cell_start_iy+1:cell_end_iy, cell_end_iz]
    m_hw_south_apx = mesh[cut_start_ix, j_start_iy+1:cell_start_iy, cut_end_iyz]
    m_hw_north_apx = mesh[cut_start_ix, cell_end_iy+1:j_end_iy:, cut_end_iyz]
    cell.m_hw = np.concatenate([m_hw_main, m_hw_south_apx, m_hw_north_apx])

    m_le_main = mesh[cell_end_ix, cell_start_iy+1:cell_end_iy, cell_start_iz]
    m_le_south_apx = mesh[cut_end_ix, j_start_iy+1:cell_start_iy, cut_start_iyz]
    m_le_north_apx = mesh[cut_end_ix, cell_end_iy+1:j_end_iy:, cut_start_iyz]
    cell.m_le = np.concatenate([m_le_main, m_le_south_apx, m_le_north_apx])

    m_he_main = mesh[cell_end_ix, cell_start_iy+1:cell_end_iy, cell_end_iz]
    m_he_south_apx = mesh[cut_end_ix, j_start_iy+1:cell_start_iy, cut_end_iyz]
    m_he_north_apx = mesh[cut_end_ix, cell_end_iy+1:j_end_iy:, cut_end_iyz]
    cell.m_he = np.concatenate([m_he_main, m_he_south_apx, m_he_north_apx])

    # Lines in z direction
    cell.m_sw = mesh[cell_start_ix, cell_start_iy, cell_start_iz+1:cell_end_iz]
    cell.m_se = mesh[cell_end_ix, cell_start_iy, cell_start_iz+1:cell_end_iz]
    cell.m_nw = mesh[cell_start_ix, cell_end_iy, cell_start_iz+1:cell_end_iz]
    cell.m_ne = mesh[cell_end_ix, cell_end_iy, cell_start_iz+1:cell_end_iz]


    # Set up indices for the membrane planes
    # Five upper and lower planes
    m_l_main = mesh[cell_start_ix+1:cell_end_ix, cell_start_iy+1:cell_end_iy, cell_start_iz]
    m_h_main = mesh[cell_start_ix+1:cell_end_ix, cell_start_iy+1:cell_end_iy, cell_end_iz]
    m_l_west = mesh[j_start_ix+1:cell_start_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz]
    m_h_west = mesh[j_start_ix+1:cell_start_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz]
    m_l_east = mesh[cell_end_ix+1:j_end_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz]
    m_h_east = mesh[cell_end_ix+1:j_end_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz]
    m_l_south = mesh[cut_start_ix+1:cut_end_ix, j_start_iy+1:cell_start_iy, cut_start_iyz]
    m_h_south = mesh[cut_start_ix+1:cut_end_ix, j_start_iy+1:cell_start_iy, cut_end_iyz]
    m_l_north = mesh[cut_start_ix+1:cut_end_ix, cell_end_iy+1:j_end_iy, cut_start_iyz]
    m_h_north = mesh[cut_start_ix+1:cut_end_ix, cell_end_iy+1:j_end_iy, cut_end_iyz]

    cell.m_l = np.concatenate([
        m_l_main,
        m_l_south,
        m_l_west,
        m_l_east,
        m_l_north,
    ])
    cell.m_h = np.concatenate([
        m_h_main,
        m_h_south,
        m_h_west,
        m_h_east,
        m_h_north,
    ])
    cell.m_l.sort()
    cell.m_h.sort()


    # Membrane planes on junction sides
    # junctions in x-direction
    m_s_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz]
    m_n_west_apx = mesh[j_start_ix+1:cell_start_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz]
    m_s_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz]
    m_n_east_apx = mesh[cell_end_ix+1:j_end_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz]

    # junctions in y-direction
    m_w_south_apx = mesh[cut_start_ix, j_start_iy+1:cell_start_iy, cut_start_iyz+1:cut_end_iyz]
    m_e_south_apx = mesh[cut_end_ix, j_start_iy+1:cell_start_iy, cut_start_iyz+1:cut_end_iyz]
    m_w_north_apx = mesh[cut_start_ix, cell_end_iy+1:j_end_iy, cut_start_iyz+1:cut_end_iyz]
    m_e_north_apx = mesh[cut_end_ix, cell_end_iy+1:j_end_iy, cut_start_iyz+1:cut_end_iyz]


    # Membrane planes on cell sides
    # We need to carve out the middle rectangular part that is covered by the appendix
    #
    # FGH
    # D*E
    # ABC
    #

    # xz planes for south and north appendices
    m_s_main = np.concatenate([
        mesh[cell_start_ix+1:cut_start_ix, cell_start_iy, cell_start_iz+1:cut_start_iyz], # A
        mesh[cut_start_ix:cut_end_ix+1, cell_start_iy, cell_start_iz+1:cut_start_iyz], # B
        mesh[cut_end_ix+1:cell_end_ix, cell_start_iy, cell_start_iz+1:cut_start_iyz], # C
        mesh[cell_start_ix+1:cut_start_ix, cell_start_iy, cut_start_iyz:cut_end_iyz+1], # D
        mesh[cut_end_ix+1:cell_end_ix, cell_start_iy, cut_start_iyz:cut_end_iyz+1], # E
        mesh[cell_start_ix+1:cut_start_ix, cell_start_iy, cut_end_iyz+1:cell_end_iz], # F
        mesh[cut_start_ix:cut_end_ix+1, cell_start_iy, cut_end_iyz+1:cell_end_iz], # G
        mesh[cut_end_ix+1:cell_end_ix, cell_start_iy, cut_end_iyz+1:cell_end_iz], # H
    ])
    m_s_main.sort()
    m_n_main = m_s_main + (mesh[0, cell_end_iy, 0] - mesh[0, cell_start_iy, 0])

    # yz planes for west and east appendices
    m_w_main = np.concatenate([
        mesh[cell_start_ix, cell_start_iy+1:cut_start_iy, cell_start_iz+1:cut_start_iyz], # A
        mesh[cell_start_ix, cut_start_iy:cut_end_iy+1, cell_start_iz+1:cut_start_iyz], # B
        mesh[cell_start_ix, cut_end_iy+1:cell_end_iy, cell_start_iz+1:cut_start_iyz], # C
        mesh[cell_start_ix, cell_start_iy+1:cut_start_iy, cut_start_iyz:cut_end_iyz+1], # D
        mesh[cell_start_ix, cut_end_iy+1:cell_end_iy, cut_start_iyz:cut_end_iyz+1], # E
        mesh[cell_start_ix, cell_start_iy+1:cut_start_iy, cut_end_iyz+1:cell_end_iz], # F
        mesh[cell_start_ix, cut_start_iy:cut_end_iy+1, cut_end_iyz+1:cell_end_iz], # G
        mesh[cell_start_ix, cut_end_iy+1:cell_end_iy, cut_end_iyz+1:cell_end_iz], # H
    ])
    m_w_main.sort()
    m_e_main = m_w_main + (mesh[cell_end_ix, 0, 0] - mesh[cell_start_ix, 0, 0])

    cell.m_s = np.concatenate([m_s_main, m_s_west_apx, m_s_east_apx])
    cell.m_n = np.concatenate([m_n_main, m_n_west_apx, m_n_east_apx])
    cell.m_w = np.concatenate([m_w_main, m_w_south_apx, m_w_north_apx])
    cell.m_e = np.concatenate([m_e_main, m_e_south_apx, m_e_north_apx])
    cell.m_s.sort()
    cell.m_n.sort()
    cell.m_w.sort()
    cell.m_e.sort()


    # Set up indices for the membrane points between cells
    cell.gx_lsw = mesh[j_start_ix, cut_start_iy, cut_start_ixz]
    cell.gx_lnw = mesh[j_start_ix, cut_end_iy, cut_start_ixz]
    cell.gx_lse = mesh[j_end_ix, cut_start_iy, cut_start_ixz]
    cell.gx_lne = mesh[j_end_ix, cut_end_iy, cut_start_ixz]
    cell.gx_hsw = mesh[j_start_ix, cut_start_iy, cut_end_ixz]
    cell.gx_hnw = mesh[j_start_ix, cut_end_iy, cut_end_ixz]
    cell.gx_hse = mesh[j_end_ix, cut_start_iy, cut_end_ixz]
    cell.gx_hne = mesh[j_end_ix, cut_end_iy, cut_end_ixz]

    cell.gy_lsw = mesh[cut_start_ix, j_start_iy, cut_start_iyz]
    cell.gy_lse = mesh[cut_end_ix, j_start_iy, cut_start_iyz]
    cell.gy_lnw = mesh[cut_start_ix, j_end_iy, cut_start_iyz]
    cell.gy_lne = mesh[cut_end_ix, j_end_iy, cut_start_iyz]
    cell.gy_hsw = mesh[cut_start_ix, j_start_iy, cut_end_iyz]
    cell.gy_hse = mesh[cut_end_ix, j_start_iy, cut_end_iyz]
    cell.gy_hnw = mesh[cut_start_ix, j_end_iy, cut_end_iyz]
    cell.gy_hne = mesh[cut_end_ix, j_end_iy, cut_end_iyz]


    # Set up indices for membrane lines between cells
    # Junctions in x direction
    cell.gx_lw = mesh[j_start_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz]
    cell.gx_hw = mesh[j_start_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz]
    cell.gx_sw = mesh[j_start_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz]
    cell.gx_nw = mesh[j_start_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz]

    cell.gx_le = mesh[j_end_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz]
    cell.gx_he = mesh[j_end_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz]
    cell.gx_se = mesh[j_end_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz]
    cell.gx_ne = mesh[j_end_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz]

    # Junctions in y direction
    cell.gy_ls = mesh[cut_start_ix+1:cut_end_ix, j_start_iy, cut_start_iyz]
    cell.gy_hs = mesh[cut_start_ix+1:cut_end_ix, j_start_iy, cut_end_iyz]
    cell.gy_sw = mesh[cut_start_ix, j_start_iy, cut_start_iyz+1:cut_end_iyz]
    cell.gy_se = mesh[cut_end_ix, j_start_iy, cut_start_iyz+1:cut_end_iyz]

    cell.gy_ln = mesh[cut_start_ix+1:cut_end_ix, j_end_iy, cut_start_iyz]
    cell.gy_hn = mesh[cut_start_ix+1:cut_end_ix, j_end_iy, cut_end_iyz]
    cell.gy_nw = mesh[cut_start_ix, j_end_iy, cut_start_iyz+1:cut_end_iyz]
    cell.gy_ne = mesh[cut_end_ix, j_end_iy, cut_start_iyz+1:cut_end_iyz]


    # Membrane planes between cells
    cell.gy_s = mesh[cut_start_ix+1:cut_end_ix, j_start_iy, cut_start_iyz+1:cut_end_iyz]
    cell.gy_n = mesh[cut_start_ix+1:cut_end_ix, j_end_iy, cut_start_iyz+1:cut_end_iyz]

    cell.gx_w = mesh[j_start_ix, cut_start_iy+1:cut_end_iy, cut_start_iyz+1:cut_end_iyz]
    cell.gx_e = mesh[j_end_ix, cut_start_iy+1:cut_end_iy, cut_start_iyz+1:cut_end_iyz]


    ############## Inner boundary of the cell ##############
    cell.i_lsw = mesh[cell_start_ix+1, cell_start_iy+1, cell_start_iz+1]
    cell.i_lse = mesh[cell_end_ix-1, cell_start_iy+1, cell_start_iz+1]
    cell.i_lnw = mesh[cell_start_ix+1, cell_end_iy-1, cell_start_iz+1]
    cell.i_lne = mesh[cell_end_ix-1, cell_end_iy-1, cell_start_iz+1]
    cell.i_hsw = mesh[cell_start_ix+1, cell_start_iy+1, cell_end_iz-1]
    cell.i_hse = mesh[cell_end_ix-1, cell_start_iy+1, cell_end_iz-1]
    cell.i_hnw = mesh[cell_start_ix+1, cell_end_iy-1, cell_end_iz-1]
    cell.i_hne = mesh[cell_end_ix-1, cell_end_iy-1, cell_end_iz-1]


    # Extra inner points between junction and cell
    cell.i_lsw = np.array((
        cell.i_lsw,
        mesh[cell_start_ix, cut_start_iy, cut_start_ixz], # west appendix
        mesh[cut_start_ix, cell_start_iy, cut_start_iyz], # south appendix
    ))
    cell.i_hsw = np.array((
        cell.i_hsw,
        mesh[cell_start_ix, cut_start_iy, cut_end_ixz], # west appendix
        mesh[cut_start_ix, cell_start_iy, cut_end_iyz], # south appendix
    ))

    cell.i_lnw = np.array((
        cell.i_lnw,
        mesh[cell_start_ix, cut_end_iy, cut_start_ixz], # west appendix
        mesh[cut_start_ix, cell_end_iy, cut_start_iyz], # south appendix
    ))
    cell.i_hnw = np.array((
        cell.i_hnw,
        mesh[cell_start_ix, cut_end_iy, cut_end_ixz], # west appendix
        mesh[cut_start_ix, cell_end_iy, cut_end_iyz], # south appendix
    ))

    cell.i_lse = np.array((
        cell.i_lse,
        mesh[cell_end_ix, cut_start_iy, cut_start_ixz], # west appendix
        mesh[cut_end_ix, cell_start_iy, cut_start_iyz], # south appendix
    ))
    cell.i_hse = np.array((
        cell.i_hse,
        mesh[cell_end_ix, cut_start_iy, cut_end_ixz], # west appendix
        mesh[cut_end_ix, cell_start_iy, cut_end_iyz], # south appendix
    ))

    cell.i_lne = np.array((
        cell.i_lne,
        mesh[cell_end_ix, cut_end_iy, cut_start_ixz], # west appendix
        mesh[cut_end_ix, cell_end_iy, cut_start_iyz], # south appendix
    ))
    cell.i_hne = np.array((
        cell.i_hne,
        mesh[cell_end_ix, cut_end_iy, cut_end_ixz], # west appendix
        mesh[cut_end_ix, cell_end_iy, cut_end_iyz], # south appendix
    ))


    # Set up indices for the inner lines
    # Lines in x direction
    i_ls_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_start_iy+1, cell_start_iz+1]
    i_ls_west_apx = mesh[j_start_ix+1+1:cell_start_ix-1, cut_start_iy+1, cut_start_ixz+1]
    i_ls_east_apx = mesh[cell_end_ix+1+1:j_end_ix-1, cut_start_iy+1, cut_start_ixz+1]
    cell.i_ls = np.concatenate([i_ls_main, i_ls_west_apx, i_ls_east_apx])

    i_hs_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_start_iy+1, cell_end_iz-1]
    i_hs_west_apx = mesh[j_start_ix+1+1:cell_start_ix-1, cut_start_iy+1, cut_end_ixz-1]
    i_hs_east_apx = mesh[cell_end_ix+1+1:j_end_ix-1, cut_start_iy+1, cut_end_ixz-1]
    cell.i_hs = np.concatenate([i_hs_main, i_hs_west_apx, i_hs_east_apx])

    i_ln_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_end_iy-1, cell_start_iz+1]
    i_ln_west_apx = mesh[j_start_ix+1+1:cell_start_ix-1, cut_end_iy-1, cut_start_ixz+1]
    i_ln_east_apx = mesh[cell_end_ix+1+1:j_end_ix-1, cut_end_iy-1, cut_start_ixz+1]
    cell.i_ln = np.concatenate([i_ln_main, i_ln_west_apx, i_ln_east_apx])

    i_hn_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_end_iy-1, cell_end_iz-1]
    i_hn_west_apx = mesh[j_start_ix+1+1:cell_start_ix-1, cut_end_iy-1, cut_end_ixz-1]
    i_hn_east_apx = mesh[cell_end_ix+1+1:j_end_ix-1, cut_end_iy-1, cut_end_ixz-1]
    cell.i_hn = np.concatenate([i_hn_main, i_hn_west_apx, i_hn_east_apx])

    # Lines in y direction
    i_lw_main = mesh[cell_start_ix+1, cell_start_iy+1+1:cell_end_iy-1, cell_start_iz+1]
    i_lw_south_apx = mesh[cut_start_ix+1, j_start_iy+1+1:cell_start_iy-1, cut_start_iyz+1]
    i_lw_north_apx = mesh[cut_start_ix+1, cell_end_iy+1+1:j_end_iy-1:, cut_start_iyz+1]
    cell.i_lw = np.concatenate([i_lw_main, i_lw_south_apx, i_lw_north_apx])

    i_hw_main = mesh[cell_start_ix+1, cell_start_iy+1+1:cell_end_iy-1, cell_end_iz-1]
    i_hw_south_apx = mesh[cut_start_ix+1, j_start_iy+1+1:cell_start_iy-1, cut_end_iyz-1]
    i_hw_north_apx = mesh[cut_start_ix+1, cell_end_iy+1+1:j_end_iy-1:, cut_end_iyz-1]
    cell.i_hw = np.concatenate([i_hw_main, i_hw_south_apx, i_hw_north_apx])

    i_le_main = mesh[cell_end_ix-1, cell_start_iy+1+1:cell_end_iy-1, cell_start_iz+1]
    i_le_south_apx = mesh[cut_end_ix-1, j_start_iy+1+1:cell_start_iy-1, cut_start_iyz+1]
    i_le_north_apx = mesh[cut_end_ix-1, cell_end_iy+1+1:j_end_iy-1:, cut_start_iyz+1]
    cell.i_le = np.concatenate([i_le_main, i_le_south_apx, i_le_north_apx])

    i_he_main = mesh[cell_end_ix-1, cell_start_iy+1+1:cell_end_iy-1, cell_end_iz-1]
    i_he_south_apx = mesh[cut_end_ix-1, j_start_iy+1+1:cell_start_iy-1, cut_end_iyz-1]
    i_he_north_apx = mesh[cut_end_ix-1, cell_end_iy+1+1:j_end_iy-1:, cut_end_iyz-1]
    cell.i_he = np.concatenate([i_he_main, i_he_south_apx, i_he_north_apx])

    # Lines in z direction
    cell.i_sw = mesh[cell_start_ix+1, cell_start_iy+1, cell_start_iz+1+1:cell_end_iz-1]
    cell.i_se = mesh[cell_end_ix-1, cell_start_iy+1, cell_start_iz+1+1:cell_end_iz-1]
    cell.i_nw = mesh[cell_start_ix+1, cell_end_iy-1, cell_start_iz+1+1:cell_end_iz-1]
    cell.i_ne = mesh[cell_end_ix-1, cell_end_iy-1, cell_start_iz+1+1:cell_end_iz-1]


    # Extra inner lines between junction and cell
    # Z direction
    cell.i_sw = np.concatenate((
        cell.i_sw,
        mesh[cell_start_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz],
        mesh[cut_start_ix, cell_start_iy, cut_start_ixz+1:cut_end_ixz],
    ))
    cell.i_nw = np.concatenate((
        cell.i_nw,
        mesh[cell_start_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz],
        mesh[cut_start_ix, cell_end_iy, cut_start_ixz+1:cut_end_ixz],
    ))
    cell.i_se = np.concatenate((
        cell.i_se,
        mesh[cell_end_ix, cut_start_iy, cut_start_ixz+1:cut_end_ixz],
        mesh[cut_end_ix, cell_start_iy, cut_start_ixz+1:cut_end_ixz],
    ))
    cell.i_ne = np.concatenate((
        cell.i_ne,
        mesh[cell_end_ix, cut_end_iy, cut_start_ixz+1:cut_end_ixz],
        mesh[cut_end_ix, cell_end_iy, cut_start_ixz+1:cut_end_ixz],
    ))

    # Y direction
    cell.i_lw = np.concatenate((
        cell.i_lw,
        mesh[cell_start_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz],
    ))
    cell.i_hw = np.concatenate((
        cell.i_hw,
        mesh[cell_start_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz],
    ))
    cell.i_le = np.concatenate((
        cell.i_le,
        mesh[cell_end_ix, cut_start_iy+1:cut_end_iy, cut_start_ixz],
    ))
    cell.i_he = np.concatenate((
        cell.i_he,
        mesh[cell_end_ix, cut_start_iy+1:cut_end_iy, cut_end_ixz],
    ))

    # X direction
    cell.i_ls = np.concatenate((
        cell.i_ls,
        mesh[cut_start_ix+1:cut_end_ix, cell_start_iy, cut_start_iyz],
    ))
    cell.i_hs = np.concatenate((
        cell.i_hs,
        mesh[cut_start_ix+1:cut_end_ix, cell_start_iy, cut_end_iyz],
    ))
    cell.i_ln = np.concatenate((
        cell.i_ln,
        mesh[cut_start_ix+1:cut_end_ix, cell_end_iy, cut_start_iyz],
    ))
    cell.i_hn = np.concatenate((
        cell.i_hn,
        mesh[cut_start_ix+1:cut_end_ix, cell_end_iy, cut_end_iyz],
    ))


    # Set up indices for the inner planes
    # Five upper and lower planes
    i_l_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_start_iy+1+1:cell_end_iy-1, cell_start_iz+1]
    i_h_main = mesh[cell_start_ix+1+1:cell_end_ix-1, cell_start_iy+1+1:cell_end_iy-1, cell_end_iz-1]

    #from IPython import embed
    #embed()
    #exit()
    i_l_west = mesh[j_start_ix+1+1:cell_start_ix, cut_start_iy+1+1:cut_end_iy-1, cut_start_ixz+1]
    i_h_west = mesh[j_start_ix+1+1:cell_start_ix, cut_start_iy+1+1:cut_end_iy-1, cut_end_ixz-1]
    i_l_east = mesh[cell_end_ix+1:j_end_ix-1, cut_start_iy+1+1:cut_end_iy-1, cut_start_ixz+1]
    i_h_east = mesh[cell_end_ix+1:j_end_ix-1, cut_start_iy+1+1:cut_end_iy-1, cut_end_ixz-1]
    i_l_south = mesh[cut_start_ix+1+1:cut_end_ix-1, j_start_iy+1+1:cell_start_iy, cut_start_iyz+1]
    i_h_south = mesh[cut_start_ix+1+1:cut_end_ix-1, j_start_iy+1+1:cell_start_iy, cut_end_iyz-1]
    i_l_north = mesh[cut_start_ix+1+1:cut_end_ix-1, cell_end_iy+1:j_end_iy-1, cut_start_iyz+1]
    i_h_north = mesh[cut_start_ix+1+1:cut_end_ix-1, cell_end_iy+1:j_end_iy-1, cut_end_iyz-1]

    cell.i_l = np.concatenate([
        i_l_main,
        i_l_south,
        i_l_west,
        i_l_east,
        i_l_north,
    ])
    cell.i_h = np.concatenate([
        i_h_main,
        i_h_south,
        i_h_west,
        i_h_east,
        i_h_north,
    ])
    cell.i_l.sort()
    cell.i_h.sort()


    # Inner planes on junction sides
    # junctions in x-direction
    i_s_west_apx = mesh[j_start_ix+1+1:cell_start_ix, cut_start_iy+1, cut_start_ixz+1:cut_end_ixz]
    i_n_west_apx = mesh[j_start_ix+1+1:cell_start_ix, cut_end_iy-1, cut_start_ixz+1:cut_end_ixz]
    i_s_east_apx = mesh[cell_end_ix+1:j_end_ix-1, cut_start_iy+1, cut_start_ixz+1:cut_end_ixz]
    i_n_east_apx = mesh[cell_end_ix+1:j_end_ix-1, cut_end_iy-1, cut_start_ixz+1:cut_end_ixz]

    # junctions in y-direction
    i_w_south_apx = mesh[cut_start_ix+1, j_start_iy+1+1:cell_start_iy, cut_start_iyz+1:cut_end_iyz]
    i_e_south_apx = mesh[cut_end_ix-1, j_start_iy+1+1:cell_start_iy, cut_start_iyz+1:cut_end_iyz]
    i_w_north_apx = mesh[cut_start_ix+1, cell_end_iy+1:j_end_iy-1, cut_start_iyz+1:cut_end_iyz]
    i_e_north_apx = mesh[cut_end_ix-1, cell_end_iy+1:j_end_iy-1, cut_start_iyz+1:cut_end_iyz]


    # Inner planes on cell sides

    # xz planes for south and north appendices
    i_s_main = np.concatenate([
        mesh[cell_start_ix+1+1:cut_start_ix, cell_start_iy+1, cell_start_iz+1+1:cut_start_iyz], # A
        mesh[cut_start_ix:cut_end_ix+1, cell_start_iy+1, cell_start_iz+1+1:cut_start_iyz], # B
        mesh[cut_end_ix+1:cell_end_ix-1, cell_start_iy+1, cell_start_iz+1+1:cut_start_iyz], # C
        mesh[cell_start_ix+1+1:cut_start_ix, cell_start_iy+1, cut_start_iyz:cut_end_iyz+1], # D
        mesh[cut_end_ix+1:cell_end_ix-1, cell_start_iy+1, cut_start_iyz:cut_end_iyz+1], # E
        mesh[cell_start_ix+1+1:cut_start_ix, cell_start_iy+1, cut_end_iyz+1:cell_end_iz-1], # F
        mesh[cut_start_ix:cut_end_ix+1, cell_start_iy+1, cut_end_iyz+1:cell_end_iz-1], # G
        mesh[cut_end_ix+1:cell_end_ix-1, cell_start_iy+1, cut_end_iyz+1:cell_end_iz-1], # H
    ])
    i_s_main.sort()
    i_n_main = i_s_main + mesh[0, cell_end_iy-1, 0] - mesh[0, cell_start_iy+1, 0]

    # yz planes for west and east appendices
    i_w_main = np.concatenate([
        mesh[cell_start_ix+1, cell_start_iy+1+1:cut_start_iy, cell_start_iz+1+1:cut_start_iyz], # A
        mesh[cell_start_ix+1, cut_start_iy:cut_end_iy+1, cell_start_iz+1+1:cut_start_iyz], # B
        mesh[cell_start_ix+1, cut_end_iy+1:cell_end_iy-1, cell_start_iz+1+1:cut_start_iyz], # C
        mesh[cell_start_ix+1, cell_start_iy+1+1:cut_start_iy, cut_start_iyz:cut_end_iyz+1], # D
        mesh[cell_start_ix+1, cut_end_iy+1:cell_end_iy-1, cut_start_iyz:cut_end_iyz+1], # E
        mesh[cell_start_ix+1, cell_start_iy+1+1:cut_start_iy, cut_end_iyz+1:cell_end_iz-1], # F
        mesh[cell_start_ix+1, cut_start_iy:cut_end_iy+1, cut_end_iyz+1:cell_end_iz-1], # G
        mesh[cell_start_ix+1, cut_end_iy+1:cell_end_iy-1, cut_end_iyz+1:cell_end_iz-1], # H
    ])
    i_w_main.sort()
    i_e_main = i_w_main + mesh[cell_end_ix-1, 0, 0] - mesh[cell_start_ix+1, 0, 0]

    cell.i_s = np.concatenate([i_s_main, i_s_west_apx, i_s_east_apx])
    cell.i_n = np.concatenate([i_n_main, i_n_west_apx, i_n_east_apx])
    cell.i_w = np.concatenate([i_w_main, i_w_south_apx, i_w_north_apx])
    cell.i_e = np.concatenate([i_e_main, i_e_south_apx, i_e_north_apx])
    cell.i_s.sort()
    cell.i_n.sort()
    cell.i_w.sort()
    cell.i_e.sort()


    # Set up indices for the membrane points between cells
    cell.gxi_lsw = mesh[j_start_ix+1, cut_start_iy+1, cut_start_ixz+1]
    cell.gxi_lnw = mesh[j_start_ix+1, cut_end_iy-1, cut_start_ixz+1]
    cell.gxi_lse = mesh[j_end_ix-1, cut_start_iy+1, cut_start_ixz+1]
    cell.gxi_lne = mesh[j_end_ix-1, cut_end_iy-1, cut_start_ixz+1]
    cell.gxi_hsw = mesh[j_start_ix+1, cut_start_iy+1, cut_end_ixz-1]
    cell.gxi_hnw = mesh[j_start_ix+1, cut_end_iy-1, cut_end_ixz-1]
    cell.gxi_hse = mesh[j_end_ix-1, cut_start_iy+1, cut_end_ixz-1]
    cell.gxi_hne = mesh[j_end_ix-1, cut_end_iy-1, cut_end_ixz-1]

    cell.gyi_lsw = mesh[cut_start_ix+1, j_start_iy+1, cut_start_iyz+1]
    cell.gyi_lse = mesh[cut_end_ix-1, j_start_iy+1, cut_start_iyz+1]
    cell.gyi_lnw = mesh[cut_start_ix+1, j_end_iy-1, cut_start_iyz+1]
    cell.gyi_lne = mesh[cut_end_ix-1, j_end_iy-1, cut_start_iyz+1]
    cell.gyi_hsw = mesh[cut_start_ix+1, j_start_iy+1, cut_end_iyz-1]
    cell.gyi_hse = mesh[cut_end_ix-1, j_start_iy+1, cut_end_iyz-1]
    cell.gyi_hnw = mesh[cut_start_ix+1, j_end_iy-1, cut_end_iyz-1]
    cell.gyi_hne = mesh[cut_end_ix-1, j_end_iy-1, cut_end_iyz-1]


    # Set up indices for membrane lines between cells
    cell.gxi_lw = mesh[j_start_ix+1, cut_start_iy+1+1:cut_end_iy-1, cut_start_ixz+1]
    cell.gxi_hw = mesh[j_start_ix+1, cut_start_iy+1+1:cut_end_iy-1, cut_end_ixz-1]
    cell.gxi_sw = mesh[j_start_ix+1, cut_start_iy+1, cut_start_ixz+1+1:cut_end_ixz-1]
    cell.gxi_nw = mesh[j_start_ix+1, cut_end_iy-1, cut_start_ixz+1+1:cut_end_ixz-1]

    cell.gxi_le = mesh[j_end_ix-1, cut_start_iy+1+1:cut_end_iy-1, cut_start_ixz+1]
    cell.gxi_he = mesh[j_end_ix-1, cut_start_iy+1+1:cut_end_iy-1, cut_end_ixz-1]
    cell.gxi_se = mesh[j_end_ix-1, cut_start_iy+1, cut_start_ixz+1+1:cut_end_ixz-1]
    cell.gxi_ne = mesh[j_end_ix-1, cut_end_iy-1, cut_start_ixz+1+1:cut_end_ixz-1]

    cell.gyi_ls = mesh[cut_start_ix+1+1:cut_end_ix-1, j_start_iy+1, cut_start_iyz+1]
    cell.gyi_hs = mesh[cut_start_ix+1+1:cut_end_ix-1, j_start_iy+1, cut_end_iyz-1]
    cell.gyi_sw = mesh[cut_start_ix+1, j_start_iy+1, cut_start_ixz+1+1:cut_end_iyz-1]
    cell.gyi_se = mesh[cut_end_ix-1, j_start_iy+1, cut_start_ixz+1+1:cut_end_iyz-1]

    cell.gyi_ln = mesh[cut_start_ix+1+1:cut_end_ix-1, j_end_iy-1, cut_start_iyz+1]
    cell.gyi_hn = mesh[cut_start_ix+1+1:cut_end_ix-1, j_end_iy-1, cut_end_iyz-1]
    cell.gyi_nw = mesh[cut_start_ix+1, j_end_iy-1, cut_start_ixz+1+1:cut_end_iyz-1]
    cell.gyi_ne = mesh[cut_end_ix-1, j_end_iy-1, cut_start_ixz+1+1:cut_end_iyz-1]


    # Membrane planes between cells
    cell.gyi_s = mesh[cut_start_ix+1+1:cut_end_ix-1, j_start_iy+1, cut_start_iyz+1+1:cut_end_iyz-1]
    cell.gyi_n = mesh[cut_start_ix+1+1:cut_end_ix-1, j_end_iy-1, cut_start_iyz+1+1:cut_end_iyz-1]

    cell.gxi_w = mesh[j_start_ix+1, cut_start_iy+1+1:cut_end_iy-1, cut_start_iyz+1+1:cut_end_iyz-1]
    cell.gxi_e = mesh[j_end_ix-1, cut_start_iy+1+1:cut_end_iy-1, cut_start_iyz+1+1:cut_end_iyz-1]


    # Set up indices for all inner points
    inner_cell = mesh[cell_start_ix:cell_end_ix+1, cell_start_iy:cell_end_iy+1, cell_start_iz:cell_end_iz+1]
    west_j = mesh[j_start_ix:cell_start_ix, cut_start_iy:cut_end_iy+1, cut_start_ixz:cut_end_ixz+1]
    east_j = mesh[cell_end_ix+1:j_end_ix+1, cut_start_iy:cut_end_iy+1, cut_start_ixz:cut_end_ixz+1]
    south_j = mesh[cut_start_ix:cut_end_ix+1, j_start_iy:cell_start_iy, cut_start_iyz:cut_end_iyz+1]
    north_j = mesh[cut_start_ix:cut_end_ix+1, cell_end_iy+1:j_end_iy+1, cut_start_iyz:cut_end_iyz+1]

    cell.i_all = np.concatenate([
        inner_cell,
        west_j,
        east_j,
        south_j,
        north_j,
    ])
    cell.i_all.sort()

    vec = np.zeros(G.N, dtype=np.int64)
    vec[cell.i_all] = 1
    not_strictly_inner = [
        cell.m_lsw, cell.m_lse, cell.m_lnw, cell.m_lne, cell.m_hsw, cell.m_hse,
        cell.m_hnw, cell.m_hne, cell.m_hw, cell.m_he, cell.m_hs, cell.m_hn,
        cell.m_lw, cell.m_le, cell.m_ls, cell.m_ln, cell.m_ne, cell.m_sw,
        cell.m_se, cell.m_nw, cell.m_w, cell.m_e, cell.m_s, cell.m_n, cell.m_h,
        cell.m_l, cell.i_lsw, cell.i_lse, cell.i_lnw, cell.i_lne, cell.i_hsw,
        cell.i_hse, cell.i_hnw, cell.i_hne, cell.i_hw, cell.i_he, cell.i_hs,
        cell.i_hn, cell.i_lw, cell.i_le, cell.i_ls, cell.i_ln, cell.i_ne,
        cell.i_sw, cell.i_se, cell.i_nw, cell.i_w, cell.i_e, cell.i_s, cell.i_n,
        cell.i_h, cell.i_l, cell.gx_w, cell.gx_e, cell.gx_hw, cell.gx_lw, cell.gx_sw,
        cell.gx_nw, cell.gx_he, cell.gx_le, cell.gx_se, cell.gx_ne, cell.gy_hs,
        cell.gy_ls, cell.gy_sw, cell.gy_se, cell.gy_hn, cell.gy_ln, cell.gy_nw,
        cell.gy_ne, cell.gy_n, cell.gy_s, cell.gx_lsw, cell.gx_lse, cell.gx_lnw,
        cell.gx_lne, cell.gx_hsw, cell.gx_hse, cell.gx_hnw, cell.gx_hne, cell.gy_lsw,
        cell.gy_lse, cell.gy_lnw, cell.gy_lne, cell.gy_hsw, cell.gy_hse, cell.gy_hnw,
        cell.gy_hne, cell.gxi_w, cell.gxi_e, cell.gxi_hw, cell.gxi_lw, cell.gxi_sw,
        cell.gxi_nw, cell.gxi_he, cell.gxi_le, cell.gxi_se, cell.gxi_ne, cell.gyi_hs,
        cell.gyi_ls, cell.gyi_sw, cell.gyi_se, cell.gyi_hn, cell.gyi_ln, cell.gyi_nw,
        cell.gyi_ne, cell.gyi_n, cell.gyi_s, cell.gxi_lsw, cell.gxi_lse, cell.gxi_lnw,
        cell.gxi_lne, cell.gxi_hsw, cell.gxi_hse, cell.gxi_hnw, cell.gxi_hne,
        cell.gyi_lsw, cell.gyi_lse, cell.gyi_lnw, cell.gyi_lne, cell.gyi_hsw,
        cell.gyi_hse, cell.gyi_hnw, cell.gyi_hne
    ]
    for index_set in not_strictly_inner:
        vec[index_set] = 0
    cell.i = np.flatnonzero(vec)


    return cell
