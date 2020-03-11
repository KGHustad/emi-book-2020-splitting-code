import numpy as np

from set_up_cell import set_up_cell, Cell
from set_up_cell_matrix import set_up_cell_matrix

WEST_MASK = 2**0
EAST_MASK = 2**1
SOUTH_MASK = 2**2
NORTH_MASK = 2**3
NUM_CELL_ARCHETYPES = 2**4

def ind_as_array(val):
    if type(val) is np.ndarray:
        return val
    elif type(val) is int:
        return np.array([val])
    else:
        raise ValueError("Unsupported type {0}".format(type(val)))

def set_up_archetype(G, mesh, i):
    neigh_in_west = bool(i & WEST_MASK)
    neigh_in_east = bool(i & EAST_MASK)
    neigh_in_south = bool(i & SOUTH_MASK)
    neigh_in_north = bool(i & NORTH_MASK)

    cell = set_up_cell(G, mesh)

    archetype = Cell()

    archetype.j_start_ix = cell.j_start_ix
    archetype.j_start_iy = cell.j_start_iy
    archetype.j_end_ix = cell.j_end_ix
    archetype.j_end_iy = cell.j_end_iy

    # Save indices for the cell
    archetype.m_lsw = ind_as_array(ind_as_array(cell.m_lsw))
    archetype.m_lse = ind_as_array(ind_as_array(cell.m_lse))
    archetype.m_lnw = ind_as_array(ind_as_array(cell.m_lnw))
    archetype.m_lne = ind_as_array(ind_as_array(cell.m_lne))
    archetype.m_hsw = ind_as_array(ind_as_array(cell.m_hsw))
    archetype.m_hse = ind_as_array(ind_as_array(cell.m_hse))
    archetype.m_hnw = ind_as_array(ind_as_array(cell.m_hnw))
    archetype.m_hne = ind_as_array(ind_as_array(cell.m_hne))

    archetype.m_hw = ind_as_array(ind_as_array(cell.m_hw))
    archetype.m_he = ind_as_array(ind_as_array(cell.m_he))
    archetype.m_hs = ind_as_array(ind_as_array(cell.m_hs))
    archetype.m_hn = ind_as_array(ind_as_array(cell.m_hn))
    archetype.m_lw = ind_as_array(ind_as_array(cell.m_lw))
    archetype.m_le = ind_as_array(ind_as_array(cell.m_le))
    archetype.m_ls = ind_as_array(ind_as_array(cell.m_ls))
    archetype.m_ln = ind_as_array(ind_as_array(cell.m_ln))
    archetype.m_ne = ind_as_array(ind_as_array(cell.m_ne))
    archetype.m_sw = ind_as_array(ind_as_array(cell.m_sw))
    archetype.m_se = ind_as_array(ind_as_array(cell.m_se))
    archetype.m_nw = ind_as_array(ind_as_array(cell.m_nw))

    archetype.m_w = ind_as_array(ind_as_array(cell.m_w))
    archetype.m_e = ind_as_array(ind_as_array(cell.m_e))
    archetype.m_s = ind_as_array(ind_as_array(cell.m_s))
    archetype.m_n = ind_as_array(ind_as_array(cell.m_n))
    archetype.m_h = ind_as_array(ind_as_array(cell.m_h))
    archetype.m_l = ind_as_array(ind_as_array(cell.m_l))

    archetype.i_lsw = ind_as_array(ind_as_array(cell.i_lsw))
    archetype.i_lse = ind_as_array(ind_as_array(cell.i_lse))
    archetype.i_lnw = ind_as_array(ind_as_array(cell.i_lnw))
    archetype.i_lne = ind_as_array(ind_as_array(cell.i_lne))
    archetype.i_hsw = ind_as_array(ind_as_array(cell.i_hsw))
    archetype.i_hse = ind_as_array(ind_as_array(cell.i_hse))
    archetype.i_hnw = ind_as_array(ind_as_array(cell.i_hnw))
    archetype.i_hne = ind_as_array(ind_as_array(cell.i_hne))

    archetype.i_hw = ind_as_array(ind_as_array(cell.i_hw))
    archetype.i_he = ind_as_array(ind_as_array(cell.i_he))
    archetype.i_hs = ind_as_array(ind_as_array(cell.i_hs))
    archetype.i_hn = ind_as_array(ind_as_array(cell.i_hn))
    archetype.i_lw = ind_as_array(ind_as_array(cell.i_lw))
    archetype.i_le = ind_as_array(ind_as_array(cell.i_le))
    archetype.i_ls = ind_as_array(ind_as_array(cell.i_ls))
    archetype.i_ln = ind_as_array(ind_as_array(cell.i_ln))
    archetype.i_ne = ind_as_array(ind_as_array(cell.i_ne))
    archetype.i_sw = ind_as_array(ind_as_array(cell.i_sw))
    archetype.i_se = ind_as_array(ind_as_array(cell.i_se))
    archetype.i_nw = ind_as_array(ind_as_array(cell.i_nw))

    archetype.i_w = ind_as_array(ind_as_array(cell.i_w))
    archetype.i_e = ind_as_array(ind_as_array(cell.i_e))
    archetype.i_s = ind_as_array(ind_as_array(cell.i_s))
    archetype.i_n = ind_as_array(ind_as_array(cell.i_n))
    archetype.i_h = ind_as_array(ind_as_array(cell.i_h))
    archetype.i_l = ind_as_array(ind_as_array(cell.i_l))
    archetype.i = ind_as_array(ind_as_array(cell.i))
    archetype.i_all = ind_as_array(ind_as_array(cell.i_all))

    archetype.gx_w = ind_as_array(ind_as_array(cell.gx_w))
    archetype.gx_e = ind_as_array(ind_as_array(cell.gx_e))
    archetype.gx_hw = ind_as_array(ind_as_array(cell.gx_hw))
    archetype.gx_lw = ind_as_array(ind_as_array(cell.gx_lw))
    archetype.gx_sw = ind_as_array(ind_as_array(cell.gx_sw))
    archetype.gx_nw = ind_as_array(ind_as_array(cell.gx_nw))
    archetype.gx_he = ind_as_array(ind_as_array(cell.gx_he))
    archetype.gx_le = ind_as_array(ind_as_array(cell.gx_le))
    archetype.gx_se = ind_as_array(ind_as_array(cell.gx_se))
    archetype.gx_ne = ind_as_array(ind_as_array(cell.gx_ne))

    archetype.gy_hs = ind_as_array(ind_as_array(cell.gy_hs))
    archetype.gy_ls = ind_as_array(ind_as_array(cell.gy_ls))
    archetype.gy_sw = ind_as_array(ind_as_array(cell.gy_sw))
    archetype.gy_se = ind_as_array(ind_as_array(cell.gy_se))
    archetype.gy_hn = ind_as_array(ind_as_array(cell.gy_hn))
    archetype.gy_ln = ind_as_array(ind_as_array(cell.gy_ln))
    archetype.gy_nw = ind_as_array(ind_as_array(cell.gy_nw))
    archetype.gy_ne = ind_as_array(ind_as_array(cell.gy_ne))
    archetype.gy_n = ind_as_array(ind_as_array(cell.gy_n))
    archetype.gy_s = ind_as_array(ind_as_array(cell.gy_s))

    archetype.gx_lsw = ind_as_array(ind_as_array(cell.gx_lsw))
    archetype.gx_lse = ind_as_array(ind_as_array(cell.gx_lse))
    archetype.gx_lnw = ind_as_array(ind_as_array(cell.gx_lne))
    archetype.gx_lne = ind_as_array(ind_as_array(cell.gx_lne))
    archetype.gx_hsw = ind_as_array(ind_as_array(cell.gx_hsw))
    archetype.gx_hse = ind_as_array(ind_as_array(cell.gx_hse))
    archetype.gx_hnw = ind_as_array(ind_as_array(cell.gx_hnw))
    archetype.gx_hne = ind_as_array(ind_as_array(cell.gx_hne))

    archetype.gy_lsw = ind_as_array(ind_as_array(cell.gy_lsw))
    archetype.gy_lse = ind_as_array(ind_as_array(cell.gy_lse))
    archetype.gy_lnw = ind_as_array(ind_as_array(cell.gy_lnw))
    archetype.gy_lne = ind_as_array(ind_as_array(cell.gy_lne))
    archetype.gy_hsw = ind_as_array(ind_as_array(cell.gy_hsw))
    archetype.gy_hse = ind_as_array(ind_as_array(cell.gy_hse))
    archetype.gy_hnw = ind_as_array(ind_as_array(cell.gy_hnw))
    archetype.gy_hne = ind_as_array(ind_as_array(cell.gy_hne))

    archetype.gxi_w = ind_as_array(ind_as_array(cell.gxi_w))
    archetype.gxi_e = ind_as_array(ind_as_array(cell.gxi_e))
    archetype.gxi_hw = ind_as_array(ind_as_array(cell.gxi_hw))
    archetype.gxi_lw = ind_as_array(ind_as_array(cell.gxi_lw))
    archetype.gxi_sw = ind_as_array(ind_as_array(cell.gxi_sw))
    archetype.gxi_nw = ind_as_array(ind_as_array(cell.gxi_nw))
    archetype.gxi_he = ind_as_array(ind_as_array(cell.gxi_he))
    archetype.gxi_le = ind_as_array(ind_as_array(cell.gxi_le))
    archetype.gxi_se = ind_as_array(ind_as_array(cell.gxi_se))
    archetype.gxi_ne = ind_as_array(ind_as_array(cell.gxi_ne))

    archetype.gyi_hs = ind_as_array(ind_as_array(cell.gyi_hs))
    archetype.gyi_ls = ind_as_array(ind_as_array(cell.gyi_ls))
    archetype.gyi_sw = ind_as_array(ind_as_array(cell.gyi_sw))
    archetype.gyi_se = ind_as_array(ind_as_array(cell.gyi_se))
    archetype.gyi_hn = ind_as_array(ind_as_array(cell.gyi_hn))
    archetype.gyi_ln = ind_as_array(ind_as_array(cell.gyi_ln))
    archetype.gyi_nw = ind_as_array(ind_as_array(cell.gyi_nw))
    archetype.gyi_ne = ind_as_array(ind_as_array(cell.gyi_ne))
    archetype.gyi_n = ind_as_array(ind_as_array(cell.gyi_n))
    archetype.gyi_s = ind_as_array(ind_as_array(cell.gyi_s))

    archetype.gxi_lsw = ind_as_array(ind_as_array(cell.gxi_lsw))
    archetype.gxi_lse = ind_as_array(ind_as_array(cell.gxi_lse))
    archetype.gxi_lnw = ind_as_array(ind_as_array(cell.gxi_lnw))
    archetype.gxi_lne = ind_as_array(ind_as_array(cell.gxi_lne))
    archetype.gxi_hsw = ind_as_array(ind_as_array(cell.gxi_hsw))
    archetype.gxi_hse = ind_as_array(ind_as_array(cell.gxi_hse))
    archetype.gxi_hnw = ind_as_array(ind_as_array(cell.gxi_hnw))
    archetype.gxi_hne = ind_as_array(ind_as_array(cell.gxi_hne))

    archetype.gyi_lsw = ind_as_array(ind_as_array(cell.gyi_lsw))
    archetype.gyi_lse = ind_as_array(ind_as_array(cell.gyi_lse))
    archetype.gyi_lnw = ind_as_array(ind_as_array(cell.gyi_lnw))
    archetype.gyi_lne = ind_as_array(ind_as_array(cell.gyi_lne))
    archetype.gyi_hsw = ind_as_array(ind_as_array(cell.gyi_hsw))
    archetype.gyi_hse = ind_as_array(ind_as_array(cell.gyi_hse))
    archetype.gyi_hnw = ind_as_array(ind_as_array(cell.gyi_hnw))
    archetype.gyi_hne = ind_as_array(ind_as_array(cell.gyi_hne))

    empty_array = np.zeros(0, dtype=np.int64)
    if not neigh_in_west:
        archetype.gx_w = empty_array
        archetype.gx_sw = empty_array
        archetype.gx_nw = empty_array
        archetype.gx_lw = empty_array
        archetype.gx_hw = empty_array
        archetype.gx_hsw = empty_array
        archetype.gx_lnw = empty_array
        archetype.gx_lsw = empty_array
        archetype.gx_hnw = empty_array

        archetype.m_w = np.concatenate([archetype.m_w, ind_as_array(cell.gx_w)])
        archetype.m_sw = np.concatenate([archetype.m_sw, ind_as_array(cell.gx_sw)])
        archetype.m_nw = np.concatenate([archetype.m_nw, ind_as_array(cell.gx_nw)])
        archetype.m_lw = np.concatenate([archetype.m_lw, ind_as_array(cell.gx_lw)])
        archetype.m_hw = np.concatenate([archetype.m_hw, ind_as_array(cell.gx_hw)])
        archetype.m_hsw = np.concatenate([archetype.m_hsw, ind_as_array(cell.gx_hsw)])
        archetype.m_lnw = np.concatenate([archetype.m_lnw, ind_as_array(cell.gx_lnw)])
        archetype.m_lsw = np.concatenate([archetype.m_lsw, ind_as_array(cell.gx_lsw)])
        archetype.m_hnw = np.concatenate([archetype.m_hnw, ind_as_array(cell.gx_hnw)])

        archetype.gxi_w = empty_array
        archetype.gxi_sw = empty_array
        archetype.gxi_nw = empty_array
        archetype.gxi_lw = empty_array
        archetype.gxi_hw = empty_array
        archetype.gxi_hsw = empty_array
        archetype.gxi_lnw = empty_array
        archetype.gxi_lsw = empty_array
        archetype.gxi_hnw = empty_array

        archetype.i_w = np.concatenate([archetype.i_w, ind_as_array(cell.gxi_w)])
        archetype.i_sw = np.concatenate([archetype.i_sw, ind_as_array(cell.gxi_sw)])
        archetype.i_nw = np.concatenate([archetype.i_nw, ind_as_array(cell.gxi_nw)])
        archetype.i_lw = np.concatenate([archetype.i_lw, ind_as_array(cell.gxi_lw)])
        archetype.i_hw = np.concatenate([archetype.i_hw, ind_as_array(cell.gxi_hw)])
        archetype.i_hsw = np.concatenate([archetype.i_hsw, ind_as_array(cell.gxi_hsw)])
        archetype.i_lnw = np.concatenate([archetype.i_lnw, ind_as_array(cell.gxi_lnw)])
        archetype.i_lsw = np.concatenate([archetype.i_lsw, ind_as_array(cell.gxi_lsw)])
        archetype.i_hnw = np.concatenate([archetype.i_hnw, ind_as_array(cell.gxi_hnw)])
    if not neigh_in_east:
        archetype.gx_e = empty_array
        archetype.gx_se = empty_array
        archetype.gx_ne = empty_array
        archetype.gx_le = empty_array
        archetype.gx_he = empty_array
        archetype.gx_hse = empty_array
        archetype.gx_lne = empty_array
        archetype.gx_lse = empty_array
        archetype.gx_hne = empty_array

        archetype.m_e = np.concatenate([archetype.m_e, ind_as_array(cell.gx_e)])
        archetype.m_se = np.concatenate([archetype.m_se, ind_as_array(cell.gx_se)])
        archetype.m_ne = np.concatenate([archetype.m_ne, ind_as_array(cell.gx_ne)])
        archetype.m_le = np.concatenate([archetype.m_le, ind_as_array(cell.gx_le)])
        archetype.m_he = np.concatenate([archetype.m_he, ind_as_array(cell.gx_he)])
        archetype.m_hse = np.concatenate([archetype.m_hse, ind_as_array(cell.gx_hse)])
        archetype.m_lne = np.concatenate([archetype.m_lne, ind_as_array(cell.gx_lne)])
        archetype.m_lse = np.concatenate([archetype.m_lse, ind_as_array(cell.gx_lse)])
        archetype.m_hne = np.concatenate([archetype.m_hne, ind_as_array(cell.gx_hne)])

        archetype.gxi_e = empty_array
        archetype.gxi_se = empty_array
        archetype.gxi_ne = empty_array
        archetype.gxi_le = empty_array
        archetype.gxi_he = empty_array
        archetype.gxi_hse = empty_array
        archetype.gxi_lne = empty_array
        archetype.gxi_lse = empty_array
        archetype.gxi_hne = empty_array

        archetype.i_e = np.concatenate([archetype.i_e, ind_as_array(cell.gxi_e)])
        archetype.i_se = np.concatenate([archetype.i_se, ind_as_array(cell.gxi_se)])
        archetype.i_ne = np.concatenate([archetype.i_ne, ind_as_array(cell.gxi_ne)])
        archetype.i_le = np.concatenate([archetype.i_le, ind_as_array(cell.gxi_le)])
        archetype.i_he = np.concatenate([archetype.i_he, ind_as_array(cell.gxi_he)])
        archetype.i_hse = np.concatenate([archetype.i_hse, ind_as_array(cell.gxi_hse)])
        archetype.i_lne = np.concatenate([archetype.i_lne, ind_as_array(cell.gxi_lne)])
        archetype.i_lse = np.concatenate([archetype.i_lse, ind_as_array(cell.gxi_lse)])
        archetype.i_hne = np.concatenate([archetype.i_hne, ind_as_array(cell.gxi_hne)])
    if not neigh_in_south:
        archetype.gy_s = empty_array
        archetype.gy_sw = empty_array
        archetype.gy_se = empty_array
        archetype.gy_ls = empty_array
        archetype.gy_hs = empty_array
        archetype.gy_hsw = empty_array
        archetype.gy_lse = empty_array
        archetype.gy_lsw = empty_array
        archetype.gy_hse = empty_array

        archetype.m_s = np.concatenate([archetype.m_s, ind_as_array(cell.gy_s)])
        archetype.m_sw = np.concatenate([archetype.m_sw, ind_as_array(cell.gy_sw)])
        archetype.m_se = np.concatenate([archetype.m_se, ind_as_array(cell.gy_se)])
        archetype.m_ls = np.concatenate([archetype.m_ls, ind_as_array(cell.gy_ls)])
        archetype.m_hs = np.concatenate([archetype.m_hs, ind_as_array(cell.gy_hs)])
        archetype.m_hsw = np.concatenate([archetype.m_hsw, ind_as_array(cell.gy_hsw)])
        archetype.m_lse = np.concatenate([archetype.m_lse, ind_as_array(cell.gy_lse)])
        archetype.m_lsw = np.concatenate([archetype.m_lsw, ind_as_array(cell.gy_lsw)])
        archetype.m_hse = np.concatenate([archetype.m_hse, ind_as_array(cell.gy_hse)])

        archetype.gyi_s = empty_array
        archetype.gyi_sw = empty_array
        archetype.gyi_se = empty_array
        archetype.gyi_ls = empty_array
        archetype.gyi_hs = empty_array
        archetype.gyi_hsw = empty_array
        archetype.gyi_lse = empty_array
        archetype.gyi_lsw = empty_array
        archetype.gyi_hse = empty_array

        archetype.i_s = np.concatenate([archetype.i_s, ind_as_array(cell.gyi_s)])
        archetype.i_sw = np.concatenate([archetype.i_sw, ind_as_array(cell.gyi_sw)])
        archetype.i_se = np.concatenate([archetype.i_se, ind_as_array(cell.gyi_se)])
        archetype.i_ls = np.concatenate([archetype.i_ls, ind_as_array(cell.gyi_ls)])
        archetype.i_hs = np.concatenate([archetype.i_hs, ind_as_array(cell.gyi_hs)])
        archetype.i_hsw = np.concatenate([archetype.i_hsw, ind_as_array(cell.gyi_hsw)])
        archetype.i_lse = np.concatenate([archetype.i_lse, ind_as_array(cell.gyi_lse)])
        archetype.i_lsw = np.concatenate([archetype.i_lsw, ind_as_array(cell.gyi_lsw)])
        archetype.i_hse = np.concatenate([archetype.i_hse, ind_as_array(cell.gyi_hse)])
    if not neigh_in_north:
        archetype.gy_n = empty_array
        archetype.gy_nw = empty_array
        archetype.gy_ne = empty_array
        archetype.gy_ln = empty_array
        archetype.gy_hn = empty_array
        archetype.gy_hnw = empty_array
        archetype.gy_lne = empty_array
        archetype.gy_lnw = empty_array
        archetype.gy_hne = empty_array

        archetype.m_n = np.concatenate([archetype.m_n, ind_as_array(cell.gy_n)])
        archetype.m_nw = np.concatenate([archetype.m_nw, ind_as_array(cell.gy_nw)])
        archetype.m_ne = np.concatenate([archetype.m_ne, ind_as_array(cell.gy_ne)])
        archetype.m_ln = np.concatenate([archetype.m_ln, ind_as_array(cell.gy_ln)])
        archetype.m_hn = np.concatenate([archetype.m_hn, ind_as_array(cell.gy_hn)])
        archetype.m_hnw = np.concatenate([archetype.m_hnw, ind_as_array(cell.gy_hnw)])
        archetype.m_lne = np.concatenate([archetype.m_lne, ind_as_array(cell.gy_lne)])
        archetype.m_lnw = np.concatenate([archetype.m_lnw, ind_as_array(cell.gy_lnw)])
        archetype.m_hne = np.concatenate([archetype.m_hne, ind_as_array(cell.gy_hne)])

        archetype.gyi_n = empty_array
        archetype.gyi_nw = empty_array
        archetype.gyi_ne = empty_array
        archetype.gyi_ln = empty_array
        archetype.gyi_hn = empty_array
        archetype.gyi_hnw = empty_array
        archetype.gyi_lne = empty_array
        archetype.gyi_lnw = empty_array
        archetype.gyi_hne = empty_array

        archetype.i_n = np.concatenate([archetype.i_n, ind_as_array(cell.gyi_n)])
        archetype.i_nw = np.concatenate([archetype.i_nw, ind_as_array(cell.gyi_nw)])
        archetype.i_ne = np.concatenate([archetype.i_ne, ind_as_array(cell.gyi_ne)])
        archetype.i_ln = np.concatenate([archetype.i_ln, ind_as_array(cell.gyi_ln)])
        archetype.i_hn = np.concatenate([archetype.i_hn, ind_as_array(cell.gyi_hn)])
        archetype.i_hnw = np.concatenate([archetype.i_hnw, ind_as_array(cell.gyi_hnw)])
        archetype.i_lne = np.concatenate([archetype.i_lne, ind_as_array(cell.gyi_lne)])
        archetype.i_lnw = np.concatenate([archetype.i_lnw, ind_as_array(cell.gyi_lnw)])
        archetype.i_hne = np.concatenate([archetype.i_hne, ind_as_array(cell.gyi_hne)])

    # Set up membrane potential
    archetype.v = np.concatenate([
        archetype.m_lsw,
        archetype.m_lse,
        archetype.m_lnw,
        archetype.m_lne,
        archetype.m_hsw,
        archetype.m_hse,
        archetype.m_hnw,
        archetype.m_hne,
        archetype.m_hw,
        archetype.m_he,
        archetype.m_hs,
        archetype.m_hn,
        archetype.m_lw,
        archetype.m_le,
        archetype.m_ls,
        archetype.m_ln,
        archetype.m_ne,
        archetype.m_sw,
        archetype.m_se,
        archetype.m_nw,
        archetype.m_w,
        archetype.m_e,
        archetype.m_s,
        archetype.m_n,
        archetype.m_h,
        archetype.m_l
    ])
    archetype.v.sort()

    # Set up all indices for the cell
    archetype.c_all = np.unique(np.concatenate([
        archetype.v,
        archetype.i,
        archetype.i_lsw,
        archetype.i_lse,
        archetype.i_lnw,
        archetype.i_lne,
        archetype.i_hsw,
        archetype.i_hse,
        archetype.i_hnw,
        archetype.i_hne,
        archetype.i_hw,
        archetype.i_he,
        archetype.i_hs,
        archetype.i_hn,
        archetype.i_lw,
        archetype.i_le,
        archetype.i_ls,
        archetype.i_ln,
        archetype.i_ne,
        archetype.i_sw,
        archetype.i_se,
        archetype.i_nw,
        archetype.i_w,
        archetype.i_e,
        archetype.i_s,
        archetype.i_n,
        archetype.i_h,
        archetype.i_l,
        archetype.gx_w,
        archetype.gx_e,
        archetype.gxi_w,
        archetype.gxi_sw,
        archetype.gxi_nw,
        archetype.gxi_lw,
        archetype.gxi_hw,
        archetype.gxi_hsw,
        archetype.gxi_lnw,
        archetype.gxi_lsw,
        archetype.gxi_hnw,
        archetype.gxi_e,
        archetype.gxi_se,
        archetype.gxi_ne,
        archetype.gxi_le,
        archetype.gxi_he,
        archetype.gxi_hse,
        archetype.gxi_lne,
        archetype.gxi_lse,
        archetype.gxi_hne,
        archetype.gyi_n,
        archetype.gyi_nw,
        archetype.gyi_ne,
        archetype.gyi_ln,
        archetype.gyi_hn,
        archetype.gyi_hnw,
        archetype.gyi_lne,
        archetype.gyi_lnw,
        archetype.gyi_hne,
        archetype.gyi_s,
        archetype.gyi_sw,
        archetype.gyi_se,
        archetype.gyi_ls,
        archetype.gyi_hs,
        archetype.gyi_hsw,
        archetype.gyi_lse,
        archetype.gyi_lsw,
        archetype.gyi_hse,
        archetype.gy_n,
        archetype.gy_s
    ]))

    # Set up all indices for the cell (except the membrane)
    archetype.i = np.unique(np.concatenate([
        archetype.i,
        archetype.i_lsw,
        archetype.i_lse,
        archetype.i_lnw,
        archetype.i_lne,
        archetype.i_hsw,
        archetype.i_hse,
        archetype.i_hnw,
        archetype.i_hne,
        archetype.i_hw,
        archetype.i_he,
        archetype.i_hs,
        archetype.i_hn,
        archetype.i_lw,
        archetype.i_le,
        archetype.i_ls,
        archetype.i_ln,
        archetype.i_ne,
        archetype.i_sw,
        archetype.i_se,
        archetype.i_nw,
        archetype.i_w,
        archetype.i_e,
        archetype.i_s,
        archetype.i_n,
        archetype.i_h,
        archetype.i_l,
        archetype.gxi_w,
        archetype.gxi_sw,
        archetype.gxi_nw,
        archetype.gxi_lw,
        archetype.gxi_hw,
        archetype.gxi_hsw,
        archetype.gxi_lnw,
        archetype.gxi_lsw,
        archetype.gxi_hnw,
        archetype.gxi_e,
        archetype.gxi_se,
        archetype.gxi_ne,
        archetype.gxi_le,
        archetype.gxi_he,
        archetype.gxi_hse,
        archetype.gxi_lne,
        archetype.gxi_lse,
        archetype.gxi_hne,
        archetype.gyi_n,
        archetype.gyi_nw,
        archetype.gyi_ne,
        archetype.gyi_ln,
        archetype.gyi_hn,
        archetype.gyi_hnw,
        archetype.gyi_lne,
        archetype.gyi_lnw,
        archetype.gyi_hne,
        archetype.gyi_s,
        archetype.gyi_sw,
        archetype.gyi_se,
        archetype.gyi_ls,
        archetype.gyi_hs,
        archetype.gyi_hsw,
        archetype.gyi_lse,
        archetype.gyi_lsw,
        archetype.gyi_hse,
    ]))

    # Set up local indices for the membrane potential
    vec = np.zeros(G.N, dtype=np.int64)
    vec[archetype.v] = 1
    vec = vec[archetype.c_all]
    archetype.v_of_c = np.flatnonzero(vec)

    # Set up local indices for the intercalated discs
    # West junction
    vec = np.zeros(G.N, dtype=np.int64)
    vec[archetype.gx_w] = 1
    vec = vec[archetype.c_all]
    archetype.gxw_of_c = np.flatnonzero(vec)

    # East junction
    vec = np.zeros(G.N, dtype=np.int64)
    vec[archetype.gx_e] = 1
    vec = vec[archetype.c_all]
    archetype.gxe_of_c = np.flatnonzero(vec)

    # South junction
    vec = np.zeros(G.N, dtype=np.int64)
    vec[archetype.gy_s] = 1
    vec = vec[archetype.c_all]
    archetype.gys_of_c = np.flatnonzero(vec)

    # North junction
    vec = np.zeros(G.N, dtype=np.int64)
    vec[archetype.gy_n] = 1
    vec = vec[archetype.c_all]
    archetype.gyn_of_c = np.flatnonzero(vec)

    # Set up intracellular system
    archetype.A = set_up_cell_matrix(G, mesh, archetype)

    return archetype

def combine_cells(G, mesh):
    # there are 2^4 cell archetypes: for each of the 4 sides, they can have a
    # neighbour or not
    archetypes = []
    for i in range(NUM_CELL_ARCHETYPES):
        archetype = set_up_archetype(G, mesh, i)
        archetypes.append(archetype)

    # create cells from archetypes
    empty_array = np.zeros(0, dtype=np.int64)
    cells = []
    for j in range(G.num_cells_y):
        for i in range(G.num_cells_x):
            n = j * G.num_cells_x + i

            c_idx = 0

            if i != 0:
                c_idx |= WEST_MASK

            if i != G.num_cells_x - 1:
                c_idx |= EAST_MASK

            if j != 0:
                c_idx |= SOUTH_MASK

            if j != G.num_cells_y - 1:
                c_idx |= NORTH_MASK

            cell = archetypes[c_idx]
            cells.append(Cell())

            offset_x = int(round((G.cell_start_x[i] - G.cell_start_x[0]) / G.dx))
            offset_y = int(round((G.cell_start_y[j] - G.cell_start_y[0]) / G.dy))

            #offset = offset_y * G.Nx + offset_x;
            offset = mesh[offset_x, offset_y, 0]

            cells[n].archetype = c_idx
            #print('Cell %d at grid pos (%d, %d) has archetype %d. Offsetting with (%d, %d)=%d\n' % (n, i, j, c_idx, offset_x, offset_y, offset))

            # copy intracellular matrix
            cells[n].A = cell.A

            cells[n].j_start_ix = cell.j_start_ix + offset_x
            cells[n].j_start_iy = cell.j_start_iy + offset_y
            cells[n].j_end_ix = cell.j_end_ix + offset_x
            cells[n].j_end_iy = cell.j_end_iy + offset_y

            # Displace all index sets by offset
            cells[n].m_lsw = cell.m_lsw + offset
            cells[n].m_lse = cell.m_lse + offset
            cells[n].m_lnw = cell.m_lnw + offset
            cells[n].m_lne = cell.m_lne + offset
            cells[n].m_hsw = cell.m_hsw + offset
            cells[n].m_hse = cell.m_hse + offset
            cells[n].m_hnw = cell.m_hnw + offset
            cells[n].m_hne = cell.m_hne + offset
            cells[n].m_hw = cell.m_hw + offset
            cells[n].m_he = cell.m_he + offset
            cells[n].m_hs = cell.m_hs + offset
            cells[n].m_hn = cell.m_hn + offset
            cells[n].m_lw = cell.m_lw + offset
            cells[n].m_le = cell.m_le + offset
            cells[n].m_ls = cell.m_ls + offset
            cells[n].m_ln = cell.m_ln + offset
            cells[n].m_ne = cell.m_ne + offset
            cells[n].m_sw = cell.m_sw + offset
            cells[n].m_se = cell.m_se + offset
            cells[n].m_nw = cell.m_nw + offset
            cells[n].m_w = cell.m_w + offset
            cells[n].m_e = cell.m_e + offset
            cells[n].m_s = cell.m_s + offset
            cells[n].m_n = cell.m_n + offset
            cells[n].m_h = cell.m_h + offset
            cells[n].m_l = cell.m_l + offset
            cells[n].i_lsw = cell.i_lsw + offset
            cells[n].i_lse = cell.i_lse + offset
            cells[n].i_lnw = cell.i_lnw + offset
            cells[n].i_lne = cell.i_lne + offset
            cells[n].i_hsw = cell.i_hsw + offset
            cells[n].i_hse = cell.i_hse + offset
            cells[n].i_hnw = cell.i_hnw + offset
            cells[n].i_hne = cell.i_hne + offset
            cells[n].i_hw = cell.i_hw + offset
            cells[n].i_he = cell.i_he + offset
            cells[n].i_hs = cell.i_hs + offset
            cells[n].i_hn = cell.i_hn + offset
            cells[n].i_lw = cell.i_lw + offset
            cells[n].i_le = cell.i_le + offset
            cells[n].i_ls = cell.i_ls + offset
            cells[n].i_ln = cell.i_ln + offset
            cells[n].i_ne = cell.i_ne + offset
            cells[n].i_sw = cell.i_sw + offset
            cells[n].i_se = cell.i_se + offset
            cells[n].i_nw = cell.i_nw + offset
            cells[n].i_w = cell.i_w + offset
            cells[n].i_e = cell.i_e + offset
            cells[n].i_s = cell.i_s + offset
            cells[n].i_n = cell.i_n + offset
            cells[n].i_h = cell.i_h + offset
            cells[n].i_l = cell.i_l + offset
            cells[n].i = cell.i + offset
            cells[n].i_all = cell.i_all + offset
            cells[n].gx_w = cell.gx_w + offset
            cells[n].gx_e = cell.gx_e + offset
            cells[n].gx_hw = cell.gx_hw + offset
            cells[n].gx_lw = cell.gx_lw + offset
            cells[n].gx_sw = cell.gx_sw + offset
            cells[n].gx_nw = cell.gx_nw + offset
            cells[n].gx_he = cell.gx_he + offset
            cells[n].gx_le = cell.gx_le + offset
            cells[n].gx_se = cell.gx_se + offset
            cells[n].gx_ne = cell.gx_ne + offset
            cells[n].gy_hs = cell.gy_hs + offset
            cells[n].gy_ls = cell.gy_ls + offset
            cells[n].gy_sw = cell.gy_sw + offset
            cells[n].gy_se = cell.gy_se + offset
            cells[n].gy_hn = cell.gy_hn + offset
            cells[n].gy_ln = cell.gy_ln + offset
            cells[n].gy_nw = cell.gy_nw + offset
            cells[n].gy_ne = cell.gy_ne + offset
            cells[n].gy_n = cell.gy_n + offset
            cells[n].gy_s = cell.gy_s + offset
            cells[n].gx_lsw = cell.gx_lsw + offset
            cells[n].gx_lse = cell.gx_lse + offset
            cells[n].gx_lnw = cell.gx_lnw + offset
            cells[n].gx_lne = cell.gx_lne + offset
            cells[n].gx_hsw = cell.gx_hsw + offset
            cells[n].gx_hse = cell.gx_hse + offset
            cells[n].gx_hnw = cell.gx_hnw + offset
            cells[n].gx_hne = cell.gx_hne + offset
            cells[n].gy_lsw = cell.gy_lsw + offset
            cells[n].gy_lse = cell.gy_lse + offset
            cells[n].gy_lnw = cell.gy_lnw + offset
            cells[n].gy_lne = cell.gy_lne + offset
            cells[n].gy_hsw = cell.gy_hsw + offset
            cells[n].gy_hse = cell.gy_hse + offset
            cells[n].gy_hnw = cell.gy_hnw + offset
            cells[n].gy_hne = cell.gy_hne + offset
            cells[n].gxi_w = cell.gxi_w + offset
            cells[n].gxi_e = cell.gxi_e + offset
            cells[n].gxi_hw = cell.gxi_hw + offset
            cells[n].gxi_lw = cell.gxi_lw + offset
            cells[n].gxi_sw = cell.gxi_sw + offset
            cells[n].gxi_nw = cell.gxi_nw + offset
            cells[n].gxi_he = cell.gxi_he + offset
            cells[n].gxi_le = cell.gxi_le + offset
            cells[n].gxi_se = cell.gxi_se + offset
            cells[n].gxi_ne = cell.gxi_ne + offset
            cells[n].gyi_hs = cell.gyi_hs + offset
            cells[n].gyi_ls = cell.gyi_ls + offset
            cells[n].gyi_sw = cell.gyi_sw + offset
            cells[n].gyi_se = cell.gyi_se + offset
            cells[n].gyi_hn = cell.gyi_hn + offset
            cells[n].gyi_ln = cell.gyi_ln + offset
            cells[n].gyi_nw = cell.gyi_nw + offset
            cells[n].gyi_ne = cell.gyi_ne + offset
            cells[n].gyi_n = cell.gyi_n + offset
            cells[n].gyi_s = cell.gyi_s + offset
            cells[n].gxi_lsw = cell.gxi_lsw + offset
            cells[n].gxi_lse = cell.gxi_lse + offset
            cells[n].gxi_lnw = cell.gxi_lnw + offset
            cells[n].gxi_lne = cell.gxi_lne + offset
            cells[n].gxi_hsw = cell.gxi_hsw + offset
            cells[n].gxi_hse = cell.gxi_hse + offset
            cells[n].gxi_hnw = cell.gxi_hnw + offset
            cells[n].gxi_hne = cell.gxi_hne + offset
            cells[n].gyi_lsw = cell.gyi_lsw + offset
            cells[n].gyi_lse = cell.gyi_lse + offset
            cells[n].gyi_lnw = cell.gyi_lnw + offset
            cells[n].gyi_lne = cell.gyi_lne + offset
            cells[n].gyi_hsw = cell.gyi_hsw + offset
            cells[n].gyi_hse = cell.gyi_hse + offset
            cells[n].gyi_hnw = cell.gyi_hnw + offset
            cells[n].gyi_hne = cell.gyi_hne + offset
            cells[n].v = cell.v + offset
            cells[n].c_all = cell.c_all + offset

            # Don't displace cell-local indices
            cells[n].v_of_c = cell.v_of_c
            cells[n].gxe_of_c = cell.gxe_of_c
            cells[n].gyn_of_c = cell.gyn_of_c
            cells[n].gxw_of_c = cell.gxw_of_c
            cells[n].gys_of_c = cell.gys_of_c

            # Set up indices to stimulate
            if G.num_stim[n] == 0:
                cells[n].to_stim = empty_array
            else:
                #[x, y, z] = meshgrid(cell.j_start_ix+G.stim_start(n):(cell.j_start_ix+G.stim_start(n)+G.num_stim(n)), 1:G.Ny, 1:G.Nz);
                #to_stim = sort(reshape(sub2ind([G.Nx,G.Ny,G.Nz], x, y, z), size(x,1)*size(x,2)*size(x, 3), 1));
                stim_start_ix = cells[n].j_start_ix + G.stim_start[n]
                stim_end_ix = stim_start_ix + G.num_stim[n] + 1
                stim_start_iy = cells[n].j_start_iy
                stim_end_iy = cells[n].j_end_iy + 1
                #print(stim_start_ix, stim_end_ix, stim_start_iy, stim_end_iy)
                to_stim = mesh[stim_start_ix:stim_end_ix, stim_start_iy:stim_end_iy, :]
                #indices = zeros(1,G. N);
                #indices(to_stim) = indices(to_stim) + 1;
                #indices(cells(n).v) = indices(cells(n).v) + 1;
                #cells(n).to_stim = round(find((indices==2)));
                cells[n].to_stim = np.intersect1d(to_stim, cells[n].v);
                #disp(cells(n).to_stim)

    # concatenate index sets for each cell into index sets for the whole mesh
    mesh.m_lsw = np.concatenate([cell.m_lsw for cell in cells])
    mesh.m_lse = np.concatenate([cell.m_lse for cell in cells])
    mesh.m_lnw = np.concatenate([cell.m_lnw for cell in cells])
    mesh.m_lne = np.concatenate([cell.m_lne for cell in cells])
    mesh.m_hsw = np.concatenate([cell.m_hsw for cell in cells])
    mesh.m_hse = np.concatenate([cell.m_hse for cell in cells])
    mesh.m_hnw = np.concatenate([cell.m_hnw for cell in cells])
    mesh.m_hne = np.concatenate([cell.m_hne for cell in cells])

    mesh.m_hw = np.concatenate([cell.m_hw for cell in cells])
    mesh.m_he = np.concatenate([cell.m_he for cell in cells])
    mesh.m_hs = np.concatenate([cell.m_hs for cell in cells])
    mesh.m_hn = np.concatenate([cell.m_hn for cell in cells])
    mesh.m_lw = np.concatenate([cell.m_lw for cell in cells])
    mesh.m_le = np.concatenate([cell.m_le for cell in cells])
    mesh.m_ls = np.concatenate([cell.m_ls for cell in cells])
    mesh.m_ln = np.concatenate([cell.m_ln for cell in cells])
    mesh.m_ne = np.concatenate([cell.m_ne for cell in cells])
    mesh.m_sw = np.concatenate([cell.m_sw for cell in cells])
    mesh.m_se = np.concatenate([cell.m_se for cell in cells])
    mesh.m_nw = np.concatenate([cell.m_nw for cell in cells])

    mesh.m_w = np.concatenate([cell.m_w for cell in cells])
    mesh.m_e = np.concatenate([cell.m_e for cell in cells])
    mesh.m_s = np.concatenate([cell.m_s for cell in cells])
    mesh.m_n = np.concatenate([cell.m_n for cell in cells])
    mesh.m_h = np.concatenate([cell.m_h for cell in cells])
    mesh.m_l = np.concatenate([cell.m_l for cell in cells])

    mesh.i_lsw = np.concatenate([cell.i_lsw for cell in cells])
    mesh.i_lse = np.concatenate([cell.i_lse for cell in cells])
    mesh.i_lnw = np.concatenate([cell.i_lnw for cell in cells])
    mesh.i_lne = np.concatenate([cell.i_lne for cell in cells])
    mesh.i_hsw = np.concatenate([cell.i_hsw for cell in cells])
    mesh.i_hse = np.concatenate([cell.i_hse for cell in cells])
    mesh.i_hnw = np.concatenate([cell.i_hnw for cell in cells])
    mesh.i_hne = np.concatenate([cell.i_hne for cell in cells])

    mesh.i_hw = np.concatenate([cell.i_hw for cell in cells])
    mesh.i_he = np.concatenate([cell.i_he for cell in cells])
    mesh.i_hs = np.concatenate([cell.i_hs for cell in cells])
    mesh.i_hn = np.concatenate([cell.i_hn for cell in cells])
    mesh.i_lw = np.concatenate([cell.i_lw for cell in cells])
    mesh.i_le = np.concatenate([cell.i_le for cell in cells])
    mesh.i_ls = np.concatenate([cell.i_ls for cell in cells])
    mesh.i_ln = np.concatenate([cell.i_ln for cell in cells])
    mesh.i_ne = np.concatenate([cell.i_ne for cell in cells])
    mesh.i_sw = np.concatenate([cell.i_sw for cell in cells])
    mesh.i_se = np.concatenate([cell.i_se for cell in cells])
    mesh.i_nw = np.concatenate([cell.i_nw for cell in cells])

    mesh.i_w = np.concatenate([cell.i_w for cell in cells])
    mesh.i_e = np.unique(np.concatenate([cell.i_e for cell in cells]))
    mesh.i_s = np.concatenate([cell.i_s for cell in cells])
    mesh.i_n = np.concatenate([cell.i_n for cell in cells])
    mesh.i_h = np.concatenate([cell.i_h for cell in cells])
    mesh.i_l = np.concatenate([cell.i_l for cell in cells])

    mesh.i = np.concatenate([cell.i for cell in cells])
    mesh.i.sort()
    mesh.i_all = np.concatenate([cell.i_all for cell in cells])
    mesh.i_all.sort()
    mesh.to_stim = np.concatenate([cell.to_stim for cell in cells])
    mesh.to_stim.sort()

    # only take the north and east junctions so that the junction points are only counted once
    # junctions in x direction
    mesh.gx_lse = np.concatenate([cell.gx_lse for cell in cells])
    mesh.gx_lne = np.concatenate([cell.gx_lne for cell in cells])
    mesh.gx_hse = np.concatenate([cell.gx_hse for cell in cells])
    mesh.gx_hne = np.concatenate([cell.gx_hne for cell in cells])

    mesh.gx_he = np.concatenate([cell.gx_he for cell in cells])
    mesh.gx_le = np.concatenate([cell.gx_le for cell in cells])
    mesh.gx_ne = np.concatenate([cell.gx_ne for cell in cells])
    mesh.gx_se = np.concatenate([cell.gx_se for cell in cells])

    mesh.gx_e = np.concatenate([cell.gx_e for cell in cells])

    # junctions in y direction
    mesh.gy_lnw = np.concatenate([cell.gy_lnw for cell in cells])
    mesh.gy_lne = np.concatenate([cell.gy_lne for cell in cells])
    mesh.gy_hnw = np.concatenate([cell.gy_hnw for cell in cells])
    mesh.gy_hne = np.concatenate([cell.gy_hne for cell in cells])

    mesh.gy_hn = np.concatenate([cell.gy_hn for cell in cells])
    mesh.gy_ln = np.concatenate([cell.gy_ln for cell in cells])
    mesh.gy_ne = np.concatenate([cell.gy_ne for cell in cells])
    mesh.gy_nw = np.concatenate([cell.gy_nw for cell in cells])

    mesh.gy_n = np.concatenate([cell.gy_n for cell in cells])

    # inner points for junctions in x direction
    mesh_gxi_lsw = np.concatenate([cell.gxi_lsw for cell in cells])
    mesh_gxi_lse = np.concatenate([cell.gxi_lse for cell in cells])
    mesh_gxi_lnw = np.concatenate([cell.gxi_lnw for cell in cells])
    mesh_gxi_lne = np.concatenate([cell.gxi_lne for cell in cells])
    mesh_gxi_hsw = np.concatenate([cell.gxi_hsw for cell in cells])
    mesh_gxi_hse = np.concatenate([cell.gxi_hse for cell in cells])
    mesh_gxi_hnw = np.concatenate([cell.gxi_hnw for cell in cells])
    mesh_gxi_hne = np.concatenate([cell.gxi_hne for cell in cells])

    mesh.gxi_lsw = mesh_gxi_lsw
    mesh.i_ls = np.concatenate([mesh.i_ls, mesh_gxi_lse])
    mesh.gxi_lnw = mesh_gxi_lnw
    mesh.i_ln = np.concatenate([mesh.i_ln, mesh_gxi_lne])
    mesh.gxi_hsw = mesh_gxi_hsw
    mesh.i_hs = np.concatenate([mesh.i_hs, mesh_gxi_hse])
    mesh.gxi_hnw = mesh_gxi_hnw
    mesh.i_hn = np.concatenate([mesh.i_hn, mesh_gxi_hne])

    mesh_gxi_lw = np.concatenate([cell.gxi_lw for cell in cells])
    mesh_gxi_le = np.concatenate([cell.gxi_le for cell in cells])
    mesh_gxi_hw = np.concatenate([cell.gxi_hw for cell in cells])
    mesh_gxi_he = np.concatenate([cell.gxi_he for cell in cells])
    mesh_gxi_sw = np.concatenate([cell.gxi_sw for cell in cells])
    mesh_gxi_se = np.concatenate([cell.gxi_se for cell in cells])
    mesh_gxi_nw = np.concatenate([cell.gxi_nw for cell in cells])
    mesh_gxi_ne = np.concatenate([cell.gxi_ne for cell in cells])

    mesh.gxi_lw = mesh_gxi_lw
    mesh.i_l = np.concatenate([mesh.i_l, mesh_gxi_le])
    mesh.gxi_hw = mesh_gxi_hw
    mesh.i_h = np.concatenate([mesh.i_h, mesh_gxi_he])
    mesh.gxi_sw = mesh_gxi_sw
    mesh.i_s = np.concatenate([mesh.i_s, mesh_gxi_se])
    mesh.gxi_nw = mesh_gxi_nw
    mesh.i_n = np.concatenate([mesh.i_n, mesh_gxi_ne])

    # inner points for junctions in y direction
    mesh_gyi_lsw = np.concatenate([cell.gyi_lsw for cell in cells])
    mesh_gyi_lse = np.concatenate([cell.gyi_lse for cell in cells])
    mesh_gyi_lnw = np.concatenate([cell.gyi_lnw for cell in cells])
    mesh_gyi_lne = np.concatenate([cell.gyi_lne for cell in cells])
    mesh_gyi_hsw = np.concatenate([cell.gyi_hsw for cell in cells])
    mesh_gyi_hse = np.concatenate([cell.gyi_hse for cell in cells])
    mesh_gyi_hnw = np.concatenate([cell.gyi_hnw for cell in cells])
    mesh_gyi_hne = np.concatenate([cell.gyi_hne for cell in cells])

    mesh.gyi_lsw = mesh_gyi_lsw
    mesh.gyi_lse = mesh_gyi_lse
    mesh.i_lw = np.concatenate([mesh.i_lw, mesh_gyi_lnw])
    mesh.i_le = np.concatenate([mesh.i_le, mesh_gyi_lne])
    mesh.gyi_hsw = mesh_gyi_hsw
    mesh.gyi_hse = mesh_gyi_hse
    mesh.i_hw = np.concatenate([mesh.i_hw, mesh_gyi_hnw])
    mesh.i_he = np.concatenate([mesh.i_he, mesh_gyi_hne])

    mesh_gyi_ls = np.concatenate([cell.gyi_ls for cell in cells])
    mesh_gyi_ln = np.concatenate([cell.gyi_ln for cell in cells])
    mesh_gyi_hs = np.concatenate([cell.gyi_hs for cell in cells])
    mesh_gyi_hn = np.concatenate([cell.gyi_hn for cell in cells])
    mesh_gyi_sw = np.concatenate([cell.gyi_sw for cell in cells])
    mesh_gyi_nw = np.concatenate([cell.gyi_nw for cell in cells])
    mesh_gyi_se = np.concatenate([cell.gyi_se for cell in cells])
    mesh_gyi_ne = np.concatenate([cell.gyi_ne for cell in cells])

    mesh.gyi_ls = mesh_gyi_ls # mesh.gxi_lw = mesh_gxi_lw
    mesh.i_l = np.concatenate([mesh.i_l, mesh_gyi_ln]) # mesh_gyi_ln
    mesh.gyi_hs = mesh_gyi_hs # mesh.gxi_hw = mesh_gxi_hw
    mesh.i_h = np.concatenate([mesh.i_h, mesh_gyi_hn]) # mesh_gyi_hn
    mesh.gyi_sw = mesh_gyi_sw # mesh.gxi_sw = mesh_gxi_sw
    mesh.i_w = np.concatenate([mesh.i_w, mesh_gyi_nw]) # mesh_gyi_nw
    mesh.gyi_se = mesh_gyi_se # mesh.gxi_nw = mesh_gxi_nw
    mesh.i_e = np.concatenate([mesh.i_e, mesh_gyi_ne]) # mesh_gyi_ne

    mesh.gxi_w = np.concatenate([cell.gxi_w for cell in cells])
    mesh.gyi_s = np.concatenate([cell.gyi_s for cell in cells])
    mesh.i = np.concatenate([mesh.i]
        + [cell.gxi_e for cell in cells]
        + [cell.gyi_n for cell in cells]
    )

    return mesh, cells
