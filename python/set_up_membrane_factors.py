import numpy as np


def set_up_membrane_factors(G, mesh, float_t=np.float64):
    """Set up factors for the membrane currents at the corners needed to make
    the linear system symmentric"""

    membrane_factors = np.zeros(G.N, dtype=float_t)
    membrane_factors[mesh.v] = 1

    # Corner points
    corner_index_sets = [
        mesh.m_lsw,
        mesh.m_lse,
        mesh.m_lnw,
        mesh.m_lne,
        mesh.m_hsw,
        mesh.m_hse,
        mesh.m_hnw,
        mesh.m_hne,
    ]
    for index_set in corner_index_sets:
        membrane_factors[index_set] = 3

    # Membrane lines
    line_index_sets = [
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
    ]
    for index_set in line_index_sets:
        membrane_factors[index_set] = 2

    membrane_factors /= G.dx
    membrane_factors = membrane_factors[mesh.v].copy()

    return membrane_factors
