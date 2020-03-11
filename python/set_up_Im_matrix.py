import numpy as np
import scipy.sparse

def set_up_Im_matrix(G, mesh):
    "Set up matrix used to extract the membrane current I_m from the intracellular potential"
    # TODO: Rewrite this to reuse the arrays better, like how it's done in set_up_symmetric_extracellular_matrix

    # Load parameters
    N = G.N
    Nx = G.Nx
    Ny = G.Ny
    dx = G.dx
    dy = G.dy
    dz = G.dz

    sigma_i = G.sigma_i

    m_lsw = mesh.m_lsw
    m_lse = mesh.m_lse
    m_lnw = mesh.m_lnw
    m_lne = mesh.m_lne
    m_hsw = mesh.m_hsw
    m_hse = mesh.m_hse
    m_hnw = mesh.m_hnw
    m_hne = mesh.m_hne

    m_hw = mesh.m_hw
    m_he = mesh.m_he
    m_hs = mesh.m_hs
    m_hn = mesh.m_hn
    m_lw = mesh.m_lw
    m_le = mesh.m_le
    m_ls = mesh.m_ls
    m_ln = mesh.m_ln
    m_ne = mesh.m_ne
    m_sw = mesh.m_sw
    m_se = mesh.m_se
    m_nw = mesh.m_nw

    m_w = mesh.m_w
    m_e = mesh.m_e
    m_s = mesh.m_s
    m_n = mesh.m_n
    m_h = mesh.m_h
    m_l = mesh.m_l

    stride_x = mesh[1, 0, 0] - mesh[0, 0, 0]
    stride_y = mesh[0, 1, 0] - mesh[0, 0, 0]
    stride_z = mesh[0, 0, 1] - mesh[0, 0, 0]

    # A: Planes

    # 1a) Set up factors for the low membrane
    index = m_l
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    M = scipy.sparse.spdiags(-sigma_i/dz*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dz*i_vec_qp, stride_z, N, N)


    # 2a) Set up factors for the high membrane
    index = m_h
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    M = M + scipy.sparse.spdiags(-sigma_i/dz*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dz*i_vec_qm, -stride_z, N, N)


    # 3a) Set up factors for the south membrane
    index = m_s
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    M = M + scipy.sparse.spdiags(-sigma_i/dy*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dy*i_vec_jp, stride_y, N, N)


    # 4a) Set up factors for the north membrane
    index = m_n
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    M = M + scipy.sparse.spdiags(-sigma_i/dy*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dy*i_vec_jm, -stride_y, N, N)

    # 5a) Set up factors for the left membrane
    index = m_w
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-sigma_i/dx*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 6a) Set up factors for the right membrane
    index = m_e
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-sigma_i/dx*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(sigma_i/dx*i_vec_km, -stride_x, N, N)


    # B: Edges

    # 1b) Set up factors for the high left membrane
    index = m_hw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_kp, stride_x, N, N)


    # 2b) Set up factors for the high right membrane
    index = m_he
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_km, -stride_x, N, N)


    # 3b) Set up factors for the high south membrane
    index = m_hs
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dy)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jp, stride_y, N, N)


    # 4b) Set up factors for the high north membrane
    index = m_hn
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dy)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jm, -stride_y, N, N)


    # 5b) Set up factors for the low left membrane
    index = m_lw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_kp, stride_x, N, N)


    # 6b) Set up factors for the low right membrane
    index = m_le
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_km, -stride_x, N, N)


    # 7b) Set up factors for the low south membrane
    index = m_ls
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dy)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jp, stride_y, N, N)


    # 8b) Set up factors for the low north membrane
    index = m_ln
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dz+sigma_i/dy)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jm, -stride_y, N, N)

    # 9b) Set up factors for the north left membrane
    index = m_nw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 10b) Set up factors for the north right membrane
    index = m_ne
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_km, -stride_x, N, N)

    # 11b) Set up factors for the south left membrane
    index = m_sw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 12b) Set up factors for the south east membrane
    index = m_se
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-0.5*(sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(0.5*sigma_i/dx*i_vec_km, -stride_x, N, N)


    # C: Corners

    # 1c) Set up factors for the lower, south, left membrane
    index = m_lsw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 2c) Set up factors for the lower, south, east membrane
    index = m_lse
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_km, -stride_x, N, N)

    # 3c) Set up factors for the lower, north, left membrane
    index = m_lnw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 4c) Set up factors for the lower, north, right membrane
    index = m_lne
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qp = np.zeros(N, dtype=np.float64)
    i_vec_qp[index+stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qp, stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_km, -stride_x, N, N)

    # 5c) Set up factors for the higher, south, left membrane
    index = m_hsw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 6c) Set up factors for the higher, south, east membrane
    index = m_hse
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jp = np.zeros(N, dtype=np.float64)
    i_vec_jp[index+stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jp, stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_km, -stride_x, N, N)

    # 7c) Set up factors for the higher, north, left membrane
    index = m_hnw
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_kp = np.zeros(N, dtype=np.float64)
    i_vec_kp[index+stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_kp, stride_x, N, N)

    # 8c) Set up factors for the higher, north, right membrane
    index = m_hne
    i_vec = np.zeros(N, dtype=np.float64)
    i_vec[index] = 1
    i_vec_qm = np.zeros(N, dtype=np.float64)
    i_vec_qm[index-stride_z] = 1
    i_vec_jm = np.zeros(N, dtype=np.float64)
    i_vec_jm[index-stride_y] = 1
    i_vec_km = np.zeros(N, dtype=np.float64)
    i_vec_km[index-stride_x] = 1
    M = M + scipy.sparse.spdiags(-1/3*(sigma_i/dz+sigma_i/dy+sigma_i/dx)*i_vec, 0, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dz*i_vec_qm, -stride_z, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dy*i_vec_jm, -stride_y, N, N)
    M = M + scipy.sparse.spdiags(1/3*sigma_i/dx*i_vec_km, -stride_x, N, N)

    # Convert to CSR
    M = M.tocsr()

    # Reshape matrix to the intracellular domain
    M = M[mesh.i_all, :]
    M = M[:, mesh.i_all]

    return M
