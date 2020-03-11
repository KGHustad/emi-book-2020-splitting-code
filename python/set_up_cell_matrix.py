import numpy as np
import scipy.sparse

def set_up_cell_matrix(G, mesh, cell):
    "Set up matrix for the finite difference equations for a cell"

    # Load parameters
    N = G.N
    Nx = G.Nx
    Ny = G.Ny
    dx = G.dx
    dy = G.dy
    dz= G.dz
    Cm = G.Cm
    dt = G.dt
    sigma_i = G.sigma_i

    # Load grid points
    m_lsw = cell.m_lsw
    m_lse = cell.m_lse
    m_lnw = cell.m_lnw
    m_lne = cell.m_lne
    m_hsw = cell.m_hsw
    m_hse = cell.m_hse
    m_hnw = cell.m_hnw
    m_hne = cell.m_hne

    m_hw = cell.m_hw
    m_he = cell.m_he
    m_hs = cell.m_hs
    m_hn = cell.m_hn
    m_lw = cell.m_lw
    m_le = cell.m_le
    m_ls = cell.m_ls
    m_ln = cell.m_ln
    m_ne = cell.m_ne
    m_sw = cell.m_sw
    m_se = cell.m_se
    m_nw = cell.m_nw

    m_w = cell.m_w
    m_e = cell.m_e
    m_s = cell.m_s
    m_n = cell.m_n
    m_h = cell.m_h
    m_l = cell.m_l

    i = cell.i

    gx_e = cell.gx_e
    gy_n = cell.gy_n
    gx_w = cell.gx_w
    gy_s = cell.gy_s

    vec = np.zeros(N, dtype=np.float64)
    vec_kp = np.zeros(N, dtype=np.float64)
    vec_km = np.zeros(N, dtype=np.float64)
    vec_jp = np.zeros(N, dtype=np.float64)
    vec_jm = np.zeros(N, dtype=np.float64)
    vec_qp = np.zeros(N, dtype=np.float64)
    vec_qm = np.zeros(N, dtype=np.float64)

    stride_x = mesh[1, 0, 0] - mesh[0, 0, 0]
    stride_y = mesh[0, 1, 0] - mesh[0, 0, 0]
    stride_z = mesh[0, 0, 1] - mesh[0, 0, 0]


    ###### INTRACELLULAR DOMAIN #######
    index = i
    vec[index] = -((sigma_i+sigma_i)/(dx*dx) + (sigma_i+sigma_i)/(dy*dy) + (sigma_i+sigma_i)/(dz*dz))
    vec_kp[index+stride_x] = sigma_i/(dx*dx)
    vec_km[index-stride_x] = sigma_i/(dx*dx)
    vec_jp[index+stride_y] = sigma_i/(dy*dy)
    vec_jm[index-stride_y] = sigma_i/(dy*dy)
    vec_qp[index+stride_z] = sigma_i/(dz*dz)
    vec_qm[index-stride_z] = sigma_i/(dz*dz)

    ###### MEMBRANE #######

    # A: Planes

    # 1a) Set up factors for the low membrane
    index = m_l
    vec[index] = sigma_i/(dz)
    vec_qp[index+stride_z] = -sigma_i/(dz)

    # 2a) Set up factors for the high membrane
    index = m_h
    vec[index] = sigma_i/(dz)
    vec_qm[index-stride_z] = -sigma_i/(dz)

    # 3a) Set up factors for the south membrane
    index = m_s
    vec[index] = sigma_i/(dy)
    vec_jp[index+stride_y] = -sigma_i/(dy)

    # 4a) Set up factors for the north membrane
    index = m_n
    vec[index] = sigma_i/(dy)
    vec_jm[index-stride_y] = -sigma_i/(dy)

    # 5a) Set up factors for the left membrane
    index = m_w
    vec[index] = sigma_i/(dx)
    vec_kp[index+stride_x] = -sigma_i/(dx)

    # 6a) Set up factors for the right membrane
    index = m_e
    vec[index] = sigma_i/(dx)
    vec_km[index-stride_x] = -sigma_i/(dx)

    # B: Lines

    # 1b) Set up factors for the high left membrane
    index = m_hw
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dx))
    vec_qm[index-stride_z] = -0.5*sigma_i/(dz)
    vec_kp[index+stride_x] = -0.5*sigma_i/(dx)

    # 2b) Set up factors for the high right membrane
    index = m_he
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dx))
    vec_qm[index-stride_z] = -0.5*sigma_i/(dz)
    vec_km[index-stride_x] = -0.5*sigma_i/(dx)

    # 3b) Set up factors for the high south membrane
    index = m_hs
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dy))
    vec_qm[index-stride_z] = -0.5*sigma_i/(dz)
    vec_jp[index+stride_y] = -0.5*sigma_i/(dy)

    # 4b) Set up factors for the high north membrane
    index = m_hn
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dy))
    vec_qm[index-stride_z] = -0.5*sigma_i/(dz)
    vec_jm[index-stride_y] = -0.5*sigma_i/(dy)

    # 5b) Set up factors for the low left membrane
    index = m_lw
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dx))
    vec_qp[index+stride_z] = -0.5*sigma_i/(dz)
    vec_kp[index+stride_x] = -0.5*sigma_i/(dx)

    # 6b) Set up factors for the low right membrane
    index = m_le
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dx))
    vec_qp[index+stride_z] = -0.5*sigma_i/(dz)
    vec_km[index-stride_x] = -0.5*sigma_i/(dx)

    # 7b) Set up factors for the low south membrane
    index = m_ls
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dy))
    vec_qp[index+stride_z] = -0.5*sigma_i/(dz)
    vec_jp[index+stride_y] = -0.5*sigma_i/(dy)

    # 8b) Set up factors for the low north membrane
    index = m_ln
    vec[index] = 0.5*(sigma_i/(dz)+sigma_i/(dy))
    vec_qp[index+stride_z] = -0.5*sigma_i/(dz)
    vec_jm[index-stride_y] = -0.5*sigma_i/(dy)

    # 9b) Set up factors for the north left membrane
    index = m_nw
    vec[index] = 0.5*(sigma_i/(dy)+sigma_i/(dx))
    vec_jm[index-stride_y] = -0.5*sigma_i/(dy)
    vec_kp[index+stride_x] = -0.5*sigma_i/(dx)

    # 10b) Set up factors for the north right membrane
    index = m_ne
    vec[index] = 0.5*(sigma_i/(dy)+sigma_i/(dx))
    vec_jm[index-stride_y] = -0.5*sigma_i/(dy)
    vec_km[index-stride_x] = -0.5*sigma_i/(dx)

    # 11b) Set up factors for the south left membrane
    index = m_sw
    vec[index] = 0.5*(sigma_i/(dy)+sigma_i/(dx))
    vec_jp[index+stride_y] = -0.5*sigma_i/(dy)
    vec_kp[index+stride_x] = -0.5*sigma_i/(dx)

    # 12b) Set up factors for the south east membrane
    index = m_se
    vec[index] = 0.5*(sigma_i/(dy)+sigma_i/(dx))
    vec_jp[index+stride_y] = -0.5*sigma_i/(dy)
    vec_km[index-stride_x] = -0.5*sigma_i/(dx)

    # C: Corners

    # 1c) Set up factors for the lower, south, left membrane
    index = m_lsw
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qp[index+stride_z] = -1/3*sigma_i/(dz)
    vec_jp[index+stride_y] = -1/3*sigma_i/(dy)
    vec_kp[index+stride_x] = -1/3*sigma_i/(dx)

    # 2c) Set up factors for the lower, south, east membrane
    index = m_lse
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qp[index+stride_z] = -1/3*sigma_i/(dz)
    vec_jp[index+stride_y] = -1/3*sigma_i/(dy)
    vec_km[index-stride_x] = -1/3*sigma_i/(dx)

    # 3c) Set up factors for the lower, north, left membrane
    index = m_lnw
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qp[index+stride_z] = -1/3*sigma_i/(dz)
    vec_jm[index-stride_y] = -1/3*sigma_i/(dy)
    vec_kp[index+stride_x] = -1/3*sigma_i/(dx)

    # 4c) Set up factors for the lower, north, right membrane
    index = m_lne
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qp[index+stride_z] = -1/3*sigma_i/(dz)
    vec_jm[index-stride_y] = -1/3*sigma_i/(dy)
    vec_km[index-stride_x] = -1/3*sigma_i/(dx)

    # 5c) Set up factors for the higher, south, left membrane
    index = m_hsw
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qm[index-stride_z] = -1/3*sigma_i/(dz)
    vec_jp[index+stride_y] = -1/3*sigma_i/(dy)
    vec_kp[index+stride_x] = -1/3*sigma_i/(dx)

    # 6c) Set up factors for the higher, south, east membrane
    index = m_hse
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qm[index-stride_z] = -1/3*sigma_i/(dz)
    vec_jp[index+stride_y] = -1/3*sigma_i/(dy)
    vec_km[index-stride_x] = -1/3*sigma_i/(dx)

    # 7c) Set up factors for the higher, north, left membrane
    index = m_hnw
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qm[index-stride_z] = -1/3*sigma_i/(dz)
    vec_jm[index-stride_y] = -1/3*sigma_i/(dy)
    vec_kp[index+stride_x] = -1/3*sigma_i/(dx)

    # 8c) Set up factors for the higher, north, right membrane
    index = m_hne
    vec[index] = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx))
    vec_qm[index-stride_z] = -1/3*sigma_i/(dz)
    vec_jm[index-stride_y] = -1/3*sigma_i/(dy)
    vec_km[index-stride_x] = -1/3*sigma_i/(dx)

    # For all the membrane points
    index = cell.v
    vec[index] = vec[index] + Cm/dt


    ####### GAP JUNCTION (Neumann boundary condition) #######
    # Gap junction in x-direction
    index = gx_w
    vec[index] = sigma_i/(dx)
    vec_kp[index+stride_x] = -sigma_i/(dx)

    index = gx_e
    vec[index] = -sigma_i/(dx)
    vec_km[index-stride_x] = sigma_i/(dx)

    # Gap junction in y-direction
    index = gy_s
    vec[index] = sigma_i/(dy)
    vec_jp[index+stride_y] = -sigma_i/(dy)

    index = gy_n
    vec[index] = -sigma_i/(dy)
    vec_jm[index-stride_y] = sigma_i/(dy)


    ####### SET UP THE MATRIX #######
    A = scipy.sparse.spdiags(vec, 0, N, N)
    A += scipy.sparse.spdiags(vec_kp, stride_x, N, N)
    A += scipy.sparse.spdiags(vec_km, -stride_x, N, N)
    A += scipy.sparse.spdiags(vec_jp, stride_y, N, N)
    A += scipy.sparse.spdiags(vec_jm, -stride_y, N, N)
    A += scipy.sparse.spdiags(vec_qp, stride_z, N, N)
    A += scipy.sparse.spdiags(vec_qm, -stride_z, N, N)

    # Convert matrix to CSR
    A = A.tocsr()

    # Reshape matrix to a single cell
    A = A[cell.c_all, :]
    A = A[:, cell.c_all]

    A.sort_indices()

    return A

