import numpy as np
import scipy.sparse

def set_up_symmetric_extracellular_matrix(G, mesh):
    "Set up a symmetric matrix for the finite difference equations of the extracellular system"

    # Load parameters
    N = G.N
    Nx = G.Nx
    Ny = G.Ny
    dx = G.dx
    dy = G.dy
    dz= G.dz
    sigma_e = G.sigma_e

    # Load mesh
    e_w = mesh.e_w
    e_e = mesh.e_e
    e_s = mesh.e_s
    e_n = mesh.e_n
    e_h = mesh.e_h
    e_l = mesh.e_l

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

    e = mesh.e

    gx_he = mesh.gx_he
    gx_le = mesh.gx_le
    gx_se = mesh.gx_se
    gx_ne = mesh.gx_ne

    gy_hn = mesh.gy_hn
    gy_ln = mesh.gy_ln
    gy_nw = mesh.gy_nw
    gy_ne = mesh.gy_ne

    gx_lse = mesh.gx_lse
    gx_lne = mesh.gx_lne
    gx_hse = mesh.gx_hse
    gx_hne = mesh.gx_hne

    gy_lnw = mesh.gy_lnw
    gy_lne = mesh.gy_lne
    gy_hnw = mesh.gy_hnw
    gy_hne = mesh.gy_hne

    stride_x = mesh[1, 0, 0] - mesh[0, 0, 0]
    stride_y = mesh[0, 1, 0] - mesh[0, 0, 0]
    stride_z = mesh[0, 0, 1] - mesh[0, 0, 0]

    """
    % Load parameters
    N = G.N;
    Nx = G.Nx;
    Ny = G.Ny;
    dx = G.dx;
    dy = G.dy;
    dz= G.dz;
    sigma_e = G.sigma_e;

    % Load mesh
    e_w = mesh.e_w;
    e_e = mesh.e_e;
    e_s = mesh.e_s;
    e_n = mesh.e_n;
    e_h = mesh.e_h;
    e_l = mesh.e_l;

    m_lsw = mesh.m_lsw;
    m_lse = mesh.m_lse;
    m_lnw = mesh.m_lnw;
    m_lne = mesh.m_lne;
    m_hsw = mesh.m_hsw;
    m_hse = mesh.m_hse;
    m_hnw = mesh.m_hnw;
    m_hne = mesh.m_hne;

    m_hw = mesh.m_hw;
    m_he = mesh.m_he;
    m_hs = mesh.m_hs;
    m_hn = mesh.m_hn;
    m_lw = mesh.m_lw;
    m_le = mesh.m_le;
    m_ls = mesh.m_ls;
    m_ln = mesh.m_ln;
    m_ne = mesh.m_ne;
    m_sw = mesh.m_sw;
    m_se = mesh.m_se;
    m_nw = mesh.m_nw;

    m_w = mesh.m_w;
    m_e = mesh.m_e;
    m_s = mesh.m_s;
    m_n = mesh.m_n;
    m_h = mesh.m_h;
    m_l = mesh.m_l;

    e = mesh.e;

    gx_he = mesh.gx_he;
    gx_le = mesh.gx_le;
    gx_se = mesh.gx_se;
    gx_ne = mesh.gx_ne;

    gy_hn = mesh.gy_hn;
    gy_ln = mesh.gy_ln;
    gy_nw = mesh.gy_nw;
    gy_ne = mesh.gy_ne;

    gx_lse = mesh.gx_lse;
    gx_lne = mesh.gx_lne;
    gx_hse = mesh.gx_hse;
    gx_hne = mesh.gx_hne;

    gy_lnw = mesh.gy_lnw;
    gy_lne = mesh.gy_lne;
    gy_hnw = mesh.gy_hnw;
    gy_hne = mesh.gy_hne;
    """

    vec = np.zeros(N, dtype=np.float64)
    vec_kp = np.zeros(N, dtype=np.float64)
    vec_km = np.zeros(N, dtype=np.float64)
    vec_jp = np.zeros(N, dtype=np.float64)
    vec_jm = np.zeros(N, dtype=np.float64)
    vec_qp = np.zeros(N, dtype=np.float64)
    vec_qm = np.zeros(N, dtype=np.float64)

    """
    vec = zeros(N, 1);
    vec_kp = zeros(N, 1);
    vec_km = zeros(N, 1);
    vec_jp = zeros(N, 1);
    vec_jm = zeros(N, 1);
    vec_qp = zeros(N, 1);
    vec_qm = zeros(N, 1);
    """

    # EXTRACELLULAR DOMAIN

    # Neumann boundary conditions
    # 1a) Set up rows for the extracellular low boundary
    index = e_l
    vec[index] = -sigma_e/(dz*dz)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # 2a) Set up rows for the extracellular high boundary
    index = e_h
    vec[index] = -sigma_e/(dz*dz)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # 3a) Set up rows for the extracellular south boundary
    index = e_s
    vec[index] = -sigma_e/(dy*dy)
    vec_jp[index+stride_y] = sigma_e/(dy*dy)

    # 4a) Set up rows for the extracellular north boundary
    index = e_n
    vec[index] = -sigma_e/(dy*dy)
    vec_jm[index-stride_y] = sigma_e/(dy*dy)

    # 5a) Set up rows for the extracellular left boundary
    index = e_w
    vec[index] = -sigma_e/(dx*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 6a) Set up rows for the extracellular right boundary
    index = e_e;
    vec[index] = -sigma_e/(dx*dx);
    vec_km[index-stride_x] = sigma_e/(dx*dx);

    """
    %%%%%%% EXTRACELLULAR DOMAIN %%%%%%%

    % Neumann boundary conditions
    % 1a) Set up rows for the extracellular low boundary
    index = e_l;
    vec(index) = -sigma_e/(dz*dz);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % 2a) Set up rows for the extracellular high boundary
    index = e_h;
    vec(index) = -sigma_e/(dz*dz);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % 3a) Set up rows for the extracellular south boundary
    index = e_s;
    vec(index) = -sigma_e/(dy*dy);
    vec_jp(index+Nx) = sigma_e/(dy*dy);

    % 4a) Set up rows for the extracellular north boundary
    index = e_n;
    vec(index) = -sigma_e/(dy*dy);
    vec_jm(index-Nx) = sigma_e/(dy*dy);

    % 5a) Set up rows for the extracellular left boundary
    index = e_w;
    vec(index) = -sigma_e/(dx*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 6a) Set up rows for the extracellular right boundary
    index = e_e;
    vec(index) = -sigma_e/(dx*dx);
    vec_km(index-1) = sigma_e/(dx*dx);
    """

    # Set up rows for the inner extracellular domain
    index = e
    vec[index] = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) + (sigma_e+sigma_e)/(dz*dz))
    vec_kp[index+stride_x] = sigma_e/(dx*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dy)
    vec_jm[index-stride_y] = sigma_e/(dy*dy)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    """
    % Set up rows for the inner extracellular domain
    index = e;
    vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                        + (sigma_e+sigma_e)/(dz*dz));
    vec_kp(index+1) = sigma_e/(dx*dx);
    vec_km(index-1) = sigma_e/(dx*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dy);
    vec_jm(index-Nx) = sigma_e/(dy*dy);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);
    """

    # MEMBRANE

    ## A: Sides

    # 1a) Set up factors for the low membrane
    index = m_l
    vec[index] = -sigma_e/(dz*dz)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # 2a) Set up factors for the high membrane
    index = m_h
    vec[index] = -sigma_e/(dz*dz)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # 3a) Set up factors for the south membrane
    index = m_s
    vec[index] = -sigma_e/(dy*dy)
    vec_jm[index-stride_y] = sigma_e/(dy*dy)

    # 4a) Set up factors for the north membrane
    index = m_n
    vec[index] = -sigma_e/(dy*dy)
    vec_jp[index+stride_y] = sigma_e/(dy*dy)

    # 5a) Set up factors for the left membrane
    index = m_w
    vec[index] = -sigma_e/(dx*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 6a) Set up factors for the right membrane
    index = m_e
    vec[index] = -sigma_e/(dx*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)


    ## B: Edges

    # 1b) Set up factors for the high left membrane
    index = m_hw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 2b) Set up factors for the high right membrane
    index = m_he
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 3b) Set up factors for the high south membrane
    index = m_hs
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)

    # 4b) Set up factors for the high north membrane
    index = m_hn
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)

    # 5b) Set up factors for the low left membrane
    index = m_lw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 6b) Set up factors for the low right membrane
    index = m_le
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 7b) Set up factors for the low south membrane
    index = m_ls
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)

    # 8b) Set up factors for the low north membrane
    index = m_ln
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)

    # 9b) Set up factors for the north left membrane
    index = m_nw
    vec[index] = -(sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 10b) Set up factors for the north right membrane
    index = m_ne
    vec[index] = -(sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 11b) Set up factors for the south left membrane
    index = m_sw
    vec[index] = -(sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 12b) Set up factors for the south east membrane
    index = m_se
    vec[index] = -(sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)


    ## C: Corners

    # 1c) Set up factors for the lower, south, left membrane
    index = m_lsw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 2c) Set up factors for the lower, south, east membrane
    index = m_lse
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 3c) Set up factors for the lower, north, left membrane
    index = m_lnw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 4c) Set up factors for the lower, north, right membrane
    index = m_lne
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qm[index-stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 5c) Set up factors for the higher, south, left membrane
    index = m_hsw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 6c) Set up factors for the higher, south, east membrane
    index = m_hse
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jm[index-stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # 7c) Set up factors for the higher, north, left membrane
    index = m_hnw
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # 8c) Set up factors for the higher, north, right membrane
    index = m_hne
    vec[index] = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx))
    vec_qp[index+stride_z] = sigma_e/(dz*dx)
    vec_jp[index+stride_y] = sigma_e/(dy*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)


    """
    %%%%%% MEMBRANE %%%%%%

    % 1a) Set up factors for the low membrane
    index = m_l;
    vec(index) = -sigma_e/(dz*dz);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % 2a) Set up factors for the high membrane
    index = m_h;
    vec(index) = -sigma_e/(dz*dz);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % 3a) Set up factors for the south membrane
    index = m_s;
    vec(index) = -sigma_e/(dy*dy);
    vec_jm(index-Nx) = sigma_e/(dy*dy);

    % 4a) Set up factors for the north membrane
    index = m_n;
    vec(index) = -sigma_e/(dy*dy);
    vec_jp(index+Nx) = sigma_e/(dy*dy);

    % 5a) Set up factors for the left membrane
    index = m_w;
    vec(index) = -sigma_e/(dx*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 6a) Set up factors for the right membrane
    index = m_e;
    vec(index) = -sigma_e/(dx*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 1b) Set up factors for the high left membrane
    index = m_hw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 2b) Set up factors for the high right membrane
    index = m_he;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 3b) Set up factors for the high south membrane
    index = m_hs;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);

    % 4b) Set up factors for the high north membrane
    index = m_hn;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);

    % 5b) Set up factors for the low left membrane
    index = m_lw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 6b) Set up factors for the low right membrane
    index = m_le;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 7b) Set up factors for the low south membrane
    index = m_ls;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);

    % 8b) Set up factors for the low north membrane
    index = m_ln;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);

    % 9b) Set up factors for the north left membrane
    index = m_nw;
    vec(index) = -(sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 10b) Set up factors for the north right membrane
    index = m_ne;
    vec(index) = -(sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 11b) Set up factors for the south left membrane
    index = m_sw;
    vec(index) = -(sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 12b) Set up factors for the south east membrane
    index = m_se;
    vec(index) = -(sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 1c) Set up factors for the lower, south, left membrane
    index = m_lsw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 2c) Set up factors for the lower, south, east membrane
    index = m_lse;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 3c) Set up factors for the lower, north, left membrane
    index = m_lnw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 4c) Set up factors for the lower, north, right membrane
    index = m_lne;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 5c) Set up factors for the higher, south, left membrane
    index = m_hsw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 6c) Set up factors for the higher, south, east membrane
    index = m_hse;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jm(index-Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % 7c) Set up factors for the higher, north, left membrane
    index = m_hnw;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % 8c) Set up factors for the higher, north, right membrane
    index = m_hne;
    vec(index) = -(sigma_e/(dz*dx)+sigma_e/(dy*dx)+sigma_e/(dx*dx));
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dx);
    vec_jp(index+Nx) = sigma_e/(dy*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);
    """

    # CORNER POINTS BETWEEN MEMBRANE AND GAP JUNCTION

    # junctions in x direction

    # South corner line
    index = gx_se
    vec[index] = -sigma_e/(dy*dy)
    vec_jm[index-stride_y] = sigma_e/(dy*dy)

    # North corner line
    index = gx_ne
    vec[index] = -sigma_e/(dy*dy)
    vec_jp[index+stride_y] = sigma_e/(dy*dy)

    # Lower corner line
    index = gx_le
    vec[index] = -sigma_e/(dz*dz)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper corner line
    index = gx_he
    vec[index] = -sigma_e/(dz*dz)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # Lower, south corner point
    index = gx_lse
    vec[index] = -(sigma_e/(dy*dy) + sigma_e/(dz*dz))
    vec_jm[index-stride_y] = sigma_e/(dy*dy)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper, south corner point
    index = gx_hse
    vec[index] = -(sigma_e/(dy*dy) + sigma_e/(dz*dz))
    vec_jm[index-stride_y] = sigma_e/(dy*dy)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # Lower, north corner point
    index = gx_lne
    vec[index] = -(sigma_e/(dy*dy) + sigma_e/(dz*dz))
    vec_jp[index+stride_y] = sigma_e/(dy*dy)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper, north corner point
    index = gx_hne
    vec[index] = -(sigma_e/(dy*dy) + sigma_e/(dz*dz))
    vec_jp[index+stride_y] = sigma_e/(dy*dy)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)


    # junctions in y direction

    # West corner line
    index = gy_nw
    vec[index] = -sigma_e/(dx*dx)
    vec_km[index-stride_x] = sigma_e/(dx*dx)

    # East corner line
    index = gy_ne
    vec[index] = -sigma_e/(dx*dx)
    vec_kp[index+stride_x] = sigma_e/(dx*dx)

    # Lower corner line
    index = gy_ln
    vec[index] = -sigma_e/(dz*dz)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper corner line
    index = gy_hn
    vec[index] = -sigma_e/(dz*dz)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # Lower, west corner point
    index = gy_lnw
    vec[index] = -(sigma_e/(dx*dx) + sigma_e/(dz*dz))
    vec_km[index-stride_x] = sigma_e/(dx*dx)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper, west corner point
    index = gy_hnw
    vec[index] = -(sigma_e/(dx*dx) + sigma_e/(dz*dz))
    vec_km[index-stride_x] = sigma_e/(dx*dx)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    # Lower, east corner point
    index = gy_lne
    vec[index] = -(sigma_e/(dx*dx) + sigma_e/(dz*dz))
    vec_kp[index+stride_x] = sigma_e/(dx*dx)
    vec_qm[index-stride_z] = sigma_e/(dz*dz)

    # Upper, east corner point
    index = gy_hne
    vec[index] = -(sigma_e/(dx*dx) + sigma_e/(dz*dz))
    vec_kp[index+stride_x] = sigma_e/(dx*dx)
    vec_qp[index+stride_z] = sigma_e/(dz*dz)

    """
    %%%%%%% CORNER POINTS BETWEEN MEMBRANE AND GAP JUNCTION %%%%%%%

    % South corner line
    index = gx_se;
    vec(index) = -sigma_e/(dy*dy);
    vec_jm(index-Nx) = sigma_e/(dy*dy);

    % North corner line
    index = gx_ne;
    vec(index) = -sigma_e/(dy*dy);
    vec_jp(index+Nx) = sigma_e/(dy*dy);

    % Lower corner line
    index = gx_le;
    vec(index) = -sigma_e/(dz*dz);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper corner line
    index = gx_he;
    vec(index) = -sigma_e/(dz*dz);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % Lower, south corner point
    index = gx_lse;
    vec(index) = -(sigma_e/(dy*dy) + sigma_e/(dz*dz));
    vec_jm(index-Nx) = sigma_e/(dy*dy);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper, south corner point
    index = gx_hse;
    vec(index) = -(sigma_e/(dy*dy) + sigma_e/(dz*dz));
    vec_jm(index-Nx) = sigma_e/(dy*dy);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % Lower, north corner point
    index = gx_lne;
    vec(index) = -(sigma_e/(dy*dy) + sigma_e/(dz*dz));
    vec_jp(index+Nx) = sigma_e/(dy*dy);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper, north corner point
    index = gx_hne;
    vec(index) = -(sigma_e/(dy*dy) + sigma_e/(dz*dz));
    vec_jp(index+Nx) = sigma_e/(dy*dy);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);


    % West corner line
    index = gy_nw;
    vec(index) = -sigma_e/(dx*dx);
    vec_km(index-1) = sigma_e/(dx*dx);

    % East corner line
    index = gy_ne;
    vec(index) = -sigma_e/(dx*dx);
    vec_kp(index+1) = sigma_e/(dx*dx);

    % Lower corner line
    index = gy_ln;
    vec(index) = -sigma_e/(dz*dz);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper corner line
    index = gy_hn;
    vec(index) = -sigma_e/(dz*dz);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % Lower, west corner point
    index = gy_lnw;
    vec(index) = -(sigma_e/(dx*dx) + sigma_e/(dz*dz));
    vec_km(index-1) = sigma_e/(dx*dx);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper, west corner point
    index = gy_hnw;
    vec(index) = -(sigma_e/(dx*dx) + sigma_e/(dz*dz));
    vec_km(index-1) = sigma_e/(dx*dx);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);

    % Lower, east corner point
    index = gy_lne;
    vec(index) = -(sigma_e/(dx*dx) + sigma_e/(dz*dz));
    vec_kp(index+1) = sigma_e/(dx*dx);
    vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

    % Upper, east corner point
    index = gy_hne;
    vec(index) = -(sigma_e/(dx*dx) + sigma_e/(dz*dz));
    vec_kp(index+1) = sigma_e/(dx*dx);
    vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
    """

    A = scipy.sparse.spdiags(vec, 0, N, N)
    A += scipy.sparse.spdiags(vec_kp, stride_x, N, N)
    A += scipy.sparse.spdiags(vec_km, -stride_x, N, N)
    A += scipy.sparse.spdiags(vec_jp, stride_y, N, N)
    A += scipy.sparse.spdiags(vec_jm, -stride_y, N, N)
    A += scipy.sparse.spdiags(vec_qp, stride_z, N, N)
    A += scipy.sparse.spdiags(vec_qm, -stride_z, N, N)
    """
    A = scipy.sparse.diags(
        #[
        #    vec,
        #    vec_km,
        #    vec_kp,
        #    vec_jm,
        #    vec_jp,
        #    vec_qm,
        #    vec_qp,
        #],
        [
            vec,
            vec_km[stride_x:],
            vec_kp[:-stride_x],
            vec_jm[stride_y:],
            vec_jp[:-stride_y],
            vec_qm[stride_z:],
            vec_qp[:-stride_z],
        ],
        offsets=np.array([
            0,
            -stride_x,
            stride_x,
            -stride_y,
            stride_y,
            -stride_z,
            stride_z,
        ]),
        shape=(N,N)
    )
    """
    #print(A.shape)

    A = A.tocsr()
    #print(A.shape)
    #print(A.nnz)
    # local numbering in mesh.e_all
    A = A[mesh.e_all, :]
    #print(A.nnz)
    A = A[:, mesh.e_all]
    #print(A.nnz)

    A = -A

    A.sort_indices()

    return A

    """
    %%%%%%% SET UP THE MATRIX %%%%%%%
    A = spdiags(vec, 0, N, N);
    A = A + spdiags(vec_kp, 1, N, N);
    A = A + spdiags(vec_km, -1, N, N);
    A = A + spdiags(vec_jp, Nx, N, N);
    A = A + spdiags(vec_jm, -Nx, N, N);
    A = A + spdiags(vec_qp, Nx*Ny, N, N);
    A = A + spdiags(vec_qm, -Nx*Ny, N, N);


    % Reshape matrix to get fewer unknowns
    A = A(mesh.e_all, mesh.e_all);
    A = -A;
    """
