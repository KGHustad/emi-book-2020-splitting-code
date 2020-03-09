function A = set_up_extracellular_matrix(G, mesh)
%A = set_up_extracellular_matrix(G, mesh) 
% Set up matrix for the finite difference equations for the extracellular
% domain
%
% Output argument:
%       A: finite difference matrix
%
% Input arguments:
%       G: object containing information about the problem
%    mesh: object containing information about the mesh

% Load parameters
N = G.N; Nx = G.Nx; Ny = G.Ny; dx = G.dx; dy = G.dy; dz= G.dz; 
sigma_e = G.sigma_e;

% Load mesh
e_lsw = mesh.e_lsw;
e_lse = mesh.e_lse;
e_lnw = mesh.e_lnw;
e_lne = mesh.e_lne;
e_hsw = mesh.e_hsw;
e_hse = mesh.e_hse;
e_hnw = mesh.e_hnw;
e_hne = mesh.e_hne;

e_hw = mesh.e_hw;
e_he = mesh.e_he;
e_hs = mesh.e_hs;
e_hn = mesh.e_hn;
e_lw = mesh.e_lw;
e_le = mesh.e_le;
e_ls = mesh.e_ls;
e_ln = mesh.e_ln;
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;

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


vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);
vec_qp = zeros(N, 1);
vec_qm = zeros(N, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        Dirichlet boundary condition in the x-direction        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = [e_w, e_e, e_hw, e_he, e_lw, e_le, e_nw, e_ne, e_sw, e_se, ...
    e_lsw, e_lse, e_lnw, e_lne, e_hsw, e_hse, e_hnw, e_hne];
vec(index) = 1/(dx*dx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Neumann boundary condition in the y, and z-directions     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a) Set up rows for the extracellular low boundary
index = e_l;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
           + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = sigma_e/(dy*dy);
vec_jm(index-Nx) = sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = 2*sigma_e/(dz*dz);

% 2a) Set up rows for the extracellular high boundary
index = e_h;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
           + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = sigma_e/(dy*dy);
vec_jm(index-Nx) = sigma_e/(dy*dy);
vec_qm(index-Nx*Ny) = 2*sigma_e/(dz*dz);

% 3a) Set up rows for the extracellular south boundary
index = e_s;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz)); 
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = 2*sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

% 4a) Set up rows for the extracellular north boundary
index = e_n;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jm(index-Nx) = 2*sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);
    
% 3b) Set up rows for the extracellular high south boundary
index = e_hs;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = 2*sigma_e/(dy*dy);
vec_qm(index-Nx*Ny) = 2*sigma_e/(dz*dz);

% 4b) Set up rows for the extracellular high north boundary
index = e_hn;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jm(index-Nx) = 2*sigma_e/(dy*dy);
vec_qm(index-Nx*Ny) = 2*sigma_e/(dz*dz);
    
% 7b) Set up rows for the extracellular low south boundary
index = e_ls;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = 2*sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = 2*sigma_e/(dz*dz);

% 8b) Set up rows for the extracellular low north boundary
index = e_ln;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jm(index-Nx) = 2*sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = 2*sigma_e/(dz*dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                  INNER EXTRACELLULAR DOMAIN                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = e;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                    + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = sigma_e/(dy*dy);
vec_jm(index-Nx) = sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                           MEMBRANE                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a) Set up factors for the low membrane
index = m_l;
vec(index) = -sigma_e/(dz);
vec_qm(index-Nx*Ny) = sigma_e/(dz);

% 2a) Set up factors for the high membrane
index = m_h;
vec(index) = -sigma_e/(dz);
vec_qp(index+Nx*Ny) = sigma_e/(dz);

% 3a) Set up factors for the south membrane
index = m_s;
vec(index) = -sigma_e/(dy);
vec_jm(index-Nx) = sigma_e/(dy);

% 4a) Set up factors for the north membrane
index = m_n;
vec(index) = -sigma_e/(dy);
vec_jp(index+Nx) = sigma_e/(dy);

% 5a) Set up factors for the left membrane
index = m_w;
vec(index) = -sigma_e/(dx);
vec_km(index-1) = sigma_e/(dx);

% 6a) Set up factors for the right membrane
index = m_e;
vec(index) = -sigma_e/(dx);
vec_kp(index+1) = sigma_e/(dx);

% 1b) Set up factors for the high left membrane
index = m_hw;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 0.5*sigma_e/(dz);
vec_km(index-1) = 0.5*sigma_e/(dx);

% 2b) Set up factors for the high right membrane
index = m_he;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 0.5*sigma_e/(dz);
vec_kp(index+1) = 0.5*sigma_e/(dx);

% 3b) Set up factors for the high south membrane
index = m_hs;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dy));
vec_qp(index+Nx*Ny) = 0.5*sigma_e/(dz);
vec_jm(index-Nx) = 0.5*sigma_e/(dy);

% 4b) Set up factors for the high north membrane
index = m_hn;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dy));
vec_qp(index+Nx*Ny) = 0.5*sigma_e/(dz);
vec_jp(index+Nx) = 0.5*sigma_e/(dy);

% 5b) Set up factors for the low left membrane
index = m_lw;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 0.5*sigma_e/(dz);
vec_km(index-1) = 0.5*sigma_e/(dx);

% 6b) Set up factors for the low right membrane
index = m_le;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 0.5*sigma_e/(dz);
vec_kp(index+1) = 0.5*sigma_e/(dx);

% 7b) Set up factors for the low south membrane
index = m_ls;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dy));
vec_qm(index-Nx*Ny) = 0.5*sigma_e/(dz);
vec_jm(index-Nx) = 0.5*sigma_e/(dy);

% 8b) Set up factors for the low north membrane
index = m_ln;
vec(index) = -0.5*(sigma_e/(dz)+sigma_e/(dy));
vec_qm(index-Nx*Ny) = 0.5*sigma_e/(dz);
vec_jp(index+Nx) = 0.5*sigma_e/(dy);

% 9b) Set up factors for the north left membrane
index = m_nw;
vec(index) = -0.5*(sigma_e/(dy)+sigma_e/(dx));
vec_jp(index+Nx) = 0.5*sigma_e/(dy);
vec_km(index-1) = 0.5*sigma_e/(dx);

% 10b) Set up factors for the north right membrane
index = m_ne;
vec(index) = -0.5*(sigma_e/(dy)+sigma_e/(dx));
vec_jp(index+Nx) = 0.5*sigma_e/(dy);
vec_kp(index+1) = 0.5*sigma_e/(dx);

% 11b) Set up factors for the south left membrane
index = m_sw;
vec(index) = -0.5*(sigma_e/(dy)+sigma_e/(dx));
vec_jm(index-Nx) = 0.5*sigma_e/(dy);
vec_km(index-1) = 0.5*sigma_e/(dx);

% 12b) Set up factors for the south east membrane
index = m_se;
vec(index) = -0.5*(sigma_e/(dy)+sigma_e/(dx));
vec_jm(index-Nx) = 0.5*sigma_e/(dy);
vec_kp(index+1) = 0.5*sigma_e/(dx);

% 1c) Set up factors for the lower, south, left membrane
index = m_lsw;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 1/3*sigma_e/(dz);
vec_jm(index-Nx) = 1/3*sigma_e/(dy);
vec_km(index-1) = 1/3*sigma_e/(dx);

% 2c) Set up factors for the lower, south, east membrane
index = m_lse;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 1/3*sigma_e/(dz);
vec_jm(index-Nx) = 1/3*sigma_e/(dy);
vec_kp(index+1) = 1/3*sigma_e/(dx);

% 3c) Set up factors for the lower, north, left membrane
index = m_lnw;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 1/3*sigma_e/(dz);
vec_jp(index+Nx) = 1/3*sigma_e/(dy);
vec_km(index-1) = 1/3*sigma_e/(dx);

% 4c) Set up factors for the lower, north, right membrane
index = m_lne;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qm(index-Nx*Ny) = 1/3*sigma_e/(dz);
vec_jp(index+Nx) = 1/3*sigma_e/(dy);
vec_kp(index+1) = 1/3*sigma_e/(dx);

% 5c) Set up factors for the higher, south, left membrane
index = m_hsw;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 1/3*sigma_e/(dz);
vec_jm(index-Nx) = 1/3*sigma_e/(dy);
vec_km(index-1) = 1/3*sigma_e/(dx);

% 6c) Set up factors for the higher, south, east membrane
index = m_hse;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 1/3*sigma_e/(dz);
vec_jm(index-Nx) = 1/3*sigma_e/(dy);
vec_kp(index+1) = 1/3*sigma_e/(dx);

% 7c) Set up factors for the higher, north, left membrane
index = m_hnw;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 1/3*sigma_e/(dz);
vec_jp(index+Nx) = 1/3*sigma_e/(dy);
vec_km(index-1) = 1/3*sigma_e/(dx);

% 8c) Set up factors for the higher, north, right membrane
index = m_hne;
vec(index) = -1/3*(sigma_e/(dz)+sigma_e/(dy)+sigma_e/(dx));
vec_qp(index+Nx*Ny) = 1/3*sigma_e/(dz);
vec_jp(index+Nx) = 1/3*sigma_e/(dy);
vec_kp(index+1) = 1/3*sigma_e/(dx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     CORNER POINTS BETWEEN MEMBRANE AND GAP JUNCTION     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% South corner line
index = gx_se;
vec(index) = 1/(dx*dx);
vec_jm(index-Nx) = -1/(dx*dx);

% North corner line
index = gx_ne;
vec(index) = 1/(dx*dx);
vec_jp(index+Nx) = -1/(dx*dx);

% Lower corner line
index = gx_le;
vec(index) = 1/(dx*dx);
vec_qm(index-Nx*Ny) = -1/(dx*dx);

% Upper corner line
index = gx_he;
vec(index) = 1/(dx*dx);
vec_qp(index+Nx*Ny) = -1/(dx*dx);

% Lower, south corner point
index = gx_lse;
vec(index) = 1/(dx*dx);
vec_jm(index-Nx) = -0.5/(dx*dx);
vec_qm(index-Nx*Ny) = -0.5/(dx*dx);

% Upper, south corner point
index = gx_hse;
vec(index) = 1/(dx*dx);
vec_jm(index-Nx) = -0.5/(dx*dx);
vec_qp(index+Nx*Ny) = -0.5/(dx*dx);

% Lower, north corner point
index = gx_lne;
vec(index) = 1/(dx*dx);
vec_jp(index+Nx) = -0.5/(dx*dx);
vec_qm(index-Nx*Ny) = -0.5/(dx*dx);

% Upper, north corner point
index = gx_hne;
vec(index) = 1/(dx*dx);
vec_jp(index+Nx) = -0.5/(dx*dx);
vec_qp(index+Nx*Ny) = -0.5/(dx*dx);


% West corner line
index = gy_nw;
vec(index) = 1/(dx*dx);
vec_km(index-1) = -1/(dx*dx);

% East corner line
index = gy_ne;
vec(index) = 1/(dx*dx);
vec_kp(index+1) = -1/(dx*dx);

% Lower corner line
index = gy_ln;
vec(index) = 1/(dx*dx);
vec_qm(index-Nx*Ny) = -1/(dx*dx);

% Upper corner line
index = gy_hn;
vec(index) = 1/(dx*dx);
vec_qp(index+Nx*Ny) = -1/(dx*dx);

% Lower, west corner point
index = gy_lnw;
vec(index) = 1/(dx*dx);
vec_km(index-1) = -0.5/(dx*dx);
vec_qm(index-Nx*Ny) = -0.5/(dx*dx);

% Upper, west corner point
index = gy_hnw;
vec(index) = 1/(dx*dx);
vec_km(index-1) = -0.5/(dx*dx);
vec_qp(index+Nx*Ny) = -0.5/(dx*dx);

% Lower, east corner point
index = gy_lne;
vec(index) = 1/(dx*dx);
vec_kp(index+1) = -0.5/(dx*dx);
vec_qm(index-Nx*Ny) = -0.5/(dx*dx);

% Upper, east corner point
index = gy_hne;
vec(index) = 1/(dx*dx);
vec_kp(index+1) = -0.5/(dx*dx);
vec_qp(index+Nx*Ny) = -0.5/(dx*dx);


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

end

