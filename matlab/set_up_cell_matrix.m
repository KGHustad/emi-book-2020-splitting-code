function A = set_up_cell_matrix(G, cell)
%A = set_up_cell_matrix(G, cell) 
% Set up matrix for the finite difference equations for a single cell
%
% Output arguments:
%       A: finite difference matrix
%       
% Input arguments:
%       G: object containing information about the problem
%    cell: object containing information about the mesh of a cell

% Load parameters
N = G.N; Nx = G.Nx; Ny = G.Ny; dx = G.dx; dy = G.dy; dz = G.dz;
Cm = G.Cm; dt = G.dt; sigma_i = G.sigma_i;

% Load mesh
m_lsw = cell.m_lsw;
m_lse = cell.m_lse;
m_lnw = cell.m_lnw;
m_lne = cell.m_lne;
m_hsw = cell.m_hsw;
m_hse = cell.m_hse;
m_hnw = cell.m_hnw;
m_hne = cell.m_hne;

m_hw = cell.m_hw;
m_he = cell.m_he;
m_hs = cell.m_hs;
m_hn = cell.m_hn;
m_lw = cell.m_lw;
m_le = cell.m_le;
m_ls = cell.m_ls;
m_ln = cell.m_ln;
m_ne = cell.m_ne;
m_sw = cell.m_sw;
m_se = cell.m_se;
m_nw = cell.m_nw;

m_w = cell.m_w;
m_e = cell.m_e;
m_s = cell.m_s;
m_n = cell.m_n;
m_h = cell.m_h;
m_l = cell.m_l;

i = cell.i;

gx_e = cell.gx_e;
gy_n = cell.gy_n;
gx_w = cell.gx_w;
gy_s = cell.gy_s;

vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);
vec_qp = zeros(N, 1);
vec_qm = zeros(N, 1);


%%%%%% PURELY INTRACELLULAR DOMAIN %%%%%%%
index = i;
vec(index) = -((sigma_i+sigma_i)/(dx*dx) + (sigma_i+sigma_i)/(dy*dy) + (sigma_i+sigma_i)/(dz*dz)); 
vec_kp(index+1) = sigma_i/(dx*dx);
vec_km(index-1) = sigma_i/(dx*dx);
vec_jp(index+Nx) = sigma_i/(dy*dy);
vec_jm(index-Nx) = sigma_i/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_i/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_i/(dz*dz);

%%%%%% MEMBRANE %%%%%%%

% 1a) Set up factors for the low membrane
index = m_l;
vec(index) = sigma_i/(dz);
vec_qp(index+Nx*Ny) = -sigma_i/(dz);

% 2a) Set up factors for the high membrane
index = m_h;
vec(index) = sigma_i/(dz);
vec_qm(index-Nx*Ny) = -sigma_i/(dz);

% 3a) Set up factors for the south membrane
index = m_s;
vec(index) = sigma_i/(dy);
vec_jp(index+Nx) = -sigma_i/(dy);

% 4a) Set up factors for the north membrane
index = m_n;
vec(index) = sigma_i/(dy);
vec_jm(index-Nx) = -sigma_i/(dy);

% 5a) Set up factors for the left membrane
index = m_w;
vec(index) = sigma_i/(dx);
vec_kp(index+1) = -sigma_i/(dx);

% 6a) Set up factors for the right membrane
index = m_e;
vec(index) = sigma_i/(dx);
vec_km(index-1) = -sigma_i/(dx);

% 1b) Set up factors for the high left membrane
index = m_hw;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -0.5*sigma_i/(dz);
vec_kp(index+1) = -0.5*sigma_i/(dx);

% 2b) Set up factors for the high right membrane
index = m_he;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -0.5*sigma_i/(dz);
vec_km(index-1) = -0.5*sigma_i/(dx);

% 3b) Set up factors for the high south membrane
index = m_hs;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dy));
vec_qm(index-Nx*Ny) = -0.5*sigma_i/(dz);
vec_jp(index+Nx) = -0.5*sigma_i/(dy);

% 4b) Set up factors for the high north membrane
index = m_hn;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dy));
vec_qm(index-Nx*Ny) = -0.5*sigma_i/(dz);
vec_jm(index-Nx) = -0.5*sigma_i/(dy);

% 5b) Set up factors for the low left membrane
index = m_lw;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -0.5*sigma_i/(dz);
vec_kp(index+1) = -0.5*sigma_i/(dx);

% 6b) Set up factors for the low right membrane
index = m_le;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -0.5*sigma_i/(dz);
vec_km(index-1) = -0.5*sigma_i/(dx);

% 7b) Set up factors for the low south membrane
index = m_ls;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dy));
vec_qp(index+Nx*Ny) = -0.5*sigma_i/(dz);
vec_jp(index+Nx) = -0.5*sigma_i/(dy);

% 8b) Set up factors for the low north membrane
index = m_ln;
vec(index) = 0.5*(sigma_i/(dz)+sigma_i/(dy));
vec_qp(index+Nx*Ny) = -0.5*sigma_i/(dz);
vec_jm(index-Nx) = -0.5*sigma_i/(dy);

% 9b) Set up factors for the north left membrane
index = m_nw;
vec(index) = 0.5*(sigma_i/(dy)+sigma_i/(dx));
vec_jm(index-Nx) = -0.5*sigma_i/(dy);
vec_kp(index+1) = -0.5*sigma_i/(dx);

% 10b) Set up factors for the north right membrane
index = m_ne;
vec(index) = 0.5*(sigma_i/(dy)+sigma_i/(dx));
vec_jm(index-Nx) = -0.5*sigma_i/(dy);
vec_km(index-1) = -0.5*sigma_i/(dx);

% 11b) Set up factors for the south left membrane
index = m_sw;
vec(index) = 0.5*(sigma_i/(dy)+sigma_i/(dx));
vec_jp(index+Nx) = -0.5*sigma_i/(dy);
vec_kp(index+1) = -0.5*sigma_i/(dx);

% 12b) Set up factors for the south east membrane
index = m_se;
vec(index) = 0.5*(sigma_i/(dy)+sigma_i/(dx));
vec_jp(index+Nx) = -0.5*sigma_i/(dy);
vec_km(index-1) = -0.5*sigma_i/(dx);

% 1c) Set up factors for the lower, south, left membrane
index = m_lsw;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -1/3*sigma_i/(dz);
vec_jp(index+Nx) = -1/3*sigma_i/(dy);
vec_kp(index+1) = -1/3*sigma_i/(dx);

% 2c) Set up factors for the lower, south, east membrane
index = m_lse;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -1/3*sigma_i/(dz);
vec_jp(index+Nx) = -1/3*sigma_i/(dy);
vec_km(index-1) = -1/3*sigma_i/(dx);

% 3c) Set up factors for the lower, north, left membrane
index = m_lnw;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -1/3*sigma_i/(dz);
vec_jm(index-Nx) = -1/3*sigma_i/(dy);
vec_kp(index+1) = -1/3*sigma_i/(dx);

% 4c) Set up factors for the lower, north, right membrane
index = m_lne;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qp(index+Nx*Ny) = -1/3*sigma_i/(dz);
vec_jm(index-Nx) = -1/3*sigma_i/(dy);
vec_km(index-1) = -1/3*sigma_i/(dx);

% 5c) Set up factors for the higher, south, left membrane
index = m_hsw;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -1/3*sigma_i/(dz);
vec_jp(index+Nx) = -1/3*sigma_i/(dy);
vec_kp(index+1) = -1/3*sigma_i/(dx);

% 6c) Set up factors for the higher, south, east membrane
index = m_hse;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -1/3*sigma_i/(dz);
vec_jp(index+Nx) = -1/3*sigma_i/(dy);
vec_km(index-1) = -1/3*sigma_i/(dx);

% 7c) Set up factors for the higher, north, left membrane
index = m_hnw;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -1/3*sigma_i/(dz);
vec_jm(index-Nx) = -1/3*sigma_i/(dy);
vec_kp(index+1) = -1/3*sigma_i/(dx);

% 8c) Set up factors for the higher, north, right membrane
index = m_hne;
vec(index) = 1/3*(sigma_i/(dz)+sigma_i/(dy)+sigma_i/(dx));
vec_qm(index-Nx*Ny) = -1/3*sigma_i/(dz);
vec_jm(index-Nx) = -1/3*sigma_i/(dy);
vec_km(index-1) = -1/3*sigma_i/(dx);

% For all the membrane points
index = cell.v;
vec(index) = vec(index) + Cm/dt;


%%%%%%% GAP JUNCTION (Neumann boundary condition) %%%%%%%
% Gap juntion in x-direction
index = gx_w;
vec(index) = sigma_i/(dx); 
vec_kp(index+1) = -sigma_i/(dx);

index = gx_e;
vec(index) = -sigma_i/(dx); 
vec_km(index-1) = sigma_i/(dx);

% Gap juntion in y-direction
index = gy_s;
vec(index) = sigma_i/(dy); 
vec_jp(index+Nx) = -sigma_i/(dy);

index = gy_n;
vec(index) = -sigma_i/(dy); 
vec_jm(index-Nx) = sigma_i/(dy);


%%%%%%% SET UP THE MATRIX %%%%%%%
A = spdiags(vec, 0, N, N);
A = A + spdiags(vec_kp, 1, N, N);
A = A + spdiags(vec_km, -1, N, N);
A = A + spdiags(vec_jp, Nx, N, N);
A = A + spdiags(vec_jm, -Nx, N, N);
A = A + spdiags(vec_qp, Nx*Ny, N, N);
A = A + spdiags(vec_qm, -Nx*Ny, N, N);


% Reshape matrix to get fewer unknowns
A = A(cell.c_all, cell.c_all);

end

