function [mesh, cells] = set_up_mesh(G)
%[mesh, cells] = set_up_mesh(G) Set up vectors containing the indices 
% of the different types of nodes
%
% Output arguments:
%      mesh: object containing information about the entire mesh
%      cells: collection of objects containing information about the 
%             individual cells
%
% Input arguments:
%      G: object containing information about the domain geometry

% Read geometry
Nx = G.Nx; Ny = G.Ny; Nz = G.Nz; N = G.N;

% Set up indices for outer boundary
% Set up indices for outer boundary (corners)
mesh.e_lsw = 1;
mesh.e_lse = Nx;
mesh.e_lnw = (Ny-1)*Nx + 1;
mesh.e_lne = Nx*Ny;
mesh.e_hsw = (Nz-1)*Ny*Nx+1;
mesh.e_hse = (Nz-1)*Ny*Nx+Nx;
mesh.e_hnw = ((Nz-1)*Ny+Ny-1)*Nx + 1;
mesh.e_hne = Nx*Ny*Nz;

% Set up indices for outer boundary (lines)
mesh.e_lw = ((2:Ny-1)-1)*Nx + 1;
mesh.e_le = ((2:Ny-1)-1)*Nx + Nx;
mesh.e_ls = 2:Nx-1;
mesh.e_ln = (Ny-1)*Nx + (2:Nx-1);
mesh.e_hw = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + 1;
mesh.e_he = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + Nx;
mesh.e_hs = (Nz-1)*Ny*Nx+(2:Nx-1);
mesh.e_hn = ((Nz-1)*Ny+Ny-1)*Nx + (2:Nx-1);
mesh.e_sw = (1:Nz-2)*Ny*Nx+1;
mesh.e_se = (1:Nz-2)*Ny*Nx+Nx;
mesh.e_nw = ((1:Nz-2)*Ny+Ny-1)*Nx+1;
mesh.e_ne = ((1:Nz-2)*Ny+Ny-1)*Nx+Nx;

% Set up indices for outer boundary (sides)
[x, y] = meshgrid(2:Nx-1, 2:Ny-1);
mesh.e_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
mesh.e_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, Nz*ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
[x, z] = meshgrid(2:Nx-1, 2:Nz-1);
mesh.e_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
mesh.e_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, Ny*ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
[y, z] = meshgrid(2:Ny-1, 2:Nz-1);
mesh.e_w = sort(reshape(sub2ind([Nx,Ny,Nz], ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';
mesh.e_e =  sort(reshape(sub2ind([Nx,Ny,Nz], Nx*ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';

% Set up indices for the cell
[mesh, cells] = combine_cells(G, mesh);

% Set up indices for the remaining extracellular cells
indices = ones(1, N);
indices([mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse, ...
            mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, ...
            mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, ...
            mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
            mesh.e_l, mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, mesh.m_hsw, mesh.m_hse, ...
            mesh.m_hnw, mesh.m_hne, mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, ...
            mesh.m_lw, mesh.m_le, mesh.m_ls, mesh.m_ln, mesh.m_ne, mesh.m_sw, ...
            mesh.m_se, mesh.m_nw, mesh.m_w, mesh.m_e, mesh.m_s, mesh.m_n, mesh.m_h, ...
            mesh.m_l, mesh.i_all]) = 0;
mesh.e = round(find(indices));

% Set up indices for the membrane potential
mesh.v = sort([mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, mesh.m_hsw, mesh.m_hse, ...
            mesh.m_hnw, mesh.m_hne, mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, ...
            mesh.m_lw, mesh.m_le, mesh.m_ls, mesh.m_ln, mesh.m_ne, mesh.m_sw, ...
            mesh.m_se, mesh.m_nw, mesh.m_w, mesh.m_e, mesh.m_s, mesh.m_n, mesh.m_h, ...
            mesh.m_l]);
        
% Set up indices for the gap junction
mesh.w = sort([mesh.gx_e, mesh.gy_n]);

% Set up internal indices for the membrane potential
for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).v) = 1;
    all = all(mesh.v);
    cells(n).v_idx = find(all);
end

% Set up indices for the nodes not in the intracellular domain
mesh.not_i = [mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, ...
    mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne, mesh.e_hw, ...
    mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le, ...
    mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, mesh.e_se, ...
    mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
    mesh.e_l, mesh.e, mesh.gx_he, mesh.gx_le, mesh.gx_se, ...
    mesh.gx_ne, mesh.gx_lse, mesh.gx_hse, mesh.gx_lne, mesh.gx_hne, ...
    mesh.gy_ne, mesh.gy_nw, mesh.gy_ln, mesh.gy_hn, mesh.gy_lnw, ...
    mesh.gy_hnw, mesh.gy_lne, mesh.gy_hne];


% Set up indices for the extracellular domain
mesh.e_all = sort([mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, mesh.e_hsw, mesh.e_hse, ...
            mesh.e_hnw, mesh.e_hne, mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, ...
            mesh.e_lw, mesh.e_le, mesh.e_ls, mesh.e_ln, mesh.e_ne, mesh.e_sw, ...
            mesh.e_se, mesh.e_nw, mesh.e_w, mesh.e_e, mesh.e_s, mesh.e_n, mesh.e_h, ...
            mesh.e_l, mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, mesh.m_hsw, mesh.m_hse, ...
            mesh.m_hnw, mesh.m_hne, mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, ...
            mesh.m_lw, mesh.m_le, mesh.m_ls, mesh.m_ln, mesh.m_ne, mesh.m_sw, ...
            mesh.m_se, mesh.m_nw, mesh.m_w, mesh.m_e, mesh.m_s, mesh.m_n, mesh.m_h, ...
            mesh.m_l,  mesh.gx_he, mesh.gx_le, mesh.gx_se, mesh.gx_ne, ...
            mesh.gy_hn, mesh.gy_ln, mesh.gy_nw, mesh.gy_ne,  mesh.gx_lse, ...
            mesh.gx_lne, mesh.gx_hse, mesh.gx_hne, mesh.gy_lnw, mesh.gy_lne, ...
            mesh.gy_hnw, mesh.gy_hne, mesh.e]);

% Set up indices for the intracellular domain
mesh.i_all = sort(unique([mesh.v, mesh.w, mesh.i, mesh.i_lsw, mesh.i_lse, mesh.i_lnw, mesh.i_lne, mesh.i_hsw, mesh.i_hse, ...
            mesh.i_hnw, mesh.i_hne, mesh.i_hw, mesh.i_he, mesh.i_hs, mesh.i_hn, ...
            mesh.i_lw, mesh.i_le, mesh.i_ls, mesh.i_ln, mesh.i_ne, mesh.i_sw, ...
            mesh.i_se, mesh.i_nw, mesh.i_w, mesh.i_e, mesh.i_s, mesh.i_n, mesh.i_h, ...
            mesh.i_l, mesh.gx_e, mesh.gy_n, mesh.gxi_w, mesh.gxi_hw, mesh.gxi_lw, ...
            mesh.gxi_sw, mesh.gxi_nw, mesh.gyi_hs, mesh.gyi_ls, mesh.gyi_sw, ...
            mesh.gyi_se, mesh.gyi_s, mesh.gxi_lsw, mesh.gxi_lnw, mesh.gxi_hsw, ...
            mesh.gxi_hnw, mesh.gyi_lsw, mesh.gyi_lse, mesh.gyi_hsw, mesh.gyi_hse]));

mesh.i_inner_all = sort(unique([mesh.i, mesh.i_lsw, mesh.i_lse, mesh.i_lnw, mesh.i_lne, mesh.i_hsw, mesh.i_hse, ...
            mesh.i_hnw, mesh.i_hne, mesh.i_hw, mesh.i_he, mesh.i_hs, mesh.i_hn, ...
            mesh.i_lw, mesh.i_le, mesh.i_ls, mesh.i_ln, mesh.i_ne, mesh.i_sw, ...
            mesh.i_se, mesh.i_nw, mesh.i_w, mesh.i_e, mesh.i_s, mesh.i_n, mesh.i_h, ...
            mesh.i_l]));

        
% Set up internal indices for the membrane
for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).v) = 1;
    all = all(mesh.e_all);
    cells(n).v_of_e = find(all);
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).v) = 1;
    all = all(mesh.i_all);
    cells(n).v_of_i = find(all);
end

% Set up internal indices for the intercalated discs
for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gx_w) = 1;
    all = all(mesh.i_all);
    cells(n).gxw_of_i = find(all);
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gy_s) = 1;
    all = all(mesh.i_all);
    cells(n).gys_of_i = find(all);
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gx_e) = 1;
    all = all(mesh.i_all);
    cells(n).gxe_of_i = find(all);
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gy_n) = 1;
    all = all(mesh.i_all);
    cells(n).gyn_of_i = find(all);
end

% Set up internal indices for the cell
for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).c_all) = 1;
    all = all(mesh.i_all);
    cells(n).c_of_i = find(all);
end

% Set up internal indices for the intercalated disc
for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gy_n) = 1;
    all = all(mesh.w);
    cells(n).gyn_of_w = find(all);
    if size(cells(n).gyn_of_w, 2) == 0
        cells(n).gyn_of_w = zeros(0,1);
    end
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gx_e) = 1;
    all = all(mesh.w);
    cells(n).gxe_of_w = find(all);
    if size(cells(n).gxe_of_w, 2) == 0
        cells(n).gxe_of_w = zeros(0,1);
    end
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gx_w) = 1;
    all = all(mesh.w);
    cells(n).gxw_of_w = find(all);
    if size(cells(n).gxw_of_w, 2) == 0
        cells(n).gxw_of_w = zeros(0,1);
    end
end

for n=1:G.num_cells
    all = zeros(G.N, 1);
    all(cells(n).gy_s) = 1;
    all = all(mesh.w);
    cells(n).gys_of_w = find(all);
    if size(cells(n).gys_of_w, 2) == 0
        cells(n).gys_of_w = zeros(0,1);
    end
end
        
% Inner membrane points
mesh.v_inner = [mesh.m_l, mesh.m_h, mesh.m_s, mesh.m_n, mesh.m_e, mesh.m_w];
all = zeros(G.N, 1);
all(mesh.v_inner) = 1;
all = all(mesh.v);
mesh.v_inner_of_v = find(all);

% Set up extraction of membrane points from the intracellular mesh
all = zeros(G.N, 1);
all(mesh.v) = 1;
all = all(mesh.i_all);
mesh.v_of_i = find(all);

% Set up extraction of membrane points from the extracellular mesh
all = zeros(G.N, 1);
all(mesh.v) = 1;
all = all(mesh.e_all);
mesh.v_of_e = find(all);

% Set up extraction of membrane points from the intracellular mesh
all = zeros(G.N, 1);
all(mesh.v_inner) = 1;
all = all(mesh.i_all);
mesh.v_inner_of_i = find(all);

% Set up extraction of membrane points from the extracellular mesh
all = zeros(G.N, 1);
all(mesh.v_inner) = 1;
all = all(mesh.e_all);
mesh.v_inner_of_e = find(all);

% Set up extraction of intercalated disc points from the intracellular mesh
all = zeros(G.N, 1);
all(mesh.w) = 1;
all = all(mesh.i_all);
mesh.w_of_i = find(all);

end