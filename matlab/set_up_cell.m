function cell = set_up_cell(G, idx)
%cell = set_up_cell(G, idx) Set up mesh for a single 3D cell
%
% Output arguments:
%       cell: object containing information about the mesh of the cell
%
% Input arguments:
%       G: object containing information about the problem
%     idx: current cell index

% Compute local indices
idx_x = rem(idx-1, G.num_cells_x) + 1;
idx_y = floor((idx-1)/G.num_cells_x) + 1;

% Load parameters
dx = G.dx; dy = G.dy; dz = G.dz; Nx = G.Nx; Ny = G.Ny; Nz = G.Nz;
cell_start_x = G.cell_start_x(idx_x); 
cell_start_y = G.cell_start_y(idx_y);
cell_start_z = G.cell_start_z;
cell_length_x = G.cell_length_x;
cell_length_y = G.cell_length_y;
cell_length_z = G.cell_length_z;
j_length_x = G.j_length_x;
j_length_y = G.j_length_y;
j_width_y_x = G.j_width_y_x;
j_width_y_z = G.j_width_y_z;
j_width_x_y = G.j_width_x_y;
j_width_x_z = G.j_width_x_z;


% Define cell indices
cell_start_ix = round(cell_start_x/dx)+1;
cell_start_iy = round(cell_start_y/dy)+1;
cell_start_iz = round(cell_start_z/dz)+1;
cell_end_ix = cell_start_ix + round(cell_length_x/dx);
cell_end_iy = cell_start_iy + round(cell_length_y/dy);
cell_end_iz = cell_start_iz + round(cell_length_z/dz);

% Number of intervals for the junction parts
j_length_x_n = round(j_length_x/dx);
j_length_y_n = round(j_length_y/dy);
j_width_x_n = round(j_width_y_x/dx); 
j_width_y_z_n = round(j_width_y_z/dz);
j_width_y_n = round(j_width_x_y/dy);
j_width_x_z_n = round(j_width_x_z/dz); 

% Number of intervals for the non-junction part
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

% Points for start and end of junction parts
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


% Set up indices for the membrane points
cell.m_lsw = ((cell_start_iz-1)*Ny + (cell_start_iy-1))*Nx + cell_start_ix;
cell.m_lse = ((cell_start_iz-1)*Ny +(cell_start_iy-1))*Nx + cell_end_ix;
cell.m_lnw = ((cell_start_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_start_ix;
cell.m_lne = ((cell_start_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_end_ix;
cell.m_hsw = ((cell_end_iz-1)*Ny + (cell_start_iy-1))*Nx + cell_start_ix;
cell.m_hse = ((cell_end_iz-1)*Ny +(cell_start_iy-1))*Nx + cell_end_ix;
cell.m_hnw = ((cell_end_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_start_ix;
cell.m_hne = ((cell_end_iz-1)*Ny + (cell_end_iy-1))*Nx + cell_end_ix;

% Set up indices for the membrane lines
cell.m_ls = ((cell_start_iz-1)*Ny + (cell_start_iy-1))*Nx + (cell_start_ix+1:cell_end_ix-1);
cell.m_ls = [cell.m_ls, ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + (j_start_ix+1:cell_start_ix-1)];
cell.m_ls = [cell.m_ls, ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + (cell_end_ix+1:j_end_ix-1)];

cell.m_hs = ((cell_end_iz-1)*Ny + (cell_start_iy-1))*Nx + (cell_start_ix+1:cell_end_ix-1);
cell.m_hs = [cell.m_hs, ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + (j_start_ix+1:cell_start_ix-1)];
cell.m_hs = [cell.m_hs, ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + (cell_end_ix+1:j_end_ix-1)];

cell.m_ln = ((cell_start_iz-1)*Ny + (cell_end_iy-1))*Nx + (cell_start_ix+1:cell_end_ix-1);
cell.m_ln = [cell.m_ln, ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + (j_start_ix+1:cell_start_ix-1)];
cell.m_ln = [cell.m_ln, ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + (cell_end_ix+1:j_end_ix-1)];

cell.m_hn = ((cell_end_iz-1)*Ny + (cell_end_iy-1))*Nx + (cell_start_ix+1:cell_end_ix-1);
cell.m_hn = [cell.m_hn, ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + (j_start_ix+1:cell_start_ix-1)];
cell.m_hn = [cell.m_hn, ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + (cell_end_ix+1:j_end_ix-1)];

cell.m_lw = ((cell_start_iz-1)*Ny + ((cell_start_iy+1:cell_end_iy-1)-1))*Nx + cell_start_ix;
cell.m_lw = [cell.m_lw, ((cut_start_iyz-1)*Ny + ((j_start_iy+1:cell_start_iy-1)-1))*Nx + cut_start_ix];
cell.m_lw = [cell.m_lw, ((cut_start_iyz-1)*Ny + ((cell_end_iy+1:j_end_iy-1)-1))*Nx + cut_start_ix];

cell.m_hw = ((cell_end_iz-1)*Ny + ((cell_start_iy+1:cell_end_iy-1)-1))*Nx + cell_start_ix;
cell.m_hw = [cell.m_hw, ((cut_end_iyz-1)*Ny + ((j_start_iy+1:cell_start_iy-1)-1))*Nx + cut_start_ix];
cell.m_hw = [cell.m_hw, ((cut_end_iyz-1)*Ny + ((cell_end_iy+1:j_end_iy-1)-1))*Nx + cut_start_ix];

cell.m_le = ((cell_start_iz-1)*Ny + ((cell_start_iy+1:cell_end_iy-1)-1))*Nx + cell_end_ix;
cell.m_le = [cell.m_le, ((cut_start_iyz-1)*Ny + ((j_start_iy+1:cell_start_iy-1)-1))*Nx + cut_end_ix];
cell.m_le = [cell.m_le, ((cut_start_iyz-1)*Ny + ((cell_end_iy+1:j_end_iy-1)-1))*Nx + cut_end_ix];

cell.m_he = ((cell_end_iz-1)*Ny + ((cell_start_iy+1:cell_end_iy-1)-1))*Nx + cell_end_ix;
cell.m_he = [cell.m_he, ((cut_end_iyz-1)*Ny + ((j_start_iy+1:cell_start_iy-1)-1))*Nx + cut_end_ix];
cell.m_he = [cell.m_he, ((cut_end_iyz-1)*Ny + ((cell_end_iy+1:j_end_iy-1)-1))*Nx + cut_end_ix];

cell.m_sw = (((cell_start_iz+1:cell_end_iz-1)-1)*Ny + (cell_start_iy-1))*Nx + cell_start_ix;
cell.m_se = (((cell_start_iz+1:cell_end_iz-1)-1)*Ny + (cell_start_iy-1))*Nx + cell_end_ix;
cell.m_nw = (((cell_start_iz+1:cell_end_iz-1)-1)*Ny + (cell_end_iy-1))*Nx + cell_start_ix;
cell.m_ne = (((cell_start_iz+1:cell_end_iz-1)-1)*Ny + (cell_end_iy-1))*Nx + cell_end_ix;

% Set up indices for the membrane planes
% Five upper and lower planes
[x, y] = meshgrid(cell_start_ix+1:cell_end_ix-1, cell_start_iy+1:cell_end_iy-1);
cell.m_l = reshape(sub2ind([Nx,Ny,Nz], x, y, cell_start_iz*ones(size(x))), size(x,1)*size(x,2), 1)';
cell.m_h = reshape(sub2ind([Nx,Ny,Nz], x, y, cell_end_iz*ones(size(x))), size(x,1)*size(x,2), 1)';
[x, y] = meshgrid(j_start_ix+1:cell_start_ix-1, cut_start_iy+1:cut_end_iy-1);
cell.m_l = [cell.m_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_start_ixz*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.m_h = [cell.m_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_end_ixz*ones(size(x))), size(x,1)*size(x,2), 1)'];
[x, y] = meshgrid(cell_end_ix+1:j_end_ix-1, cut_start_iy+1:cut_end_iy-1);
cell.m_l = [cell.m_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_start_ixz*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.m_h = [cell.m_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_end_ixz*ones(size(x))), size(x,1)*size(x,2), 1)'];
[x, y] = meshgrid(cut_start_ix+1:cut_end_ix-1, j_start_iy+1:cell_start_iy-1);
cell.m_l = [cell.m_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_start_iyz*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.m_h = [cell.m_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_end_iyz*ones(size(x))), size(x,1)*size(x,2), 1)'];
[x, y] = meshgrid(cut_start_ix+1:cut_end_ix-1, cell_end_iy+1:j_end_iy-1);
cell.m_l = [cell.m_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_start_iyz*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.m_h = [cell.m_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, cut_end_iyz*ones(size(x))), size(x,1)*size(x,2), 1)'];

% Membrane planes on junctions planes
[x, z] = meshgrid(j_start_ix+1:cell_start_ix-1, cut_start_ixz+1:cut_end_ixz-1);
cell.m_n = reshape(sub2ind([Nx,Ny,Nz], x, cut_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
cell.m_s = reshape(sub2ind([Nx,Ny,Nz], x, cut_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_end_ix+1:j_end_ix-1, cut_start_ixz+1:cut_end_ixz-1);
cell.m_n = [cell.m_n, ...
    reshape(sub2ind([Nx,Ny,Nz], x, cut_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)'];
cell.m_s = [cell.m_s, ...
    reshape(sub2ind([Nx,Ny,Nz], x, cut_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)'];
[y, z] = meshgrid(j_start_iy+1:cell_start_iy-1, cut_start_iyz+1:cut_end_iyz-1);
cell.m_w = reshape(sub2ind([Nx,Ny,Nz], cut_start_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
cell.m_e = reshape(sub2ind([Nx,Ny,Nz], cut_end_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_end_iy+1:j_end_iy-1, cut_start_iyz+1:cut_end_iyz-1);
cell.m_w = [cell.m_w, ...
    reshape(sub2ind([Nx,Ny,Nz], cut_start_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)'];
cell.m_e = [cell.m_e, ...
    reshape(sub2ind([Nx,Ny,Nz], cut_end_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)'];

% Membrane planes on cell
[x, z] = meshgrid(cut_start_ix:cut_end_ix, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], x, cell_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_start_ix+1:cell_end_ix-1, cell_start_iz+1:cell_end_iz-1);
entire = reshape(sub2ind([Nx,Ny,Nz], x, cell_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.m_n = [cell.m_n, find(all_ind)'];

[x, z] = meshgrid(cut_start_ix:cut_end_ix, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], x, cell_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_start_ix+1:cell_end_ix-1, cell_start_iz+1:cell_end_iz-1);
entire = reshape(sub2ind([Nx,Ny,Nz], x, cell_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.m_s = [cell.m_s, find(all_ind)'];

[y, z] = meshgrid(cut_start_iy:cut_end_iy, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], cell_start_ix*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_start_iy+1:cell_end_iy-1, cell_start_iz+1:cell_end_iz-1);
entire = reshape(sub2ind([Nx,Ny,Nz], cell_start_ix*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.m_w = [cell.m_w, find(all_ind)'];

[y, z] = meshgrid(cut_start_iy:cut_end_iy, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], cell_end_ix*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_start_iy+1:cell_end_iy-1, cell_start_iz+1:cell_end_iz-1);
entire = reshape(sub2ind([Nx,Ny,Nz], cell_end_ix*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.m_e = [cell.m_e, find(all_ind)'];


% Set up indices for the membrane points between cells
cell.gx_lsw = ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + j_start_ix;
cell.gx_lnw = ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + j_start_ix;
cell.gx_lse = ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + j_end_ix;
cell.gx_lne = ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + j_end_ix;
cell.gx_hsw = ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + j_start_ix;
cell.gx_hnw = ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + j_start_ix;
cell.gx_hse = ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + j_end_ix;
cell.gx_hne = ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + j_end_ix;

cell.gy_lsw = ((cut_start_iyz-1)*Ny + (j_start_iy-1))*Nx + cut_start_ix;
cell.gy_lse = ((cut_start_iyz-1)*Ny + (j_start_iy-1))*Nx + cut_end_ix;
cell.gy_lnw = ((cut_start_iyz-1)*Ny + (j_end_iy-1))*Nx + cut_start_ix;
cell.gy_lne = ((cut_start_iyz-1)*Ny + (j_end_iy-1))*Nx + cut_end_ix;
cell.gy_hsw = ((cut_end_iyz-1)*Ny + (j_start_iy-1))*Nx + cut_start_ix;
cell.gy_hse = ((cut_end_iyz-1)*Ny + (j_start_iy-1))*Nx + cut_end_ix;
cell.gy_hnw = ((cut_end_iyz-1)*Ny + (j_end_iy-1))*Nx + cut_start_ix;
cell.gy_hne = ((cut_end_iyz-1)*Ny + (j_end_iy-1))*Nx + cut_end_ix;

% Set up indices for membrane lines between cells
cell.gx_lw = ((cut_start_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + j_start_ix;
cell.gx_hw = ((cut_end_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + j_start_ix;
cell.gx_sw = (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + cut_start_iy-1)*Nx + j_start_ix;
cell.gx_nw = (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + cut_end_iy-1)*Nx + j_start_ix;

cell.gx_le = ((cut_start_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + j_end_ix;
cell.gx_he = ((cut_end_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + j_end_ix;
cell.gx_se = (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + cut_start_iy-1)*Nx + j_end_ix;
cell.gx_ne = (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + cut_end_iy-1)*Nx + j_end_ix;

cell.gy_ls = ((cut_start_iyz-1)*Ny + (j_start_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1);
cell.gy_hs = ((cut_end_iyz-1)*Ny + (j_start_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1);
cell.gy_sw = (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + j_start_iy-1)*Nx + cut_start_ix;
cell.gy_se = (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + j_start_iy-1)*Nx + cut_end_ix;

cell.gy_ln = ((cut_start_iyz-1)*Ny + (j_end_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1);
cell.gy_hn = ((cut_end_iyz-1)*Ny + (j_end_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1);
cell.gy_nw = (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + j_end_iy-1)*Nx + cut_start_ix;
cell.gy_ne = (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + j_end_iy-1)*Nx + cut_end_ix;

% Membrane planes between cells
[x, z] = meshgrid(cut_start_ix+1:cut_end_ix-1, cut_start_iyz+1:cut_end_iyz-1);
cell.gy_s = reshape(sub2ind([Nx,Ny,Nz], x, j_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
cell.gy_n = reshape(sub2ind([Nx,Ny,Nz], x, j_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[y, z] = meshgrid(cut_start_iy+1:cut_end_iy-1, cut_start_iyz+1:cut_end_iyz-1);
cell.gx_w = reshape(sub2ind([Nx,Ny,Nz], j_start_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
cell.gx_e = reshape(sub2ind([Nx,Ny,Nz], j_end_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';


     %%%%%%%%%%%%%% Inner boundary of the cell %%%%%%%%%%%%%%

% Set up indices for the inner points
cell.i_lsw = ((cell_start_iz+1-1)*Ny + (cell_start_iy+1-1))*Nx + cell_start_ix+1;
cell.i_lse = ((cell_start_iz+1-1)*Ny +(cell_start_iy+1-1))*Nx + cell_end_ix-1;
cell.i_lnw = ((cell_start_iz+1-1)*Ny + (cell_end_iy-1-1))*Nx + cell_start_ix+1;
cell.i_lne = ((cell_start_iz+1-1)*Ny + (cell_end_iy-1-1))*Nx + cell_end_ix-1;
cell.i_hsw = ((cell_end_iz-1-1)*Ny + (cell_start_iy-1+1))*Nx + cell_start_ix+1;
cell.i_hse = ((cell_end_iz-1-1)*Ny +(cell_start_iy-1+1))*Nx + cell_end_ix-1;
cell.i_hnw = ((cell_end_iz-1-1)*Ny + (cell_end_iy-1-1))*Nx + cell_start_ix+1;
cell.i_hne = ((cell_end_iz-1-1)*Ny + (cell_end_iy-1-1))*Nx + cell_end_ix-1;

% Extra inner points between junction and cell
cell.i_lsw = [cell.i_lsw, ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + cell_start_ix];
cell.i_lsw = [cell.i_lsw, ((cut_start_iyz-1)*Ny + (cell_start_iy-1))*Nx + cut_start_ix];
cell.i_hsw = [cell.i_hsw, ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + cell_start_ix];
cell.i_hsw = [cell.i_hsw, ((cut_end_iyz-1)*Ny + (cell_start_iy-1))*Nx + cut_start_ix];
cell.i_lnw = [cell.i_lnw, ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + cell_start_ix];
cell.i_lnw = [cell.i_lnw, ((cut_start_iyz-1)*Ny + (cell_end_iy-1))*Nx + cut_start_ix];
cell.i_hnw = [cell.i_hnw, ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + cell_start_ix];
cell.i_hnw = [cell.i_hnw, ((cut_end_iyz-1)*Ny + (cell_end_iy-1))*Nx + cut_start_ix];
cell.i_lse = [cell.i_lse, ((cut_start_ixz-1)*Ny + (cut_start_iy-1))*Nx + cell_end_ix];
cell.i_lse = [cell.i_lse, ((cut_start_iyz-1)*Ny + (cell_start_iy-1))*Nx + cut_end_ix];
cell.i_hse = [cell.i_hse, ((cut_end_ixz-1)*Ny + (cut_start_iy-1))*Nx + cell_end_ix];
cell.i_hse = [cell.i_hse, ((cut_end_iyz-1)*Ny + (cell_start_iy-1))*Nx + cut_end_ix];
cell.i_lne = [cell.i_lne, ((cut_start_ixz-1)*Ny + (cut_end_iy-1))*Nx + cell_end_ix];
cell.i_lne = [cell.i_lne, ((cut_start_iyz-1)*Ny + (cell_end_iy-1))*Nx + cut_end_ix];
cell.i_hne = [cell.i_hne, ((cut_end_ixz-1)*Ny + (cut_end_iy-1))*Nx + cell_end_ix];
cell.i_hne = [cell.i_hne, ((cut_end_iyz-1)*Ny + (cell_end_iy-1))*Nx + cut_end_ix];

% Set up indices for the inner lines
cell.i_ls = ((cell_start_iz-1+1)*Ny + (cell_start_iy-1+1))*Nx + (cell_start_ix+1+1:cell_end_ix-1-1);
cell.i_ls = [cell.i_ls, ((cut_start_ixz-1+1)*Ny + (cut_start_iy-1+1))*Nx + (j_start_ix+1+1:cell_start_ix-1)];
cell.i_ls = [cell.i_ls, ((cut_start_ixz-1+1)*Ny + (cut_start_iy-1+1))*Nx + (cell_end_ix+1:j_end_ix-1-1)];

cell.i_hs = ((cell_end_iz-1-1)*Ny + (cell_start_iy-1+1))*Nx + (cell_start_ix+1+1:cell_end_ix-1-1);
cell.i_hs = [cell.i_hs, ((cut_end_ixz-1-1)*Ny + (cut_start_iy-1+1))*Nx + (j_start_ix+1+1:cell_start_ix-1)];
cell.i_hs = [cell.i_hs, ((cut_end_ixz-1-1)*Ny + (cut_start_iy-1+1))*Nx + (cell_end_ix+1:j_end_ix-1-1)];

cell.i_ln = ((cell_start_iz-1+1)*Ny + (cell_end_iy-1-1))*Nx + (cell_start_ix+1+1:cell_end_ix-1-1);
cell.i_ln = [cell.i_ln, ((cut_start_ixz-1+1)*Ny + (cut_end_iy-1-1))*Nx + (j_start_ix+1+1:cell_start_ix-1)];
cell.i_ln = [cell.i_ln, ((cut_start_ixz-1+1)*Ny + (cut_end_iy-1-1))*Nx + (cell_end_ix+1:j_end_ix-1-1)];

cell.i_hn = ((cell_end_iz-1-1)*Ny + (cell_end_iy-1-1))*Nx + (cell_start_ix+1+1:cell_end_ix-1-1);
cell.i_hn = [cell.i_hn, ((cut_end_ixz-1-1)*Ny + (cut_end_iy-1-1))*Nx + (j_start_ix+1+1:cell_start_ix-1)];
cell.i_hn = [cell.i_hn, ((cut_end_ixz-1-1)*Ny + (cut_end_iy-1-1))*Nx + (cell_end_ix+1:j_end_ix-1-1)];

cell.i_lw = ((cell_start_iz-1+1)*Ny + ((cell_start_iy+1+1:cell_end_iy-1-1)-1))*Nx + cell_start_ix+1;
cell.i_lw = [cell.i_lw, ((cut_start_iyz-1+1)*Ny + ((j_start_iy+1+1:cell_start_iy-1)-1))*Nx + cut_start_ix+1];
cell.i_lw = [cell.i_lw, ((cut_start_iyz-1+1)*Ny + ((cell_end_iy+1:j_end_iy-1-1)-1))*Nx + cut_start_ix+1];

cell.i_hw = ((cell_end_iz-1-1)*Ny + ((cell_start_iy+1+1:cell_end_iy-1-1)-1))*Nx + cell_start_ix+1;
cell.i_hw = [cell.i_hw, ((cut_end_iyz-1-1)*Ny + ((j_start_iy+1+1:cell_start_iy-1)-1))*Nx + cut_start_ix+1];
cell.i_hw = [cell.i_hw, ((cut_end_iyz-1-1)*Ny + ((cell_end_iy+1:j_end_iy-1-1)-1))*Nx + cut_start_ix+1];

cell.i_le = ((cell_start_iz-1+1)*Ny + ((cell_start_iy+1+1:cell_end_iy-1-1)-1))*Nx + cell_end_ix-1;
cell.i_le = [cell.i_le, ((cut_start_iyz-1+1)*Ny + ((j_start_iy+1+1:cell_start_iy-1)-1))*Nx + cut_end_ix-1];
cell.i_le = [cell.i_le, ((cut_start_iyz-1+1)*Ny + ((cell_end_iy+1:j_end_iy-1-1)-1))*Nx + cut_end_ix-1];

cell.i_he = ((cell_end_iz-1-1)*Ny + ((cell_start_iy+1+1:cell_end_iy-1-1)-1))*Nx + cell_end_ix-1;
cell.i_he = [cell.i_he, ((cut_end_iyz-1-1)*Ny + ((j_start_iy+1+1:cell_start_iy-1)-1))*Nx + cut_end_ix-1];
cell.i_he = [cell.i_he, ((cut_end_iyz-1-1)*Ny + ((cell_end_iy+1:j_end_iy-1-1)-1))*Nx + cut_end_ix-1];

cell.i_sw = (((cell_start_iz+1+1:cell_end_iz-1-1)-1)*Ny + (cell_start_iy-1+1))*Nx + cell_start_ix+1;
cell.i_se = (((cell_start_iz+1+1:cell_end_iz-1-1)-1)*Ny + (cell_start_iy-1+1))*Nx + cell_end_ix-1;
cell.i_nw = (((cell_start_iz+1+1:cell_end_iz-1-1)-1)*Ny + (cell_end_iy-1-1))*Nx + cell_start_ix+1;
cell.i_ne = (((cell_start_iz+1+1:cell_end_iz-1-1)-1)*Ny + (cell_end_iy-1-1))*Nx + cell_end_ix-1;

% Extra inner lines between junction and cell
cell.i_sw = [cell.i_sw, (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + (cut_start_iy-1))*Nx + cell_start_ix];
cell.i_sw = [cell.i_sw, (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + (cell_start_iy-1))*Nx + cut_start_ix];
cell.i_nw = [cell.i_nw, (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + (cut_end_iy-1))*Nx + cell_start_ix];
cell.i_nw = [cell.i_nw, (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + (cell_end_iy-1))*Nx + cut_start_ix];
cell.i_se = [cell.i_se, (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + (cut_start_iy-1))*Nx + cell_end_ix];
cell.i_se = [cell.i_se, (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + (cell_start_iy-1))*Nx + cut_end_ix];
cell.i_ne = [cell.i_ne, (((cut_start_ixz+1:cut_end_ixz-1)-1)*Ny + (cut_end_iy-1))*Nx + cell_end_ix];
cell.i_ne = [cell.i_ne, (((cut_start_iyz+1:cut_end_iyz-1)-1)*Ny + (cell_end_iy-1))*Nx + cut_end_ix];

cell.i_lw = [cell.i_lw, ((cut_start_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + cell_start_ix];
cell.i_hw = [cell.i_hw, ((cut_end_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + cell_start_ix];
cell.i_le = [cell.i_le, ((cut_start_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + cell_end_ix];
cell.i_he = [cell.i_he, ((cut_end_ixz-1)*Ny + ((cut_start_iy+1:cut_end_iy-1)-1))*Nx + cell_end_ix];

cell.i_ls = [cell.i_ls, ((cut_start_iyz-1)*Ny + (cell_start_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1)];
cell.i_hs = [cell.i_hs, ((cut_end_iyz-1)*Ny + (cell_start_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1)];
cell.i_ln = [cell.i_ln, ((cut_start_iyz-1)*Ny + (cell_end_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1)];
cell.i_hn = [cell.i_hn, ((cut_end_iyz-1)*Ny + (cell_end_iy-1))*Nx + (cut_start_ix+1:cut_end_ix-1)];

% Set up indices for the inner planes
% Five upper and lower planes
% Main cell
[x, y] = meshgrid(cell_start_ix+1+1:cell_end_ix-1-1, cell_start_iy+1+1:cell_end_iy-1-1);
cell.i_l = reshape(sub2ind([Nx,Ny,Nz], x, y, (cell_start_iz+1)*ones(size(x))), size(x,1)*size(x,2), 1)';
cell.i_h = reshape(sub2ind([Nx,Ny,Nz], x, y, (cell_end_iz-1)*ones(size(x))), size(x,1)*size(x,2), 1)';
% Left junction
[x, y] = meshgrid(j_start_ix+1+1:cell_start_ix-1, cut_start_iy+1+1:cut_end_iy-1-1);
cell.i_l = [cell.i_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_start_ixz+1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.i_h = [cell.i_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_end_ixz-1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
% Right junction
[x, y] = meshgrid(cell_end_ix+1:j_end_ix-1-1, cut_start_iy+1+1:cut_end_iy-1-1);
cell.i_l = [cell.i_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_start_ixz+1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.i_h = [cell.i_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_end_ixz-1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
% South junction
[x, y] = meshgrid(cut_start_ix+1+1:cut_end_ix-1-1, j_start_iy+1+1:cell_start_iy-1);
cell.i_l = [cell.i_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_start_iyz+1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.i_h = [cell.i_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_end_iyz-1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
% North junction
[x, y] = meshgrid(cut_start_ix+1+1:cut_end_ix-1-1, cell_end_iy+1:j_end_iy-1-1);
cell.i_l = [cell.i_l, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_start_iyz+1)*ones(size(x))), size(x,1)*size(x,2), 1)'];
cell.i_h = [cell.i_h, ...
    reshape(sub2ind([Nx,Ny,Nz], x, y, (cut_end_iyz-1)*ones(size(x))), size(x,1)*size(x,2), 1)'];

% Membrane planes on junctions planes
[x, z] = meshgrid(j_start_ix+1+1:cell_start_ix-1, cut_start_ixz+1+1:cut_end_ixz-1-1);
cell.i_n = reshape(sub2ind([Nx,Ny,Nz], x, (cut_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
cell.i_s = reshape(sub2ind([Nx,Ny,Nz], x, (cut_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_end_ix+1:j_end_ix-1-1, cut_start_ixz+1+1:cut_end_ixz-1-1);
cell.i_n = [cell.i_n, ...
    reshape(sub2ind([Nx,Ny,Nz], x, (cut_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1)'];
cell.i_s = [cell.i_s, ...
    reshape(sub2ind([Nx,Ny,Nz], x, (cut_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1)'];
[y, z] = meshgrid(j_start_iy+1+1:cell_start_iy-1, cut_start_iyz+1+1:cut_end_iyz-1-1);
cell.i_w = reshape(sub2ind([Nx,Ny,Nz], (cut_start_ix+1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
cell.i_e = reshape(sub2ind([Nx,Ny,Nz], (cut_end_ix-1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_end_iy+1:j_end_iy-1-1, cut_start_iyz+1+1:cut_end_iyz-1-1);
cell.i_w = [cell.i_w, ...
    reshape(sub2ind([Nx,Ny,Nz], (cut_start_ix+1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)'];
cell.i_e = [cell.i_e, ...
    reshape(sub2ind([Nx,Ny,Nz], (cut_end_ix-1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)'];

% Membrane planes on cell
[x, z] = meshgrid(cut_start_ix:cut_end_ix, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], x, (cell_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_start_ix+1+1:cell_end_ix-1-1, cell_start_iz+1+1:cell_end_iz-1-1);
entire = reshape(sub2ind([Nx,Ny,Nz], x, (cell_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.i_n = [cell.i_n, find(all_ind)'];

[x, z] = meshgrid(cut_start_ix:cut_end_ix, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], x, (cell_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[x, z] = meshgrid(cell_start_ix+1+1:cell_end_ix-1-1, cell_start_iz+1+1:cell_end_iz-1-1);
entire = reshape(sub2ind([Nx,Ny,Nz], x, (cell_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.i_s = [cell.i_s, find(all_ind)'];

[y, z] = meshgrid(cut_start_iy:cut_end_iy, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], (cell_start_ix+1)*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_start_iy+1+1:cell_end_iy-1-1, cell_start_iz+1+1:cell_end_iz-1-1);
entire = reshape(sub2ind([Nx,Ny,Nz], (cell_start_ix+1)*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.i_w = [cell.i_w, find(all_ind)'];

[y, z] = meshgrid(cut_start_iy:cut_end_iy, cut_start_iyz:cut_end_iyz);
to_remove = reshape(sub2ind([Nx,Ny,Nz], (cell_end_ix-1)*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
[y, z] = meshgrid(cell_start_iy+1+1:cell_end_iy-1-1, cell_start_iz+1+1:cell_end_iz-1-1);
entire = reshape(sub2ind([Nx,Ny,Nz], (cell_end_ix-1)*ones(size(y)), y,  z), size(y,1)*size(y,2), 1)';
all_ind = zeros(G.N, 1);
all_ind(entire) = 1;
all_ind(to_remove) = 0;
cell.i_e = [cell.i_e, find(all_ind)'];

% Set up indices for the membrane points between cells
cell.gxi_lsw = ((cut_start_ixz-1+1)*Ny + (cut_start_iy-1+1))*Nx + j_start_ix+1;
cell.gxi_lnw = ((cut_start_ixz-1+1)*Ny + (cut_end_iy-1-1))*Nx + j_start_ix+1;
cell.gxi_lse = ((cut_start_ixz-1+1)*Ny + (cut_start_iy-1+1))*Nx + j_end_ix-1;
cell.gxi_lne = ((cut_start_ixz-1+1)*Ny + (cut_end_iy-1-1))*Nx + j_end_ix-1;
cell.gxi_hsw = ((cut_end_ixz-1-1)*Ny + (cut_start_iy-1+1))*Nx + j_start_ix+1;
cell.gxi_hnw = ((cut_end_ixz-1-1)*Ny + (cut_end_iy-1-1))*Nx + j_start_ix+1;
cell.gxi_hse = ((cut_end_ixz-1-1)*Ny + (cut_start_iy-1+1))*Nx + j_end_ix-1;
cell.gxi_hne = ((cut_end_ixz-1-1)*Ny + (cut_end_iy-1-1))*Nx + j_end_ix-1;

cell.gyi_lsw = ((cut_start_iyz-1+1)*Ny + (j_start_iy-1+1))*Nx + cut_start_ix+1;
cell.gyi_lse = ((cut_start_iyz-1+1)*Ny + (j_start_iy-1+1))*Nx + cut_end_ix-1;
cell.gyi_lnw = ((cut_start_iyz-1+1)*Ny + (j_end_iy-1-1))*Nx + cut_start_ix+1;
cell.gyi_lne = ((cut_start_iyz-1+1)*Ny + (j_end_iy-1-1))*Nx + cut_end_ix-1;
cell.gyi_hsw = ((cut_end_iyz-1-1)*Ny + (j_start_iy-1+1))*Nx + cut_start_ix+1;
cell.gyi_hse = ((cut_end_iyz-1-1)*Ny + (j_start_iy-1+1))*Nx + cut_end_ix-1;
cell.gyi_hnw = ((cut_end_iyz-1-1)*Ny + (j_end_iy-1-1))*Nx + cut_start_ix+1;
cell.gyi_hne = ((cut_end_iyz-1-1)*Ny + (j_end_iy-1-1))*Nx + cut_end_ix-1;

% Set up indices for membrane lines between cells
cell.gxi_lw = ((cut_start_ixz-1+1)*Ny + ((cut_start_iy+1+1:cut_end_iy-1-1)-1))*Nx + j_start_ix+1;
cell.gxi_hw = ((cut_end_ixz-1-1)*Ny + ((cut_start_iy+1+1:cut_end_iy-1-1)-1))*Nx + j_start_ix+1;
cell.gxi_sw = (((cut_start_ixz+1+1:cut_end_ixz-1-1)-1)*Ny + cut_start_iy-1+1)*Nx + j_start_ix+1;
cell.gxi_nw = (((cut_start_ixz+1+1:cut_end_ixz-1-1)-1)*Ny + cut_end_iy-1-1)*Nx + j_start_ix+1;

cell.gxi_le = ((cut_start_ixz-1+1)*Ny + ((cut_start_iy+1+1:cut_end_iy-1-1)-1))*Nx + j_end_ix-1;
cell.gxi_he = ((cut_end_ixz-1-1)*Ny + ((cut_start_iy+1+1:cut_end_iy-1-1)-1))*Nx + j_end_ix-1;
cell.gxi_se = (((cut_start_ixz+1+1:cut_end_ixz-1-1)-1)*Ny + cut_start_iy-1+1)*Nx + j_end_ix-1;
cell.gxi_ne = (((cut_start_ixz+1+1:cut_end_ixz-1-1)-1)*Ny + cut_end_iy-1-1)*Nx + j_end_ix-1;

cell.gyi_ls = ((cut_start_iyz-1+1)*Ny + (j_start_iy-1+1))*Nx + (cut_start_ix+1+1:cut_end_ix-1-1);
cell.gyi_hs = ((cut_end_iyz-1-1)*Ny + (j_start_iy-1+1))*Nx + (cut_start_ix+1+1:cut_end_ix-1-1);
cell.gyi_sw = (((cut_start_iyz+1+1:cut_end_iyz-1-1)-1)*Ny + j_start_iy-1+1)*Nx + cut_start_ix+1;
cell.gyi_se = (((cut_start_iyz+1+1:cut_end_iyz-1-1)-1)*Ny + j_start_iy-1+1)*Nx + cut_end_ix-1;

cell.gyi_ln = ((cut_start_iyz-1+1)*Ny + (j_end_iy-1-1))*Nx + (cut_start_ix+1+1:cut_end_ix-1-1);
cell.gyi_hn = ((cut_end_iyz-1-1)*Ny + (j_end_iy-1-1))*Nx + (cut_start_ix+1+1:cut_end_ix-1-1);
cell.gyi_nw = (((cut_start_iyz+1+1:cut_end_iyz-1-1)-1)*Ny + j_end_iy-1-1)*Nx + cut_start_ix+1;
cell.gyi_ne = (((cut_start_iyz+1+1:cut_end_iyz-1-1)-1)*Ny + j_end_iy-1-1)*Nx + cut_end_ix-1;

% Membrane planes between cells
[x, z] = meshgrid(cut_start_ix+1+1:cut_end_ix-1-1, cut_start_iyz+1+1:cut_end_iyz-1-1);
cell.gyi_s = reshape(sub2ind([Nx,Ny,Nz], x, (j_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
cell.gyi_n = reshape(sub2ind([Nx,Ny,Nz], x, (j_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1)';
[y, z] = meshgrid(cut_start_iy+1+1:cut_end_iy-1-1, cut_start_iyz+1+1:cut_end_iyz-1-1);
cell.gxi_w = reshape(sub2ind([Nx,Ny,Nz], (j_start_ix+1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';
cell.gxi_e = reshape(sub2ind([Nx,Ny,Nz], (j_end_ix-1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1)';

% Set up indices for all inner points
[x, y, z] = meshgrid(cell_start_ix:cell_end_ix, cell_start_iy:cell_end_iy, cell_start_iz:cell_end_iz);
inner_cell = reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1)';
[x, y, z] = meshgrid(j_start_ix:cell_start_ix, cut_start_iy:cut_end_iy, cut_start_ixz:cut_end_ixz);
west_j = reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1)';
[x, y, z] = meshgrid(cell_end_ix:j_end_ix, cut_start_iy:cut_end_iy, cut_start_ixz:cut_end_ixz);
east_j = reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1)';
[x, y, z] = meshgrid(cut_start_ix:cut_end_ix, j_start_iy:cell_start_iy, cut_start_iyz:cut_end_iyz);
south_j = reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1)';
[x, y, z] = meshgrid(cut_start_ix:cut_end_ix, cell_end_iy:j_end_iy, cut_start_iyz:cut_end_iyz);
north_j = reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1)';

all_i = zeros(G.N, 1);
all_i([inner_cell, west_j, east_j, south_j, north_j]) = 1;
cell.i_all = find(all_i)';
all_i([cell.m_lsw, cell.m_lse, cell.m_lnw, cell.m_lne, cell.m_hsw, cell.m_hse, ...
    cell.m_hnw, cell.m_hne, cell.m_hw, cell.m_he, cell.m_hs, cell.m_hn, ...
    cell.m_lw, cell.m_le, cell.m_ls, cell.m_ln, cell.m_ne, cell.m_sw, ...
    cell.m_se, cell.m_nw, cell.m_w, cell.m_e, cell.m_s, cell.m_n, cell.m_h, ...
    cell.m_l, cell.i_lsw, cell.i_lse, cell.i_lnw, cell.i_lne, cell.i_hsw, ...
    cell.i_hse, cell.i_hnw, cell.i_hne, cell.i_hw, cell.i_he, cell.i_hs, ...
    cell.i_hn, cell.i_lw, cell.i_le, cell.i_ls, cell.i_ln, cell.i_ne, ...
    cell.i_sw, cell.i_se, cell.i_nw, cell.i_w, cell.i_e, cell.i_s, cell.i_n, ...
    cell.i_h, cell.i_l, cell.gx_w, cell.gx_e, cell.gx_hw, cell.gx_lw, cell.gx_sw, ...
    cell.gx_nw, cell.gx_he, cell.gx_le, cell.gx_se, cell.gx_ne, cell.gy_hs, ...
    cell.gy_ls, cell.gy_sw, cell.gy_se, cell.gy_hn, cell.gy_ln, cell.gy_nw, ...
    cell.gy_ne, cell.gy_n, cell.gy_s, cell.gx_lsw, cell.gx_lse, cell.gx_lnw, ...
    cell.gx_lne, cell.gx_hsw, cell.gx_hse, cell.gx_hnw, cell.gx_hne, cell.gy_lsw, ...
    cell.gy_lse, cell.gy_lnw, cell.gy_lne, cell.gy_hsw, cell.gy_hse, cell.gy_hnw, ...
    cell.gy_hne, cell.gxi_w, cell.gxi_e, cell.gxi_hw, cell.gxi_lw, cell.gxi_sw, ...
    cell.gxi_nw, cell.gxi_he, cell.gxi_le, cell.gxi_se, cell.gxi_ne, cell.gyi_hs, ...
    cell.gyi_ls, cell.gyi_sw, cell.gyi_se, cell.gyi_hn, cell.gyi_ln, cell.gyi_nw, ...
    cell.gyi_ne, cell.gyi_n, cell.gyi_s, cell.gxi_lsw, cell.gxi_lse, cell.gxi_lnw, ...
    cell.gxi_lne, cell.gxi_hsw, cell.gxi_hse, cell.gxi_hnw, cell.gxi_hne, ...
    cell.gyi_lsw, cell.gyi_lse, cell.gyi_lnw, cell.gyi_lne, cell.gyi_hsw, ...
    cell.gyi_hse, cell.gyi_hnw, cell.gyi_hne]) = 0;
cell.i = find(all_i)';

end

