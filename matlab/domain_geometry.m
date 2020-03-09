function G = domain_geometry(num_cells_x, num_cells_y, num_stim_x, ...
    num_stim_y, stim_start_x, stim_start_y, dx, dy, dz)
%G = domain_geometry(num_cells_x, num_cells_y, num_stim_x, num_stim_y, ...
% stim_start_x, stim_start_y, dx, dy, dz)
% Generate the geometry of a 3D domain with an number of connected cells
%
% Output arguments:
%     G: object containing information about the EMI problem to be solved
%
% Input arguments:
%     num_cells_x:  Number of cells in the x-direction
%     num_cells_y:  Number of cells in the y-direction
%     num_stim_x:   Number of cells to stimulate in the x-direction
%     num_stim_y:   Number of cells to stimulate in the y-direction
%     stim_start_x: First cell to stimulate in the x-direction
%     stim_start_y: First cell to stimulate in the y-direction
%     dx, dy, dz:   Discretization parameters

% Set up discretization parameters
G.dx = dx;
G.dy = dy;
G.dz = dz;
G.num_cells_x = num_cells_x;
G.num_cells_y = num_cells_y;
G.num_cells = G.num_cells_x*G.num_cells_y;

% Set up the cell geometry
G.cell_length_x = 100e-4;
G.cell_length_y = 18e-4;
G.cell_length_z = 18e-4;
G.j_length_x = 4e-4;   % length of Omega_w and Omega_e
G.j_length_y = 4e-4;   % length of Omega_s and Omega_n
G.j_width_y_x = 10e-4; % width (in x-direction) of gap junction in y-direction
G.j_width_y_z = 10e-4; % width (in z-direction) of gap junction in y-direction
G.j_width_x_y = 10e-4; % width (in y-direction) of gap junction in x-direction
G.j_width_x_z = 10e-4; % width (in z-direction) of gap junction in x-direction


% Place the cells in space
distance_to_boundary = 10e-4; % mimimal distance between the intracellular domain 
                             % and the boundary of the extracellular space
distance_to_boundary_z = 4e-4;
G.cell_start_x = distance_to_boundary + G.j_length_x + (0:num_cells_x-1)*(G.cell_length_x+2*G.j_length_x);
G.cell_start_y = distance_to_boundary + G.j_length_y + (0:num_cells_y-1)*(G.cell_length_y+2*G.j_length_y);
G.cell_start_z = distance_to_boundary_z;

% Calculate the domain size
G.Lx = G.cell_start_x(end) + G.cell_length_x + G.j_length_x + distance_to_boundary;
G.Ly = G.cell_start_y(end) + G.cell_length_y + G.j_length_y + distance_to_boundary;
G.Lz = G.cell_start_z(end) + G.cell_length_z + distance_to_boundary_z;

% Compute number of nodes
G.Nx = round(G.Lx/dx)+1;
G.Ny = round(G.Ly/dy)+1;
G.Nz = round(G.Lz/dz)+1;
G.N = G.Nx*G.Ny*G.Nz;

% Fraction of cell length to stimulate for each cell
p_stim = zeros(G.num_cells_x, G.num_cells_y);
p_stim(stim_start_x+(0:num_stim_x-1), stim_start_y+(0:num_stim_y-1)) = 1;  % Stimulate entire cells
p_stim = p_stim(1:G.num_cells_x, 1:G.num_cells_y);
p_stim = reshape(p_stim, G.num_cells, 1);

% Fraction of cell lenght to start stimulating
s_stim = zeros(1, G.num_cells);  % Start stimulating at the beginning of each cell

% Compute stimulation location
G.num_stim = round(p_stim.*(G.cell_length_x+2*G.j_length_x)/G.dx);
G.stim_start = round(s_stim.*(G.cell_length_x+2*G.j_length_x)/G.dx);

% Set up conductivity
G.sigma_i = 4;   % mS/cm 
G.sigma_e = 20;  % mS/cm 

% Define gap junction model (passive)
G.junc_model = 'passive';
G.Cg = 0.5;       % uF/cm^2
G.Rg = 0.0045;    % kOhm cm^2
G.V_gap = 0;      % mV

% Define ionic model (Grandi)
G.ion_model = 'Grandi';
G.Cm = 1;
G.V_idx = 39;
G.init_states = @grandi_init_states;
G.init_parameters = @grandi_init_parameters;
G.rhs = @grandi_rhs;
G.rhs_single = @grandi_rhs_single;
load('grandi_init.mat', 'states')
G.updated_initial_states = states;
load('grandi_init_Ko10.mat', 'states')
G.special_initial_states = states;
G.special_P_value = 10;
G.special_P_idx = 95; % Extracellular potassium concentration

% Define boundary condition
% Dirichlet in x-direction, Neumann in y-, and z-directions
G.bc = 'Neumann_zero_yz'; 

end
