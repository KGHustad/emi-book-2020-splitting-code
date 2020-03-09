% Run a simulation with an area of ischemia in the center of the domain

% Set up parameters
Tstop = 20; % Simulation time
num_cells_x = 22; % Number of cells in the x-direction
num_cells_y = 12; % Number of cells in the y-direction

% Set up discretization parameters
dt = 0.02;        % Time step (ms) 
dx = 2e-4;        % cm
dy = 2e-4;        % cm
dz = 2e-4;        % cm

% Set up stimulation parameters
num_stim_x = 2;   % Number of cells to stimulate in the x-direction
num_stim_y = num_cells_y;   % Number of cells to stimulate in the y-direction
stim_start_x = 1; % Index of the first cell to stimulate in the x-direction
stim_start_y = 1; % Index of the first cell to stimulate in the y-direction

% Generate the domain geometry from the geometry parameters
G = domain_geometry(num_cells_x, num_cells_y, num_stim_x, num_stim_y, ...
    stim_start_x, stim_start_y, dx, dy, dz);

% Specify stimulation
G.stim_start_time = 5; % Time point for stimulation
G.stim_amp = 2;        % Factor with which the original stimulus amplitude is multipled

% Specify time parameters
G.dt = dt; G.dt_ode = 0.001; G.Tstop = Tstop;
G.nt = round(G.dt/G.dt_ode); G.Nt = round(G.Tstop/G.dt);

% Set up mesh
fprintf('Setting up mesh. This may take some time...\n')
[mesh, cells] = set_up_mesh(G);
G.nv = length(mesh.v); G.nw = length(mesh.w);
G.ni = length(mesh.i_all); G.ne = length(mesh.e_all);

% Specify number of splitting iterations
G.N_it = 2; % Outer interations between the intracellular and extracellular domain
G.M_it = 2; % Inner iterations between cells

% Set up ischemic area
ischemic_cells_y = [4,5,6,7,8,9];
ischemic_cells_x = [8,9,10,11,12,13,14,15];
mesh.special_v = [];
for j=1:length(ischemic_cells_y)
    for i=1:length(ischemic_cells_x)
        x_idx = ischemic_cells_x(i);
        y_idx = ischemic_cells_y(j);
        cell_idx = (y_idx-1)*G.num_cells_x + x_idx;
        mesh.special_v = [mesh.special_v; cells(cell_idx).v_idx];
    end
end

% Solve the EMI system using the splitting method
[U, V, W] = solve_system(G, mesh, cells);

% Save the solution for some selected time points
t_points = [7,9,11,13,15];
idx_t = round(t_points/G.dt) + 1;
U = U(:, idx_t);
V = V(:, idx_t);

% Extract the intracellular and extracellular solutions in a sheet in the
% center of the domain
[Ui_all, Ue_all] = extract_Ui_and_Ue(U, V, G, mesh);

save('solution.mat', 'G', 'mesh', 'Ui_all', 'Ue_all', 't_points')
