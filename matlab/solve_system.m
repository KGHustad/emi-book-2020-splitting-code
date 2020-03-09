function [U_emi, V_emi, W_emi] = solve_system(G, mesh, cells)
%[U_emi, V_emi, W_emi] = solve_system(G, mesh, cells) 
% Run a simulation of the EMI model using the splitting procedure
%
% Output arguments:
%      U_emi: intracellular and extracellular solutions
%      V_emi: membrane potential solutions
%      W_emi: gap junction potential solutions

% Set up initial conditions
if isfield(G, 'updated_initial_states')
    states = G.updated_initial_states*ones(1, G.nv);
else
    states = G.init_states()*ones(1, G.nv);
end
W = zeros(G.nw, 1);
t = 0;

% Set up cell model parameters
P = set_up_stim_param(G.init_parameters, G, mesh);

% Adjust the initial states and parameters in a part of the domain
if isfield(G, 'special_initial_states')
    states(:, mesh.special_v) = G.special_initial_states*ones(1, length(mesh.special_v));
    P(G.special_P_idx, mesh.special_v) = G.special_P_value*ones(1, length(mesh.special_v));
end

% Extract the membrane potential
V = states(G.V_idx,:)';

% Set up matrices and rhs-vectors for the intracellular problems
fprintf('Setting up intracellular system. This may also take some time...\n')
Ui = zeros(G.ni, 1);
Ig = zeros(G.ni, 1);
for n=1:G.num_cells
    cells(n).Ai = set_up_cell_matrix(G, cells(n));
    cells(n).di = zeros(length(cells(n).c_all), 1);
    Ui(cells(n).c_of_i) = mean(V(cells(n).v_idx));
end

% Set up matrices and rhs-vectors for the extracellular problem
fprintf('Setting up extracellular system.\n')
Ae = set_up_extracellular_matrix(G, mesh);
Ue = zeros(G.ne, 1);
M = set_up_Im_matrix(G, mesh);

% Set up Dirichlet boundary vector
de = zeros(G.ne, 1);

% Set up parameters for saving the solution
if isfield(G, 'DT')
    G.num_save = round(G.Tstop/G.DT) + 1;
    G.save_step = round(G.DT/G.dt);
else
    G.num_save = G.Nt+1;
    G.save_step = 1;
    G.DT = G.dt;
end

if isfield(G, 'DT_V')
    G.num_save_v = round(G.Tstop/G.DT_V) + 1;
    G.save_step_v = round(G.DT_V/G.dt);
else
    G.num_save_v = G.num_save;
    G.save_step_v = 1;
    G.DT_V = G.DT;
end

% Set up matrices for saving the solution
U_emi = zeros(G.N, G.num_save);
V_emi = zeros(G.nv, G.num_save_v);
V_emi(:,1) = V;
W_emi = zeros(G.nw, G.num_save);
W_emi(:,1) = W;

% Set up gap junction resistance
Rg_factor = G.Rg*ones(G.nw, 1);

% Run simulation
print_step = G.dt;
n_print = max(round(print_step/G.dt), 1);
t1 = tic;
fprintf('Started simulation.\n')
for n = 1:G.Nt

    % Step 1: Solve ode-system for membrane model 
    for k=1:G.nt
        states = states + G.dt_ode*G.rhs(t, states, P);
        t = t + G.dt_ode;
    end 
    V = states(G.V_idx, :)';
    if any(isnan(V)) || any(isinf(V))
        fprintf('ODE solution is inf or nan.\n')
        break;
    end
    
    % Step 2: Intracellular PDEs
    W_prev = W;
    for q=1:G.N_it
        for j=1:G.M_it
            Ig(mesh.w_of_i) = (1./Rg_factor).*W + G.Cg*(W-W_prev)/G.dt;
            for k=1:G.num_cells
                % Update right-hand side vector
                cells(k).di(:) = 0;
                cells(k).di(cells(k).v_of_c) = G.Cm/G.dt*(V(cells(k).v_idx)+Ue(cells(k).v_of_e));
                cells(k).di(cells(k).gxw_of_c) = Ig(cells(k).gxw_of_i);
                cells(k).di(cells(k).gys_of_c) = Ig(cells(k).gys_of_i);
                cells(k).di(cells(k).gxe_of_c) = Ig(cells(k).gxe_of_i);
                cells(k).di(cells(k).gyn_of_c) = Ig(cells(k).gyn_of_i);
                
                % Solve system
                X = cells(k).Ai\cells(k).di;
                
                % Update solution
                W(cells(k).gxw_of_w) = Ui(cells(k).gxw_of_i) - X(cells(k).gxw_of_c);
                W(cells(k).gys_of_w) = Ui(cells(k).gys_of_i) - X(cells(k).gys_of_c);
                Ui(cells(k).c_of_i) = X;
            end
            
        end

        % Step 3: Compute extracellular potential (EMI model)
        
        % Update right-hand side vector
        Im = M*Ui;
        de(mesh.v_of_e) = -Im(mesh.v_of_i);

        % Solve the system
        Ue = Ae\de;

    end
  
    % Step 4: Update V
    V = Ui(mesh.v_of_i) - Ue(mesh.v_of_e);
    
    % Save the solution
    if rem(n, G.save_step) == 0
        U_full = zeros(G.N, 1);
        U_full(mesh.i_all) = Ui;
        U_full(mesh.e_all) = Ue;
        U_emi(:, round(n/G.save_step)+1) = U_full;
        W_emi(:, round(n/G.save_step)+1) = W;
    end
    
    if rem(n, G.save_step_v) == 0
        V_emi(:, round(n/G.save_step_v)+1) = V;
    end

    % Update the membrane potential in the cell model
    states(G.V_idx, :) = V';
    
    % Estimate remaining simulation time
    if rem(n, n_print) == 0
        % Print current point in time
        fprintf('t = %g ms. ', t);

        % Print estimated simulation time
        t2 = toc(t1);                % Time usage for n_print time steps
        t_rem = t2*(G.Nt-n)/n_print; % Estimated remaining simulation time
        fprintf('Estimated remaining simulation time: ');
        if t_rem > 86400*2
            fprintf('%.1f days \n', t_rem/86400);
        elseif t_rem > 3600
            fprintf('%.1f h \n', t_rem/3600);
        elseif t_rem > 60
            fprintf('%.1f min \n', t_rem/60);
        else
            fprintf('%.1f sec \n', t_rem);
        end
        t1 = tic;
        
    end
end


end


