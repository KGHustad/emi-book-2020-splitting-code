function P = set_up_stim_param(func, G, mesh)
%P = set_up_stim_param(func, G, mesh) Set up parameters for a EMI model 
% with stimulation of the correct grid points
%
% Output arguments:
%      P: matrix for the AP model parameters
%
% Input arguments:
%      func: function form defining the AP model parameters
%      G: object containing information about the problem
%      mesh: object containing information about the mesh

[param, param_mesh] = func();
[~, stim_idx] = ismember('stim_amplitude', param_mesh);
stim_amplitude = param(stim_idx);
P = param*ones(1, G.nv);
p_tmp = zeros(G.N, 1);
p_tmp(mesh.to_stim) = stim_amplitude*G.stim_amp;
p_tmp = p_tmp(mesh.v)';
P(stim_idx, :) = p_tmp;

% Stim duration
if isfield(G,  'stim_duration')
   [~, stim_time_idx] = ismember('stim_duration', param_mesh); 
   P(stim_time_idx, :) = G.stim_duration; 
end

% Stim start
if isfield(G,  'stim_start_time')
   [~, stim_time_idx] = ismember('stim_start', param_mesh); 
   P(stim_time_idx, :) = G.stim_start_time; 
end
end

