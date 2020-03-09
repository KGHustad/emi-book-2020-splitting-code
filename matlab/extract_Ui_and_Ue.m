function [Ui, Ue] = extract_Ui_and_Ue(U, V, G, mesh)
%[Ui, Ue] = extract_Ui_and_Ue(U, V, G, mesh)
% Extract the center solutions in the intracellular and extracellular
% domains
%
% Output arguments:
%     Ui, Ue: Intracellular and extracellular solutions in the center of
%             the domain in the z-direction
%
% Input arguments:
%     U: Intracellular and extracellular solutions from the simulation
%     V: Membrane potential solutions from the simulation
%     G: object containing information about the EMI problem
%     mesh: object containing information about the mesh

idx_u = get_center_idx_u(G);

% Set up extracellular part of the solution
Ue = nan*U;
Ue(mesh.e_all, :) = U(mesh.e_all,:);

% Set up intracellular part of the solution
Ui = U;
Ui(mesh.v,:) = Ui(mesh.v,:) + V;
Ui(mesh.not_i,:) = nan;

% Extract only the center of the domain
Ue = Ue(idx_u,:);
Ui = Ui(idx_u,:);

end

