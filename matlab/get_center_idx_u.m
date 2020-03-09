function idx = get_center_idx_u(G)
%idx = get_center_idx_u(G) Get the indices for the center x,y-sheet of the
%domain
%
% Output argument:
%      idx: indices for the center x,y-sheet of the domain
%
% Input argument:
%     G: object containing information about the EMI problem

z_idx = round(G.Nz/2);
U = zeros(G.Nx*G.Ny, G.Nz);
U(:, z_idx) = 1;
U_new = reshape(U, G.Nx*G.Ny*G.Nz, 1);
idx = find(U_new);

end