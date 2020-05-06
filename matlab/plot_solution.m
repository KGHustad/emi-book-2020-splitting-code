function plot_solution(U, V, G, names)
%plot_solution(U, V, G, names) 
% Plot the intracellular, extracellular and membrane potential

figure('Units','centimeters', 'Position', [10 10 25 25], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [25, 25]);

idx_z = round(G.Nz/2);


% Plot the intracellular potential
V_reshape = zeros(G.N, 1);
V_reshape(names.v) = V;

% Remove extracellular points
U_tmp = U;
U_tmp([names.e_lsw, names.e_lse, names.e_lnw, names.e_lne, names.e_hsw, names.e_hse, ...
        names.e_hnw, names.e_hne, names.e_hw, names.e_he, names.e_hs, names.e_hn, ...
        names.e_lw, names.e_le, names.e_ls, names.e_ln, names.e_ne, names.e_sw, ...
        names.e_se, names.e_nw, names.e_w, names.e_e, names.e_s, names.e_n, names.e_h, ...
        names.e_l, names.e, names.gx_he, names.gx_le, names.gx_se, ...
        names.gx_ne, names.gx_lse, names.gx_hse, names.gx_lne, names.gx_hne, ...
        names.gy_ne, names.gy_nw, names.gy_ln, names.gy_hn, names.gy_lnw, ...
        names.gy_hnw, names.gy_lne, names.gy_hne]) = nan;

% Calculate Ui on the membrane
U_tmp([names.v]) = U([names.v]) + V_reshape([names.v]);

subplot(3,1,1)
u_tmp = reshape(U_tmp, G.Nx*G.Ny, G.Nz);
u = reshape(u_tmp(:,idx_z), G.Nx, G.Ny);
[x,y] = meshgrid(1e4*(0:G.dx:G.Lx)', 1e4*(0:G.dy:G.Ly)');
pc = pcolor(x,y,(u)');
set(pc, 'edgecolor', 'none')
xlim([0 G.Lx*1e4]) 
ylim([0 G.Ly*1e4])
set(gca, 'fontsize', 14)
ylabel('y (\mum)')
colormap jet
c = colorbar;
ylabel(c, 'mV')
grid off
view(2)
title('Intracellular potential')

% Plot the extracellular potential
U_tmp = U;
U_tmp([names.i_all]) = nan;

subplot(3,1,2)
u_tmp = reshape(U_tmp, G.Nx*G.Ny, G.Nz);
u = reshape(u_tmp(:,idx_z), G.Nx, G.Ny);
[x,y] = meshgrid(1e4*(0:G.dx:G.Lx)', 1e4*(0:G.dy:G.Ly)');
pc = pcolor(x,y,(u)');
set(pc, 'edgecolor', 'none')
xlim([0 G.Lx*1e4]) 
ylim([0 G.Ly*1e4])
set(gca, 'fontsize', 14)
ylabel('y (\mum)')
title('Extracellular potential')
colormap jet
c = colorbar;
ylabel(c, 'mV')
grid off
view(2)

% Plot the membrane potential
U_tmp = nan*ones(G.N,1);
U_tmp(names.v) = V;

subplot(3,1,3)
u_tmp = reshape(U_tmp, G.Nx*G.Ny, G.Nz);
u = reshape(u_tmp(:,idx_z), G.Nx, G.Ny);
[x,y] = meshgrid(1e4*(0:G.dx:G.Lx)', 1e4*(0:G.dy:G.Ly)');
pc = pcolor(x,y,(u)');
set(pc, 'edgecolor', 'none')
xlim([0 G.Lx*1e4]) 
ylim([0 G.Ly*1e4])
set(gca, 'fontsize', 14)
xlabel('x (\mum)')
ylabel('y (\mum)')
colormap jet
c = colorbar;
ylabel(c, 'mV')
grid off
view(2)
title('Membrane potential')


end

