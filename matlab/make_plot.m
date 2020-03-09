% Plot the solution 

% Set up figure
figure('Units','centimeters', 'Position', [10 10 24 20], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [24, 20]);

% Load solution
load('solution.mat')
[x,y] = meshgrid(1e4*(0:G.dx:G.Lx)', 1e4*(0:G.dy:G.Ly)');

for n=1:length(t_points)-1
    
    % Load and reshape current solution
    Ui = Ui_all(:, n);
    Ue = Ue_all(:, n);
    Ui = reshape(Ui, G.Nx, G.Ny);
    Ue = reshape(Ue, G.Nx, G.Ny);
    
    subplot(length(t_points)-1, 2, (n-1)*2+1)
    pc = pcolor(x,y,(Ue)');
    set(pc, 'edgecolor', 'none')
    xlim([0 G.Lx*1e4]) 
    ylim([0 G.Ly*1e4])
    set(gca, 'fontsize', 14)
    ylabel('y (\mum)')
    pos = get(gca, 'Position');
    if n==1
        title('Extracellular potential', 'fontsize', 20)
        colormap jet
        c = colorbar;
        ylabel(c, 'u_e (mV)')
        set(c, 'Position', [0.45, 0.2, 0.022, 0.6])
    end
    if n==length(t_points)-1
        xlabel('x (\mum)')
    else
        set(gca, 'XTickLabels', {})
    end
    caxis([-0.8, 0.8])
    text(-0.28, 0.5, sprintf('t = %g ms', t_points(n)-G.stim_start_time), 'Units', 'normalized', ...
        'rotation', 90, 'fontsize', 18, 'HorizontalAlignment', 'center', ...
        'fontweight', 'bold')
    grid off
    pos(1) = pos(1) - 0.02;
    pos(3) = pos(3) - 0.02;
    set(gca, 'Position', pos)
    
    subplot(length(t_points)-1, 2, (n-1)*2+2)
    pc = pcolor(x,y,(Ui)');
    set(pc, 'edgecolor', 'none')
    xlim([0 G.Lx*1e4]) 
    ylim([0 G.Ly*1e4])
    set(gca, 'fontsize', 14, 'YTickLabels', {})
    pos = get(gca, 'Position');
    if n==1
        title('Intracellular potential', 'fontsize', 20)
        colormap jet
        c = colorbar;
        ylabel(c, 'u_i (mV)')
        set(c, 'Position', [0.9, 0.2, 0.022, 0.6])
    end
    if n==length(t_points)-1
        xlabel('x (\mum)')
    else
        set(gca, 'XTickLabels', {})
    end
    caxis([-81, 40])
    grid off
    pos(1) = pos(1) - 0.01;
    pos(3) = pos(3) - 0.02;
    set(gca, 'Position', pos)
    
end

print('-r300', '-dpng', 'ischemia.png')
