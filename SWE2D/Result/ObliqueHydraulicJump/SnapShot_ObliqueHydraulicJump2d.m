function SnapShot_ObliqueHydraulicJump2d

%% Parameter
T        = 15;
meshtype = 'tri';
filename = {'SWE2D_ObliqueHydraulicJump_tri_50.nc'};
order    = 1;
Postpro  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = [(1e-12:0.5:1)*T, T];

for i = 1:numel(time)
    figure('color', 'w');
    cmap = colormap('winter');
    Postpro.Snapshot2D('h', time(i), fileID, ...
        'value', cmap(end,:), cmap(1, :), ...
        [1, 1.5], 'EdgeColor', 'k',... % 其他参数，包括边颜色，透明度等
        'FaceAlpha', 0.8);
    xlim([0, 40]);
    ylim([0,30]);
    view(16, 52);
    zlabel('$h \rm{(m)}$', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$x \rm{(m)}$', 'Interpreter', 'latex','FontSize', 14);
    ylabel('$y \rm{(m)}$', 'Interpreter', 'latex','FontSize', 14);
    box on;
    grid on;
    colorbar('northoutside');
end% for
end% func