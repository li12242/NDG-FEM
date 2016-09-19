function SnapShot_TsunamiRunup

%% Parameter
T        = 22;
meshtype = 'quad';
filename = {'SWE2D.nc'};
order    = 1;
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = linspace(1e-12, T, 4);
for i = 1:numel(time)
    subplot(2,2,i); hold on;
    bh  = PostproTri.SnapshotConst2D('bot',fileID);
    set(bh, 'FaceColor', [0.4, 0.4, .4]);
    bot = get(bh, 'Vertices'); % get topography
    ph  = PostproTri.Snapshot2D('h', time(i), fileID);
    h   = get(ph, 'Vertices');
    ind = h(:, 3) < 1e-2;
    h(ind, 3) = nan;
    h(:, 3)   = h(:, 3) + bot(:, 3);
    set(ph, 'Vertices', h);
%     zlim([0, 5]);
    view(-140, 50);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    title(['time ', num2str(i), 'T/4'])
    box on
end% for
end% func