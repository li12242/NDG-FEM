function FlowOver3Bumps_SnapShot

%% Parameter
T        = 16;
meshtype = 'quad';
filename = {'SWE2D.nc'};
order    = 1;
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = (0:0.1:1)*T;
for i = 1:numel(time)
    figure; hold on;
    bh = PostproTri.SnapshotConst2D('bot',fileID);
    bot = get(bh, 'Verti')
    ph = PostproTri.Snapshot2D('h', time(i), fileID);
    set(bh, 'FaceColor', [0.4, 0.4, .4]);
    zlim([0, 5]);
    view(18, 24);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    box on
end% for
end% func