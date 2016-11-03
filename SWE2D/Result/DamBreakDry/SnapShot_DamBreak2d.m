function SnapShot_DamBreak2d

%% Parameter
T        = 20;
meshtype = 'tri';
filename = {'SWE2D_DamBreakWet_VBVA_tri_300.nc'};
order    = 1;
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = (0:0.1:1)*T;

for i = 1:numel(time)
    figure
    PostproTri.Snapshot2D('h', time(i), fileID);
    zlim([0, 11]);
    view(30, 32);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    box on
end% for
end% func