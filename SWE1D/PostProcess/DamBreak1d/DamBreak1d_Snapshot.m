function DamBreak1d_Snapshot
%DAMBREAK1D_SNAPSHOT Snapshots of one dimensional DamBreak problem
%   The numerical solution is compared with the exact solution at various
%   time.

% parameter
T = 20;
meshType = 'line';
filename = {'SWE1D_DamBreakWet_400.nc'};
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, 1);
fileID   = 1;
time     = (0:0.2:1)*T;

% draw results
for i = 1:numel(time)
    figure
    PostproTri.Snapshot1D('h', time(i), fileID);
    zlim([0, 11]);
    view(30, 32);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    box on
end% for
end% func