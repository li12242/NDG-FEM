function AnimaDamBreakDry2d
%ANIMADAMBREAKDRY Summary of this function goes here
%   Detailed explanation goes here

%% Parameter
meshtype = 'quad';
filename = {'SWE2D.nc'};
order    = 1;
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = PostproTri.NcFile(1).GetVarData('time');

for i = 1:numel(time)
    PostproTri.Snapshot2D('h', time(i), fileID); drawnow;
    zlim([0, 11]);
    view(14, 16);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    box on;
    cla;
end% for

end