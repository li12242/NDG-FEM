function Anima_DamBreakDry2d
%ANIMADAMBREAKDRY Summary of this function goes here
%   Detailed explanation goes here

%% Parameter
meshtype = 'quad';
filename = {'SWE2D_DamBreakDry_quad_80.nc'};
order    = 1;
PostproTri  = Utilities.PostProcess.Postprocess(filename, meshtype, order);
fileID   = 1;
time     = PostproTri.NcFile(1).GetVarData('time');

%% Create animationi
videoName = 'DamBreakDry.avi';
writerObj = VideoWriter(videoName);
writerObj.FrameRate=15;
open(writerObj)

figure('Color', 'w');
for i = 1:numel(time)-1
    PostproTri.Snapshot2D('h', time(i), fileID); drawnow;
    zlim([0, 11]);
    view(14, 16);
    zlabel('ˮλ (m)','FontSize', 14);
    xlabel('x (m)','FontSize', 14);
    ylabel('y (m)','FontSize', 14);
    grid on;
    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    
    cla;
end% for

close(writerObj);

end