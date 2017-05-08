function TransferProcess
filename = 'Convection1D_4_8.nc';

time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'loc');

p_h = point1D(filename, x, time);
end% func

function p_h = point1D(filename, x, time)
% ind = 1:6:numel(x);
itime = 1;
var = ncread(filename, 'var', [1, 1, itime],[inf, inf, 1]);

p_h = plot(x(:), var(:), '.-');
ylim([-0.5, 1.2]);

camera_on = 0;
if camera_on
    writerObj = VideoWriter('Triangle2_40.avi');
    writerObj.FrameRate = 60;
    open(writerObj);
end

for itime = 1:1:numel(time)
    
    var = ncread(filename, 'var', [1, 1, itime],[inf, inf 1]);
    
    set(p_h, 'YData', var(:));
    drawnow;
    
    if camera_on
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    fprintf('Processing: %f ...\n', itime/numel(time))
end% for

if camera_on
    close(writerObj);
end
end% func