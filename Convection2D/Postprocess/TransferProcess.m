function TransferProcess
filename = 'Convection2D_2_40.nc';

time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x');
y = ncread(filename, 'y');



p_h = point3D(filename, x, y, time);
end% func

function p_h = point3D(filename, x, y, time)
ind = 1:1:numel(x);
itime = 1;
var = ncread(filename, 'h', [1, itime],[inf, 1]);
p_h = plot3(x(ind), y(ind), var(ind), '.');
zlim([-0.5, 1.2]);
view(-22.7, 57.2);

camera_on = 1;
if camera_on
    writerObj = VideoWriter('Triangle2_40.avi');
    writerObj.FrameRate = 60;
    open(writerObj);
end

for itime = 1:1:numel(time)
    var = ncread(filename, 'h', [1, itime],[inf, 1]);
    
    set(p_h, 'ZData', var(ind));
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