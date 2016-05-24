function TransferProcess
filename = 'Convection2D_1_60.nc';

p_h = point3D(filename);
% p_h = contour2D(filename);
end% func

function p_h = point3D(filename)

time = ncread(filename, 'time');
x = ncread(filename, 'x');
y = ncread(filename, 'y');

% plot 3D points
ind = 1:1:numel(x);
itime = 1;
var = ncread(filename, 'var', [1, itime],[inf, 1]);

p_h = plot3(x(ind), y(ind), var(ind), '.');
zlim([-0.5, 1.2]); view(-22.7, 57.2);

camera_on = 1;
if camera_on
    writerObj = VideoWriter('Triangle1_60.avi');
    writerObj.FrameRate = 60;
    open(writerObj);
end

for itime = 1:1:numel(time)
    
    var = ncread(filename, 'var', [1, itime],[inf, 1]);
%     dis = ncread(filename, 'temp', [1, itime], [inf, 1]);
    
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

function p_h = contour2D(filename)

% get grid points
time = ncread(filename, 'time');
x = ncread(filename, 'x');
y = ncread(filename, 'y');

% set interpolation mesh
np = 100;
t = linspace(min(x), max(x), np);
[X, Y] = meshgrid(t, t);

% set levels
v = [-0.1:0.1:1.1];
for itime = 1:numel(time)
    var = ncread(filename, 'var', [1, itime],[inf, 1]);

    % Interpolation
    Interp = TriScatteredInterp(x, y, var, 'linear');
    Z = Interp(X, Y);
    
    contourf(X, Y, Z, v);
    drawnow;
end
end% func