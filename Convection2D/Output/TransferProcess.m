function TransferProcess
filename = 'Convection2D_4_60.nc';

time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x'); nx = numel(x);
y = ncread(filename, 'y');

itime = 1;
var = ncread(filename, 'h', [1, itime],[inf, 1]);

p_h = point3D(filename, x, y, var, time);

end% func

function p_h = point3D(filename, x, y, var, time)
ind = 1:3:numel(x);

p_h = plot3(x(ind), y(ind), var(ind), '.');
zlim([-0.5, 1.2])

for itime = 1:numel(time)
    var = ncread(filename, 'h', [1, itime],[inf, 1]);
    
    set(p_h, 'ZData', var(ind));
    drawnow;

    fprintf('Processing: %f ...\n', itime/numel(time))
end% for
end% func