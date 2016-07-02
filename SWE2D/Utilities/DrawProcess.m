function DrawProcess
close all;
DrawSurface;

end% func

function DrawSection
filename = 'SWE2D_ParabolicBowl.nc';
varname  = 'qx';
time     = ncread(filename, 'time');
x        = ncread(filename, 'x');
y        = ncread(filename, 'y');
bot      = ncread(filename, 'bot');

[np, ne] = size(x);
ntime    = numel(time);

%% Read result and interpolation
figure('Position', [114, 375, 660, 429]);

ist = 1;
var = ncread(filename, varname, [1,1,ist], [np, ne, 1]);
p   = plot3(x(:), y(:), bot(:) + var(:), '.');

for ist = 1:ntime
    var = ncread(filename, varname, [1,1,ist], [np, ne, 1]);
    set(p, 'ZData', bot(:) + var(:));
    drawnow; 
    fprintf('Processing %f...\n', ist/ntime);
end
end% func

function DrawSurface
%% Parameters
filename = 'SWE2D.nc';
varname  = 'qx';
time     = ncread(filename, 'time');
x        = ncread(filename, 'x');
y        = ncread(filename, 'y');
bot      = ncread(filename, 'bot');

[np, ne] = size(x);
ntime    = numel(time);

%% Read result and draw pics
figure('Position', [430, 375, 660, 429]);

ist = 1;
var = ncread(filename, varname, [1,1,ist], [np, ne, 1]);
p   = plot3(x(:), y(:), bot(:) + var(:), '.');

for ist = 1:1:ntime
    var = ncread(filename, varname, [1,1,ist], [np, ne, 1]);
    set(p, 'ZData', bot(:) + var(:));
    drawnow;
    fprintf('Processing %f...\n', ist/ntime);
end
end% func