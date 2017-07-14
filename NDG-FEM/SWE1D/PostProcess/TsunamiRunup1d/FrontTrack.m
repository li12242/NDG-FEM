function [time, rp, rm] = FrontTrack
% get the wet/dry front position

hm = 1e-1;
hp = hm;

% output file
filename = 'SWE1DTsunamiRunup.nc';

% rebuild mesh
x = ncread(filename, 'x'); x = reshape(x, 2, numel(x)/2);

% get results
time = ncread(filename, 'time');
rp = zeros(size(time));
rm = zeros(size(time));
xw = zeros(size(time));

parfor itime = 1:numel(time)
    h = ncread(filename, 'h', [1, itime],[inf, 1]); h = reshape(h, 2, numel(x)/2);
    xw(itime) = getFrontX(x, h, hm);
    rp(itime) = getFrontP(x, h, hp);
    rm(itime) = getFrontM(x, h, hm);
    
    
%     fprintf('Processing: %f ...\n', itime/numel(time))
end% for

% load exact interface position
load('Shoreline.mat');
figure
plot(pt, px, 'k.'); hold on;
% plot(time, rp, 'b.', time, rm, 'b.', time, xw, 'r-')
% t = legend('Exact',  '$r^+$', '$r^-$', 'HREF' );
plot(time, xw, 'r-', 'LineWidth', 1)
t = legend('Exact', 'HREF' );
xlabel('t', 'Interpreter', 'Latex');
ylabel('x', 'Interpreter', 'Latex');
set(t, 'Box', 'off',  'Interpreter', 'Latex');
end% func

function xw = getFrontX(x, h, hm)

wetId = find( h > hm, 1 ) - 1;
xw = x(wetId) + (hm - h(wetId))./(h(wetId + 1) - h(wetId))*(x(wetId + 1) - x(wetId));

end% func

function rp = getFrontP(x, h, hDry)
index = h > hDry;
rp = min(x(index));
end

function rm = getFrontM(x, h, hDry)
index = h <= hDry;
rm = max(x(index));
end% func
