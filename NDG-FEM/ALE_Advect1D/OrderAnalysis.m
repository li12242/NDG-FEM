function OrderAnalysis
[p_ori, E_ori] = L2_error('Advec1D');
[p_ale, E_ale] = L2_error('ALE_Advec1D');
% [p_app, E_app] = L2_error('ALE_Approxi_Advec1D');

fprintf('The order of Advec1D: %f\n', p_ori);
fprintf('The order of ALE_Advec1D: %f\n', p_ale);
% fprintf('The order of ALE_Approxi_Advec1D: %f\n', p_app);


x = [1, 0.5, 0.5^2, 0.5^3];
figure; hold on
plot(log10(x), log10(E_ori), 'o-', ...
    log10(x), log10(E_ale), 's-')%,...
%     log10(x), log10(E_app), '^-');
% axis equal;
% legend('fix grid', 'ALE', 'Approximate ALE', 'Location', 'NorthWest')
legend('fix grid', 'ALE', 'Location', 'NorthWest')
end% func

function [p, E] = L2_error(ncFileName)

% h space grid
file1_Name = [ncFileName, '.nc'];
ncid = netcdf.open(file1_Name);
[~, nx] = netcdf.inqDim(ncid, 0);
netcdf.close(ncid);

time = ncread(file1_Name, 'time');
timeNum = numel(time);

x = ncread(file1_Name, 'x', [1, timeNum], [inf, 1], [1, 1]);
u1 = ncread(file1_Name, 'u', [1, timeNum], [inf, 1], [1, 1]);
u1_exat = sin(x - time(end));
E1 = sqrt( sum((u1 - u1_exat).^2) /nx);

% h/2 space grid
file2_Name = [ncFileName, '_2.nc'];
ncid = netcdf.open(file2_Name);
[~, nx] = netcdf.inqDim(ncid, 0);
netcdf.close(ncid);

time = ncread(file2_Name, 'time');
timeNum = numel(time);

x = ncread(file2_Name, 'x', [1, timeNum], [inf, 1], [1, 1]);
u2 = ncread(file2_Name, 'u', [1, timeNum], [inf, 1], [1, 1]);
u2_exat = sin(x - time(end));
E2 = sqrt( sum((u2 - u2_exat).^2) /nx);

% h/4 space grid
file1_Name = [ncFileName, '_3.nc'];
ncid = netcdf.open(file1_Name);
[~, nx] = netcdf.inqDim(ncid, 0);
netcdf.close(ncid);

time = ncread(file1_Name, 'time');
timeNum = numel(time);

x = ncread(file1_Name, 'x', [1, timeNum], [inf, 1], [1, 1]);
u1 = ncread(file1_Name, 'u', [1, timeNum], [inf, 1], [1, 1]);
u1_exat = sin(x - time(end));
E3 = sqrt( sum((u1 - u1_exat).^2) /nx);

% h/8 space grid
file1_Name = [ncFileName, '_4.nc'];
ncid = netcdf.open(file1_Name);
[~, nx] = netcdf.inqDim(ncid, 0);
netcdf.close(ncid);

time = ncread(file1_Name, 'time');
timeNum = numel(time);

x = ncread(file1_Name, 'x', [1, timeNum], [inf, 1], [1, 1]);
u1 = ncread(file1_Name, 'u', [1, timeNum], [inf, 1], [1, 1]);
u1_exat = sin(x - time(end));
E4 = sqrt( sum((u1 - u1_exat).^2) /nx);

% cal accuracy
E = [E1, E2, E3, E4];
% x = [1, 0.5, 0.5^2, 0.5^3];
% plot(log10(x), log10(E), 'ro-');
p = log10(E1/E2)/log10(2);

end% func