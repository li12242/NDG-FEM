function MovingMound1d_SinglePoint

% coarse girds
coarsefile = {'SWE1D_WiderMovingMound_1000.nc'};
% fine grids
finefile   = {'SWE1D_LocalMovingMound_SpongeLayer_1000.nc'};
eletype    = 'line';
order      = 1;
coarseStru = Utilities.PostProcess.Postprocess(coarsefile, eletype, order);
fineStru   = Utilities.PostProcess.Postprocess(finefile, eletype, order);


timeR = coarseStru.NcFile.GetVarData('time');
xb    = 0;

nt    = 500;
timeC = linspace(timeR(1), timeR(end), nt);
% interpolation result
hbC   = zeros(nt, 1);
hbF   = zeros(nt, 1);
fileID = 1;
dt    = 197.0821;
for i = 1:nt
    time = timeC(i);
    hC   = coarseStru.GetVarData('h', time1, fileID);
    hF   = fineStru.GetVarData('h', time, fileID);
    hbC(i) = coarseStru.Interp1d(hC, xb, fileID);
    hbF(i) = fineStru.Interp1d(hF, xb, fileID);
    fprintf('Interpolating %f...\n', i/nt);
end% for

%% draw figure
% Initial plot
figure('Position', [157   247   820   455]);
subplot(2,2,1); hold on;
plot(timeC, hbF, 'r.');
plot(timeC, hbC, 'b.');

dh = abs(hbF - hbC);
subplot(2,2,2); hold on;
plot(timeC, dh, 'g.');
end% func