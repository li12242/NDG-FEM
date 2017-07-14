%% MovingMound1d_CompareResult
function MovingMound1d_CompareResult
% coarse girds
coarfile = {'SWE1D_WiderMovingMound_1000.nc'};
% fine grids
finefile   = {'SWE1D_LocalMovingMound_SpongeLayer_333.nc'};
eletype    = 'line';
order      = 1;
coarStru = Utilities.PostProcess.Postprocess(coarfile, eletype, order);
fineStru = Utilities.PostProcess.Postprocess(finefile, eletype, order);

% parameters
fileID = 1;
timeC  = coarStru.NcFile.GetVarData('time');
timeF  = fineStru.NcFile.GetVarData('time');
xC     = coarStru.NcFile.GetVarData('x');
xF     = fineStru.NcFile.GetVarData('x');

% Initial plot
figure('Position', [157   247   820   455]);
% parameters
bMarkerSize = 10; 
sMakerSize  = 3;
gre = [0    0.4980         0];
% plot of water depth
subplot(2,2,1); hold on; grid on;
hC = coarStru.GetVarData('h', timeC(1), fileID);
hF = fineStru.GetVarData('h', timeC(1), fileID);

ph1 = plot(xC(:), hC(:), 'rx', 'MarkerSize', sMakerSize); 
ph2 = plot(xF(:), hF(:), 'b.', 'MarkerSize', bMarkerSize); 
title('$h$', 'Interpreter', 'latex')
% error of water depth
hCF = coarStru.Interp1d(hC, xF, fileID); % interpolate to fine grid
dh  = hCF - hF;
subplot(2,2,2); 
pdh = plot(xF(:), dh(:), '-.',...
    'MarkerSize', bMarkerSize, ...
    'Color', gre); 
grid on;
title('$Error\, of\, h$', 'Interpreter', 'latex')
% plot of water flux
subplot(2,2,3);  hold on; grid on;
qC = coarStru.GetVarData('q', timeC(1), fileID);
qF = fineStru.GetVarData('q', timeC(1), fileID);
pq1 = plot(xC(:), qC(:), 'rx', 'MarkerSize', sMakerSize); 
pq2 = plot(xF(:), qF(:), 'b.', 'MarkerSize', bMarkerSize); 
title('$q$', 'Interpreter', 'latex')
% error of water flux
qCF = coarStru.Interp1d(qC, xF, fileID);
dq  = qCF - qF;
subplot(2,2,4); 
pdq = plot(xF(:), dq(:), '-.', ...
    'MarkerSize', bMarkerSize, ...
    'Color', gre); grid on;
title('$Error\, of\, q$', 'Interpreter', 'latex')
% iteration
% dt = 919.5409;
dt = 0;
for i = 1:20:numel(timeF);
    time1 = timeF(i) - dt;
    if time1<0
        time1 = 0;
    end
    hC = coarStru.GetVarData('h', time1, fileID);
    hF = fineStru.GetVarData('h', timeF(i), fileID);
    set(ph1, 'YData', hC(:)); set(ph2, 'YData', hF(:));
    
    hCF = coarStru.Interp1d(hC, xF, fileID);
    dh  = hCF - hF;
    set(pdh, 'YData', dh(:));
    
    qC = coarStru.GetVarData('q', time1, fileID);
    qF = fineStru.GetVarData('q', timeF(i), fileID);
    set(pq1, 'YData', qC(:)); set(pq2, 'YData', qF(:));
    
    qCF = coarStru.Interp1d(qC, xF, fileID);
    dq  = qCF - qF;
    set(pdq, 'YData', dq(:));
    drawnow;
end% for

end% func