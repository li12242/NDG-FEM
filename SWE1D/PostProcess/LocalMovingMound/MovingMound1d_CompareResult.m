%% MovingMound1d_CompareResult
% 
% 
function MovingMound1d_CompareResult
% coarse girds
coarsefile = 'SWE1D_WiderMovingMound_1000.nc';
% fine grids
finefile   = 'SWE1D_LocalMovingMound_SpongeLayer_333.nc';
coarseStru = Utilities.PostProcess.ResultFile(coarsefile);
fineStru   = Utilities.PostProcess.ResultFile(finefile);

timeC = ncread(coarsefile, 'time');
timeF = ncread(finefile, 'time');
xC    = ncread(coarsefile, 'x');
xF    = ncread(finefile, 'x');

% Initial plot
figure('Position', [157   247   820   455]);
% parameters
bMarkerSize = 10; sMakerSize = 3;
gre = [0    0.4980         0];
% plot of water depth
subplot(2,2,1); hold on; grid on;
hC = GetInterTimeStep(coarseStru, 'h', timeC, timeC(1));
hF = fineStru.GetTimeVarData('h', 1);
ph1 = plot(xC(:), hC(:), 'rx', 'MarkerSize', sMakerSize); 
ph2 = plot(xF(:), hF(:), 'b.', 'MarkerSize', bMarkerSize); 
title('$h$', 'Interpreter', 'latex')
% error of water depth
hCF = GetSpatialInterpolatedResult(xC, hC, xF);
dh  = hCF - hF;
subplot(2,2,2); 
pdh = plot(xF(:), dh(:), '.', 'MarkerSize', bMarkerSize, ...
    'Color', gre); grid on;
title('$Error\, of\, h$', 'Interpreter', 'latex')
% plot of water flux
subplot(2,2,3);  hold on; grid on;
qC = GetInterTimeStep(coarseStru, 'q', timeC, timeC(1));
qF = fineStru.GetTimeVarData('q', 1);
pq1 = plot(xC(:), qC(:), 'rx', 'MarkerSize', sMakerSize); 
pq2 = plot(xF(:), qF(:), 'b.', 'MarkerSize', bMarkerSize); 
title('$q$', 'Interpreter', 'latex')
% error of water flux
qCF = GetSpatialInterpolatedResult(xC, qC, xF);
dq  = qCF - qF;
subplot(2,2,4); 
pdq = plot(xF(:), dq(:), 'c.', 'MarkerSize', bMarkerSize, ...
    'Color', gre); grid on;
title('$Error\, of\, q$', 'Interpreter', 'latex')
% iteration
dt = 963.9818;
for i = 1:20:numel(timeF);
    time1 = timeF(i) - dt;
    if time1<0
        time1 = 0;
    end
    hC = GetInterTimeStep(coarseStru, 'h', timeC, time1);
    hF = fineStru.GetTimeVarData('h', i);
    set(ph1, 'YData', hC(:)); set(ph2, 'YData', hF(:));
    
    hCF = GetSpatialInterpolatedResult(xC, hC, xF);
    dh  = hCF - hF;
    set(pdh, 'YData', dh(:));
    
    qC = GetInterTimeStep(coarseStru, 'q', timeC, time1);
    qF = fineStru.GetTimeVarData('q', i);
    set(pq1, 'YData', qC(:)); set(pq2, 'YData', qF(:));
    
    qCF = GetSpatialInterpolatedResult(xC, qC, xF);
    dq  = qCF - qF;
    set(pdq, 'YData', dq(:));
    drawnow;
end% for

end% func

%% GetSpatialInterpolatedResult
function yb = GetSpatialInterpolatedResult(x, y, xb)
yb = zeros(size(xb));
for i = 1:numel(xb)
    try
        yb(i) = disInterp(x, y, xb(i));
    catch
        keyboard
    end
end% for
end% func

%% disInterp
% Interpolation function for single node
function yb = disInterp(x, y, xb)
TOTAL = 1e-2;
% the boundary point is on vertex
ind   = abs(x-xb)<TOTAL;
if any(any( ind ))
    yb = mean(y(ind));
    return;
end% if
% find which element xb belong to
% if (xb<x(1,1))
%     yb = y(1, 1);
%     return
% elseif (xb>x(end, end))
%     yb = y(end, end);
%     return
% end% if
dx  = (x(1, :) - xb).*(x(end, :) - xb);
ind = find(dx < 0);
% interpolation with linear reconstruction
coef1 = (x(end, ind) - xb)/(x(end, ind) - x(1  , ind));
coef2 = (x(1  , ind) - xb)/(x(1  , ind) - x(end, ind));

yb    = y(1, ind)*coef1 + y(end, ind)*coef2;
end% func

%% GetInterTimeStep
% Get the interpolated result at spicific time.
% Input:
%   filename - NC file name
%   varname  - variable name
%   time     - vector of time
%   stime    - spicific time
% Output:
%   h        - linear interpolated result at spicific time
function h = GetInterTimeStep(filePostStruct,varname,time,stime)
% check for input time
if (stime > time(end))
    error('The input time %f is out of computation time!\n', stime);
% elseif (stime < time(1))
%     error('The input time %f is out of computation time!\n', stime);
end% if

sk = find(time>=stime, 1);
if (sk==1)
    skm1 = 1;
    w1 = 1;
    w2 = 0;
else 
    skm1 = sk - 1;
    dt1 = abs(stime - time(skm1));
    dt2 = abs(stime - time(sk));
    w1  = dt2/(dt1+dt2);
    w2  = dt1/(dt1+dt2);
end% if

h1 = filePostStruct.GetTimeVarData(varname, skm1);
h2 = filePostStruct.GetTimeVarData(varname, sk);

h  = h1*w1 + h2*w2;
end% func