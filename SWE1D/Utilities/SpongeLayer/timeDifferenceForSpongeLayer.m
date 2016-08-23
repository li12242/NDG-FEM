function dt = timeDifferenceForSpongeLayer

% coarse girds
coarsefile = {'SWE1D_WiderMovingMound_1000.nc'};
% fine grids
finefile   = {'SWE1D_LocalMovingMound_SpongeLayer_333.nc'};
eletype    = 'line';
order      = 1;
coarseStru = Utilities.PostProcess.Postprocess(coarsefile, eletype, order);
fineStru   = Utilities.PostProcess.Postprocess(finefile, eletype, order);

% parameters
fileID = 1;
stime  = 12*3600;
hC     = coarseStru.GetVarData('h', stime, fileID);

timeF  = fineStru.NcFile.GetVarData('time');
nt     = numel(timeF);
xF     = fineStru.NcFile.GetVarData('x');
err    = zeros(nt, 1);

% interp from coarse gird to fine grid
hCF    = coarseStru.Interp1d(hC, xF, fileID);

% loop of the fine grid result, find the minimum error and time step
for i = 1:nt
    hF     = fineStru.GetVarData('h', timeF(i), fileID);
    err(i) = Linf(hF, hCF);
    fprintf('Processing %f ...\n', i/nt);
end% for

[~, ind] = min(err);

dt = timeF(ind) - stime;
end% func

%% subroutine function
function err = Linf(y1, y2)
err = max2( abs(y1-y2) );
end% func

function ymax = max2(y)
ymax = max(max(y));
end

function err = L2(y1, y2)
err = sqrt(sum2( (y1-y2).^2 ));
end% func

function ysum = sum2(y1)
ysum = sum(sum(y1));
end