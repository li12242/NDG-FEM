function GenBCfile( resultFileName, BCfileName, xb )
%MAKEBCFILE Generate boundary condition files.
% Make boundary files from the result nc file.
% Usages:
%   GenBCfile('SWE1D_1000.nc', 'SWE1D_1000_BC2',-510e3);
% 
resultfile = Utilities.PostProcess.ResultFile(resultFileName);
time  = resultfile.GetVarData('time');
ntime = numel(time);
nBV   = numel(xb);
% generate bc file
bcfile = CreateBCFile(BCfileName, ntime, nBV, time, xb);
% make interpolation at each time step
x     = resultfile.GetVarData('x');
for i = 1:ntime
    t  = resultfile.GetTimeVarData('h',i);
    tb = disInterp(x, t, xb);
    bcfile.putVarPart('hb', [0, i-1], [nBV, 1], tb);
    t  = resultfile.GetTimeVarData('q',i);
    tb = disInterp(x, t, xb);
    bcfile.putVarPart('qb', [0, i-1], [nBV, 1], tb);
end% for
bcfile.CloseFile;
end

function yb = disInterp(x, y, xb)
TOTAL = 1e-4;
% the boundary point is on vertex
ind   = abs(x-xb)<TOTAL;
if any( ind )
    yb = mean(y(ind));
    return;
end% if
% find which element xb belong to
dx  = (x(1, :) - xb).*(x(end, :) - xb);
ind = find(dx < 0);
if(isempty(ind))
    error('The boundary point %f is not in the computation domain !\n', xb);
end% if
% interpolation with linear reconstruction
coef1 = (x(end, ind) - xb)/(x(end, ind) - x(1  , ind));
coef2 = (x(1  , ind) - xb)/(x(1  , ind) - x(end, ind));

yb    = y(1, ind)*coef1 + y(end, ind)*coef2;
end% func

function bcfile = CreateBCFile(filename, ntime, nBV, time, xb)
% Define dimensions
tdim = Utilities.NetcdfClass.NcDim('time', ntime); % unlimited dimensions
node = Utilities.NetcdfClass.NcDim('nBV', nBV);

% Define variables
x    = Utilities.NetcdfClass.NcVar('xb',   node, 'double');
tvar = Utilities.NetcdfClass.NcVar('time', tdim, 'double');
h    = Utilities.NetcdfClass.NcVar('hb', [node, tdim], 'double');
q    = Utilities.NetcdfClass.NcVar('qb', [node, tdim], 'double');

% Define nc file
bcfile = Utilities.NetcdfClass.NcFile(filename, [node, tdim],...
    [x, tvar, h, q]);

bcfile.CreateFile;
bcfile.putVarPart('xb',   0, nBV,   xb);
bcfile.putVarPart('time', 0, ntime, time);
end% func

