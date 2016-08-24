function GenBCfile( resultFileName, BCfileName, x_spl )
%MAKEBCFILE Generate boundary condition files.
% Make boundary files from the result nc file.
% Usages:
%   GenBCfile('SWE1D_1000.nc', 'SWE1D_1000_BC2',-510e3);
% 

% parameters
eletype    = 'line';
N          = 1;
fileID     = 1;
resultfile = Utilities.PostProcess.Postprocess({resultFileName}, eletype, N);
time       = resultfile.NcFile(fileID).GetVarData('time');
ntime      = numel(time);
[np, ne]   = size(x_spl);
% generate bc file
bcfile = CreateBCFile(BCfileName, ntime, np, ne, time, x_spl);
% make interpolation at each time step
fprintf(['Generation of boundary file ', BCfileName, '\n']);
for i = 1:ntime
    t  = resultfile.NcFile(fileID).GetTimeVarData('h',i);
    tb = resultfile.Interp1d(t, x_spl, fileID);
    bcfile.putVarPart('h_spl', [0, 0, i-1], [np, ne, 1], tb);
    t  = resultfile.NcFile(fileID).GetTimeVarData('q',i);
    tb = resultfile.Interp1d(t, x_spl, fileID);
    bcfile.putVarPart('q_spl', [0, 0, i-1], [np, ne, 1], tb);
end% for
bcfile.CloseFile;
end

function bcfile = CreateBCFile(filename, ntime, np, ne, time, x_spl)
% Define dimensions
tdim = Utilities.NetcdfClass.NcDim('time', ntime); % unlimited dimensions
node = Utilities.NetcdfClass.NcDim('np', np);
nele = Utilities.NetcdfClass.NcDim('ne', ne);

% Define variables
x    = Utilities.NetcdfClass.NcVar('x_spl',[node, nele], 'double');
tvar = Utilities.NetcdfClass.NcVar('time', tdim, 'double');
h    = Utilities.NetcdfClass.NcVar('h_spl', [node, nele, tdim], 'double');
q    = Utilities.NetcdfClass.NcVar('q_spl', [node, nele, tdim], 'double');

% Define nc file
bcfile = Utilities.NetcdfClass.NcFile(filename, [node, nele, tdim],...
    [x, tvar, h, q]);

bcfile.CreateFile;
bcfile.putVarPart('x_spl',[0,0], [np, ne], x_spl);
bcfile.putVarPart('time', 0, ntime, time);
end% func

