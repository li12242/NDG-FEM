function file = CreateOutputFile(filename, phys)

mesh = phys.mesh;
%% Define dimensions
time = Utilities.netcdf.NcDim('time', 0); % unlimited dimensions
node = Utilities.netcdf.NcDim('node', mesh.Shape.nNode);
nele = Utilities.netcdf.NcDim('nele', mesh.nElement);

%% Define variables
x    = Utilities.netcdf.NcVar('x', [node, nele], 'double');
t    = Utilities.netcdf.NcVar('time', time, 'double');
h    = Utilities.netcdf.NcVar('h', [node, nele, time], 'double');
q    = Utilities.netcdf.NcVar('q', [node, nele, time], 'double');
bot  = Utilities.netcdf.NcVar('bot', [node, nele], 'double');

%% Define nc file
file = Utilities.netcdf.NcFile(filename, [node, nele, time],...
    [x, t, h, q, bot]);

%% Create file and put coordinate and bottom level
file.CreateFile;

file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);
end% func
