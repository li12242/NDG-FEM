function file = CreateOutputFile(filename, phys)

mesh = phys.mesh;
%% Define dimensions
time = Utilities.Netcdf.NcDim('time', 0); % unlimited dimensions
node = Utilities.Netcdf.NcDim('node', mesh.Shape.nNode);
nele = Utilities.Netcdf.NcDim('nele', mesh.nElement);

%% Define variables
x    = Utilities.Netcdf.NcVar('x', [node, nele], 'double');
t    = Utilities.Netcdf.NcVar('time', time, 'double');
h    = Utilities.Netcdf.NcVar('h', [node, nele, time], 'double');
q    = Utilities.Netcdf.NcVar('q', [node, nele, time], 'double');
bot  = Utilities.Netcdf.NcVar('bot', [node, nele], 'double');

%% Define nc file
file = Utilities.Netcdf.NcFile(filename, [node, nele, time],...
    [x, t, h, q, bot]);

%% Create file and put coordinate and bottom level
file.CreateFile;

file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);
end% func
