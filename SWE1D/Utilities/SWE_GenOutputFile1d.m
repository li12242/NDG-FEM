function file = SWE_GenOutputFile1d(filename, phys)

mesh = phys.mesh;
%% Define dimensions
time = Utilities.NetcdfClass.NcDim('time', 0); % unlimited dimensions
node = Utilities.NetcdfClass.NcDim('node', mesh.Shape.nNode);
nele = Utilities.NetcdfClass.NcDim('nele', mesh.nElement);

%% Define variables
x    = Utilities.NetcdfClass.NcVar('x', [node, nele], 'double');
t    = Utilities.NetcdfClass.NcVar('time', time, 'double');
h    = Utilities.NetcdfClass.NcVar('h', [node, nele, time], 'double');
q    = Utilities.NetcdfClass.NcVar('q', [node, nele, time], 'double');
bot  = Utilities.NetcdfClass.NcVar('bot', [node, nele], 'double');

%% Define nc file
file = Utilities.NetcdfClass.NcFile(filename, [node, nele, time],...
    [x, t, h, q, bot]);

%% Create file and put coordinate and bottom level
file.CreateFile;

file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);
end% func
