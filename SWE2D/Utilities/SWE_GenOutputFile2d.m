function file = SWE_GenOutputFile2d(filename, phys)
% Define NetCDF output file for SWE2D
% Input:
%   filename    - name of NetCDF output file
%   phys        - phys structure
% Output:
%   file        - netcdf file structure
% 

mesh = phys.mesh;
time = Utilities.NetcdfClass.NcDim('time', 0); % unlimited dimensions
node = Utilities.NetcdfClass.NcDim('node', mesh.Shape.nNode);
nele = Utilities.NetcdfClass.NcDim('nele', mesh.nElement);
one  = Utilities.NetcdfClass.NcDim('one', 1);

x    = Utilities.NetcdfClass.NcVar('x', [node, nele], 'double');
y    = Utilities.NetcdfClass.NcVar('y', [node, nele], 'double');
t    = Utilities.NetcdfClass.NcVar('time', time, 'double');
h    = Utilities.NetcdfClass.NcVar('h', [node, nele, time], 'double');
qx   = Utilities.NetcdfClass.NcVar('qx', [node, nele, time], 'double');
qy   = Utilities.NetcdfClass.NcVar('qy', [node, nele, time], 'double');
bot  = Utilities.NetcdfClass.NcVar('bot', [node, nele], 'double');
elapseTime = Utilities.NetcdfClass.NcVar('elapsedTime', one, 'double');

file = Utilities.NetcdfClass.NcFile(filename, [node, nele, time, one],...
    [x, y, t, h, qx, qy, bot, elapseTime]);

% initialize output file
file.CreateFile;

% set vertex location value
file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('y',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.y);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);

% initial value
step = 0;
Np   = mesh.Shape.nNode; % No. of nodes in each element
Ne   = mesh.nElement; % No. of element
file.putVarPart('time', step, 1, 0);
file.putVarPart('h',  [0,0,step], [Np,Ne,1], phys.h);
file.putVarPart('qx', [0,0,step], [Np,Ne,1], phys.qx);
file.putVarPart('qy', [0,0,step], [Np,Ne,1], phys.qy);
end% func