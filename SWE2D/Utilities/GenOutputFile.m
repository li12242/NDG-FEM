function file = GenOutputFile(filename, phys)
% Define NetCDF output file for SWE2D
% Input:
%   filename    - name of NetCDF output file
%   phys        - phys structure
% Output:
%   file        - netcdf file structure
% 

mesh = phys.mesh;
time = Utilities.Netcdf.NcDim('time', 0); % unlimited dimensions
node = Utilities.Netcdf.NcDim('node', mesh.Shape.nNode);
nele = Utilities.Netcdf.NcDim('nele', mesh.nElement);

x    = Utilities.Netcdf.NcVar('x', [node, nele], 'double');
y    = Utilities.Netcdf.NcVar('y', [node, nele], 'double');
t    = Utilities.Netcdf.NcVar('time', time, 'double');
h    = Utilities.Netcdf.NcVar('h', [node, nele, time], 'double');
qx   = Utilities.Netcdf.NcVar('qx', [node, nele, time], 'double');
qy   = Utilities.Netcdf.NcVar('qy', [node, nele, time], 'double');
bot  = Utilities.Netcdf.NcVar('bot', [node, nele], 'double');

file = Utilities.Netcdf.NcFile(filename, [node, nele, time],...
    [x, y, t, h, qx, qy, bot]);

% initialize output file
file.CreateFile;

% set vertex location value
file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('y',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.y);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);

end% func