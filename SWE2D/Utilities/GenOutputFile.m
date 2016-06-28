function file = GenOutputFile(filename, phys)
% Define NetCDF output file for SWE2D
% Input:
%   filename    - name of NetCDF output file
%   phys        - phys structure
% Output:
%   file        - netcdf file structure
% 

mesh = phys.mesh;
time = Utilities.netcdf.NcDim('time', 0); % unlimited dimensions
node = Utilities.netcdf.NcDim('node', mesh.Shape.nNode);
nele = Utilities.netcdf.NcDim('nele', mesh.nElement);

x    = Utilities.netcdf.NcVar('x', [node, nele], 'double');
y    = Utilities.netcdf.NcVar('y', [node, nele], 'double');
t    = Utilities.netcdf.NcVar('time', time, 'double');
h    = Utilities.netcdf.NcVar('h', [node, nele, time], 'double');
qx   = Utilities.netcdf.NcVar('qx', [node, nele, time], 'double');
qy   = Utilities.netcdf.NcVar('qy', [node, nele, time], 'double');
bot  = Utilities.netcdf.NcVar('bot', [node, nele], 'double');

file = Utilities.netcdf.NcFile(filename, [node, nele, time],...
    [x, y, t, h, qx, qy, bot]);

% initialize output file
file.CreateFile;

% set vertex location value
file.putVarPart('x',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('y',   [0, 0], [mesh.Shape.nNode, mesh.nElement], mesh.y);
file.putVarPart('bot', [0, 0], [mesh.Shape.nNode, mesh.nElement], phys.bot);

end% func