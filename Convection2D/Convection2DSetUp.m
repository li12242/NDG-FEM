function [mesh, var] = Convection2DSetUp(N, M)
% 2D convection problem
% dc/dt + d(uc)/dx + d(vc)/dy = 0
% Input:
%   N - degree of polynomial
%   M - No. of elements on each edge
% Output:
%   mesh - mesh object
%   var - scalar variable

% mesh = quadSolver(N, M);
mesh = triSolver(N, M);
var = ConvectionInit(mesh);

w = 5*pi/6;
% flow rate in domain, [u, v]
u = -w.*mesh.y; v = w.*mesh.x;

FinalTime = 2.4;
filename = ['Convection2D_', num2str(N),'_',num2str(M)];
outfile = CreateNetcdf(filename, mesh);

var = Convection2DSolver(mesh, var, FinalTime, u, v, outfile);
% postprocess(mesh, var)
end% func

function file = CreateNetcdf(filename, mesh)

time = Utilities.NetcdfClass.NcDim('time', 0); % unlimited dimensions
node = Utilities.NetcdfClass.NcDim('node', mesh.nNode);

x    = Utilities.NetcdfClass.NcVar('x', node, 'double');
y    = Utilities.NetcdfClass.NcVar('y', node, 'double');
t    = Utilities.NetcdfClass.NcVar('time', time, 'double');
var  = Utilities.NetcdfClass.NcVar('var', [node, time], 'double');

file = Utilities.NetcdfClass.NcFile(filename,[node, time],[x, y, t, var]);

% initialize output file
file.CreateFile;

% set vertex location value
file.putVarPart('x', 0, mesh.nNode, mesh.x);
file.putVarPart('y', 0, mesh.nNode, mesh.y);

end% func

function var = ConvectionInit(mesh)
% % Guass profile
% sigma = 125*1e3/33^2; 
% xc = 0; yc = 3/5;
% var = exp(-sigma.*( (mesh.x - xc).^2 + (mesh.y - yc).^2) );

% square distribution
var = zeros(size(mesh.x));
b   = 1/12;
x0  = 0;   xc = mean(mesh.x); 
y0  = 3/4; yc = mean(mesh.y);

flag = ( abs(xc - x0)<b & abs(yc - y0)<b );
var(:, flag) = 1;
end% func

function mesh = triSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('Convection2D/mesh/triangle');
np   = M + 1;
[VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(np, np, -1, 1, -1, 1, 1);

tri  = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderQuad('Convection2D/mesh/quad');

% uniform mesh
np   = M + 1;
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(np, np, -1, 1, -1, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func

function postprocess(mesh, var)
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM))
end% func

