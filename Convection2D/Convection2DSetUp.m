function [mesh, var] = Convection2DSetUp
% 2D convection problem
% dc/dt + d(uc)/dx + d(vc)/dy = 0
% Input:
%   N - degree of polynomial
%   M - No. of elements on each edge
% 

N = 2; M = 40;

mesh = quadSolver(N, M);
% mesh = triSolver(N, M);
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

function outfile = CreateNetcdf(filename, mesh)

time = Utilities.Netcdf.dimobj('time', 0); % unlimited dimensions
node = Utilities.Netcdf.dimobj('node', mesh.nNode);

x = Utilities.Netcdf.varobj('x', node, 'double');
y = Utilities.Netcdf.varobj('y', node, 'double');
t = Utilities.Netcdf.varobj('time', time, 'double');
var = Utilities.Netcdf.varobj('var', [node, time], 'double');

outfile = Utilities.Netcdf.fileobj(filename, [node, time], [x, y, t, var]);

% initialize output file
outfile.createFile;

% set vertex location value
outfile.putVarPart('x', 0, mesh.nNode, mesh.x);
outfile.putVarPart('y', 0, mesh.nNode, mesh.y);

end% func

function var = ConvectionInit(mesh)
% var = ones(size(mesh.x));
% var = mesh.x;
% xc = mean(mesh.x);
% left = xc < 0.5; right = xc > 0.5;
% var = sin(pi*mesh.x);%.*sin(2*pi*mesh.y);
% var = zeros(size(mesh.x));
% var(:,left) = 1; var(:,right) = 0;

sigma = 125*1e3/33^2; 
xc = 0; yc = 3/5;
var = exp(-sigma.*( (mesh.x - xc).^2 + (mesh.y - yc).^2) );
end% func

function mesh = triSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('Convection2D/mesh/triangle');

[VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(M, -1, 1, 1);

tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
end% func

function mesh = quadSolver(N, M)

% read triangle mesh
% [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderQuad('Convection2D/mesh/quad');

% uniform mesh
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(M+1, -1, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);

end% func

function postprocess(mesh, var)
plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM))
end% func

