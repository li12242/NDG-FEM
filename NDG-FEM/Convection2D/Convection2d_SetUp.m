function var = Convection2d_SetUp(meshtype, N, M)
% 2D convection problem
% dc/dt + d(uc)/dx + d(vc)/dy = 0
% INPUT:
%   meshtype - element type of mesh grid
%   N        - degree of polynomial
%   M        - No. of elements on each edge
% OUTPUT:
%   var - scalar variable

% name of different test case
% casename = 'Rotation';
casename = 'AdvectionDiffusion';

phys.casename = casename;
phys.meshtype = meshtype;

switch meshtype
    case 'quad'
        [mesh, phys.dx] = quadSolver(N, M);
    case 'tri'
        [mesh, phys.dx] = triSolver(N, M);
end% switch
phys.mesh = mesh;
% initial field
phys = Convection2d_Init(phys);

% output file
filename   = ['Convection2D_', meshtype, '_', num2str(N),'_',num2str(M)];
ncfile     = Convection2d_output(filename, mesh);
phys.file  = ncfile;

var = Convection2d_Solver(phys);
end% func

function file = Convection2d_output(filename, mesh)

time = Utilities.NetcdfClass.NcDim('time', 0); % unlimited dimensions
np   = Utilities.NetcdfClass.NcDim('np', mesh.Shape.nNode);
ne   = Utilities.NetcdfClass.NcDim('ne', mesh.nElement);

x    = Utilities.NetcdfClass.NcVar('x', [np, ne], 'double');
y    = Utilities.NetcdfClass.NcVar('y', [np, ne], 'double');
t    = Utilities.NetcdfClass.NcVar('time', time, 'double');
var  = Utilities.NetcdfClass.NcVar('var', [np, ne, time], 'double');

file = Utilities.NetcdfClass.NcFile...
    (filename,[np, ne, time],[x, y, t, var]);

% initialize output file
file.CreateFile;

% set vertex location value
file.putVarPart('x', [0,0], [mesh.Shape.nNode, mesh.nElement], mesh.x);
file.putVarPart('y', [0,0], [mesh.Shape.nNode, mesh.nElement], mesh.y);

end% func

%% Initial functions
function phys = Convection2d_Init(phys)
mesh = phys.mesh;
switch phys.casename
    case 'Rotation'
        [var,u,v,ftime,Dx,Dy] = Rotation_Init(mesh);
    case 'AdvectionDiffusion'
        [var,u,v,ftime,Dx,Dy] = AdvectDiff_Init(mesh);
end% switch
phys.var = var;
phys.u = u;
phys.v = v;
phys.ftime = ftime;
phys.Dx = Dx;
phys.Dy = Dy;
end% func

function [var, u, v, ftime, Dx, Dy] = Rotation_Init(mesh)
% Guass profile
var = zeros(size(mesh.x));
r0  = 0.15;
x0  = 0.25; 
y0  = 0.5;
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = (1+cos(r2(ind)*pi))./2;
% velocity field
w      = 5*pi/6;
u      = w.*(0.5-mesh.y); 
v      = w.*(mesh.x-0.5);
% final time
ftime  = 2.4;
% diffusion parameter
Dx = 0; 
Dy = 0;
end% func

function [var, u, v, ftime, Dx, Dy] = AdvectDiff_Init(mesh)
% diffusion parameter
Dx = 0.01; 
Dy = 0.01;
% velocity field
u = 0.5.*ones(size(mesh.x));
v = 0.5.*ones(size(mesh.x));
% Guass profile
x0  = -0.5;
y0  = -0.5;
t = -(mesh.x-x0).^2/Dx -(mesh.y-y0).^2/Dy;
var = exp(t);
% final time
ftime = 0.5;
end

%% solver for different type of element
function [mesh, dx] = triSolver(N, M)
% uniform triangle mesh
np   = M + 1;
[VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D...
    (np, np, -1, 1, -1, 1, false);

tri  = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
wv = sum(tri.M);
vol = (wv*mesh.J);
dx = min(sqrt(vol/pi));
end% func

function [mesh, dx] = quadSolver(N, M)
% uniform quadrilaterial mesh
np   = M + 1;
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
    (np, np, -1, 1, -1, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);
wv = sum(quad.M);
vol = (wv*mesh.J);
dx = min(sqrt(vol/pi));
end% func

