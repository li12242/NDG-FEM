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
% casename = 'GaussMount';
casename = 'SquareMount';
% casename = 'SolidBody';

phys.casename = casename;
phys.meshtype = meshtype;

switch meshtype
    case 'quad'
        [mesh, phys.dx] = quadSolver(N, M);
    case 'tri'
        [mesh, phys.dx] = triSolver(N, M);
end% switch
phys.mesh = mesh;
% velocity field and simulation time
w      = 5*pi/6;
u      = w.*(0.5-mesh.y); 
v      = w.*(mesh.x-0.5);
phys.u = u;
phys.v = v;
FinalTime  = 2.4;
phys.ftime = FinalTime;

% initial field
phys = Convection2d_Init(phys);

% output file
filename   = ['Convection2D_', meshtype, '_', num2str(N),'_',num2str(M)];
ncfile     = CreateNetcdf(filename, mesh);
phys.file  = ncfile;

var = Convection2d_Solver(phys);
end% func

function file = CreateNetcdf(filename, mesh)

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
    case 'GaussMount'
        var = GaussMount(mesh);
    case 'SquareMount'
        var = SquareMount(mesh);
    case 'SolidBody'
        var = SolidBody(mesh);
end% switch
phys.var = var;
end% func

function var = SolidBody(mesh)
var = zeros(size(mesh.x));
% slotted cylinder
x0  = 0.5; 
y0  = 0.75; 
r0  = 0.15;
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
ind = ind & ((abs(mesh.x - x0)>=0.025) | (mesh.y >= 0.85));
var(ind) = 1.0;
% cone
x0  = 0.5; 
y0  = 0.25;
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = 1-r2(ind);
% hump
x0  = 0.25;
y0  = 0.5;
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = (1+cos(r2(ind)*pi))./4;
end% func

function var = GaussMount(mesh)
% Guass profile
var = zeros(size(mesh.x));
r0  = 0.15;
x0  = 0.25; 
y0  = 0.5;
r2  = sqrt((mesh.x-x0).^2+(mesh.y-y0).^2)./r0;
ind = ( r2<=1.0);
var(ind) = (1+cos(r2(ind)*pi))./2;
end% func

function var = SquareMount(mesh)
% square distribution
var = zeros(size(mesh.x));
b   = 1/12;
x0  = 0.50; xc = mean(mesh.x); 
y0  = 0.75; yc = mean(mesh.y);

flag = ( abs(xc - x0)<b & abs(yc - y0)<b );
var(:, flag) = 1;
end% func

%% solver for different type of element
function [mesh, dx] = triSolver(N, M)
% uniform triangle mesh
np   = M + 1;
[VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D...
    (np, np, 0, 1, 0, 1, false);

tri  = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTri(tri, EToV, VX, VY);
dx   = 2/np/(N+1);
end% func

function [mesh, dx] = quadSolver(N, M)
% uniform quadrilaterial mesh
np   = M + 1;
[EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
    (np, np, 0, 1, 0, 1);

quad = StdRegions.Quad(N);
mesh = MultiRegions.RegionQuad(quad, EToV, VX, VY);
dx   = 2/np/(N+1);
end% func

