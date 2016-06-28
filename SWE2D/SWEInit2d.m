function phys = SWEInit2d(phys)
%% Function SWEInit2d
% Init mesh and other parameters for differents test cases.
% Input: phys - structure variable contaions, contains
%           |
%           | - casename : test case name
%           | - n : order of polynomials
%           | - ne : number of elements on each edge
%           | - meshType : 'tri' for triangle or 'quad' for quadrilaterals
% 
% Output: phys
%           |
%           | - ftime : final time
%           | - h : water depth
%           | - q : flow flux
%           | - mesh : mesh object
%           | - dt : delta time
% 

%% Get parameters
% test case name
casename = phys.casename;
% order of polynomials used for approximation
N = phys.n;
% number of elements on each edge
Ne = phys.ne;
meshType = phys.meshType;
%% Initialization
% init mesh and initial condition of test case
switch casename
    case 'DamBreakDry'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            DamBreakDry(N, Ne, meshType);
    case 'DamBreakWet'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            DamBreakWet(N, Ne, meshType);
    case 'ParabolicBowl'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            ParabolicBowl(phys, N, Ne, meshType);
end% switch

%% Assignments
% Assign the obtained variables to structure variable phys
phys.h     = h;
phys.qx    = qx;
phys.qy    = qy;
phys.dt    = dt;
phys.dx    = dx;
phys.ftime = ftime;
phys.bot   = botLevel; % bottom elevation
phys.mesh  = mesh;

end% func

function [mesh, h, qx, qy, bot, ftime, dt, dx] = ParabolicBowl(phys, N, Ne, meshType)
%% Parameters
g     = phys.gra;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
X     = 1;
Y     = -0.41884;

rmin  = -4000; 
rmax  =  4000;
%% Initialize mesh
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(Ne,rmin,rmax,0);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        np = Ne + 1;
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(np, rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch
%% Initialize bottom level
r2   = (mesh.x.^2 + mesh.y.^2);
bot = alpha*r2;

%% Initialize variables
h      = 1/(X+Y) + alpha*(Y^2 - X^2).*r2/(X+Y)^2 - bot;
h(h<0) = 0;
qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1;
ftime  = T/2;
dx     = (rmax - rmin)./Ne/(N+1);
end% func

function [mesh, h, qx, qy, botLevel, ftime, dt, dx] = DamBreakDry(N, Ne, meshType)
%% Initialize the mesh grid
% The grid range is and simulation ends at ftime (seconds).
% The elements of mesh grid can be triangles or quadrialterals and 
% be obtains a uniform mesh grid.
rmin = 0; rmax = 1000;
ftime = 20;
damPosition = 500;
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(Ne,rmin,rmax,0);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        np = Ne + 1;
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(np, rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h  = zeros(size(mesh.x));
qx = zeros(size(mesh.x));
qy = zeros(size(mesh.x));
botLevel = zeros(size(mesh.x));

xc        = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-1;
dx        = (rmax - rmin)./Ne/(N+1);
end% func

function [mesh, h, qx, qy, botLevel, ftime, dt, dx] = DamBreakWet(N, Ne, meshType)
%% Initialize the mesh grid
% The grid range is [-1, 1] and simulation ends at ftime (seconds).
% The elements of mesh grid can be triangles or quadrialterals and 
% be obtains a uniform mesh grid.
rmin = 0; rmax = 1000;
ftime = 20;
damPosition = 500;
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(Ne,rmin,rmax,0);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        np = Ne + 1;
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D(np, rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h  = 2*ones(size(mesh.x));
qx = zeros(size(mesh.x));
qy = zeros(size(mesh.x));
botLevel = zeros(size(mesh.x));

xc        = mean(mesh.x); 
ind       = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-1;
dx        = (rmax - rmin)./Ne/(N+1);
end% func
