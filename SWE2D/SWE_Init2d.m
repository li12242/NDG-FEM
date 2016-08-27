function phys = SWE_Init2d(phys)
%% Function SWEInit2d
% Init mesh and other parameters for differents test cases.
% Input: phys - structure variable contaions, contains
%           |
%           | - casename : test case name
%           | - N : order of polynomials
%           | - Nx & Ny : number of points on each dimension
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
Nx = phys.nx;
Ny = phys.ny;
meshType = phys.meshType;
%% Initialization
% init mesh and initial condition of test case
switch casename
    case 'DamBreakDry'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            DamBreakDry(N, Nx, Ny, meshType);
    case 'DamBreakWet'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            DamBreakWet(N, Nx, Ny, meshType);
    case 'ParabolicBowl'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            ParabolicBowl(phys, N, Nx, Ny, meshType);
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

function [mesh, h, qx, qy, bot, ftime, dt, dx] = ParabolicBowl(phys, N, Nx, Ny, meshType)
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
        [VX,VY,EToV] = ...
            Utilities.Mesh.MeshGenTriangle2D(Nx,Ny,rmin,rmax,rmin,rmax,false);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = ...
            Utilities.Mesh.MeshGenRectangle2D(Nx,Ny,rmin,rmax,rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch
%% Initialize bottom level
r2    = (mesh.x.^2 + mesh.y.^2);
bot   = alpha*r2;
r2ext = (X+Y)/(alpha*(X^2 - Y^2));
%% Initialize variables
h      = 1/(X+Y) + alpha*(Y^2 - X^2).*r2/(X+Y)^2;
h(r2>r2ext) = 0;
qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1e-4;
ftime  = T;
dx     = min((rmax - rmin)./Nx/(N+1), (rmax - rmin)./Ny/(N+1));
end% func

function [mesh, h, qx, qy, bot, ftime, dt, dx] = DamBreakDry(N, Nx, Ny, meshType)
%% Initialize the mesh grid
% The grid range is and simulation ends at ftime (seconds).
% The elements of mesh grid can be triangles or quadrialterals and 
% be obtains a uniform mesh grid.
rmin  = 0; 
rmax  = 1000;
width = 200;
ftime = 20;
damPosition = 500;
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2,false);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h   = zeros(size(mesh.x));
qx  = zeros(size(mesh.x));
qy  = zeros(size(mesh.x));
bot = zeros(size(mesh.x));

xc        = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-2;
dx        = min((rmax - rmin)./Nx/(N+1), width./Ny/(N+1));

end% func

function [mesh, h, qx, qy, bot, ftime, dt, dx] = DamBreakWet(N, Nx, Ny, meshType)
%% Initialize the mesh grid
% The grid range is [-1, 1] and simulation ends at ftime (seconds).
% The elements of mesh grid can be triangles or quadrialterals and 
% be obtains a uniform mesh grid.
rmin  = 0; 
rmax  = 1000;
width = 200;
ftime = 20;
damPosition = 500;
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [VX,VY,EToV] = ...
            Utilities.Mesh.MeshGenTriangle2D(Nx,Ny,rmin,rmax,-width/2,width/2,false);
        mesh = MultiRegions.RegionTri(shape, EToV, VX, VY);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = ...
            Utilities.Mesh.MeshGenRectangle2D(Nx,Ny,rmin,rmax,-width/2,width/2);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h   = 2*ones(size(mesh.x));
qx  = zeros(size(mesh.x));
qy  = zeros(size(mesh.x));
bot = zeros(size(mesh.x));

xc        = mean(mesh.x); 
ind       = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-1;
dx        = min((rmax - rmin)./Nx/(N+1), width./Ny/(N+1));

end% func
