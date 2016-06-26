function phys = SWEInit2d(phys)
% Input: phys - structure variable contaions
%           | - casename : test case name
%           | - n : order of polynomials
%           | - ne : number of elements on each edge
%           | - meshType : 'tri' for triangle or 'quad' for quadrilaterals
% Output: phys
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
        [mesh, h, qx, qy, ftime, dt] = DamBreakDry(N, Ne, meshType);
    case 'DamBreakWet'
        [mesh, h, qx, qy, ftime, dt] = DamBreakWet(N, Ne, meshType);
    case 'ParabolicBowl'
        [mesh, h, qx, qy, ftime, dt] = ParabolicBowl(N, Ne, meshType);
end% switch

%% Assignments
% Assign the obtained variables to structure variable phys
phys.mesh  = mesh;
phys.h     = h;
phys.qx    = qx;
phys.qy    = qy;
phys.dt    = dt;
phys.ftime = ftime;

end% func

function [mesh, h, q, ftime, dt] = DamBreakDry(N, Ne, meshType)
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
        [EToV, VX, VY] = MeshGenRectangle2D(np, rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h = zeros(size(mesh.x));
q = zeros(size(mesh.x));

xc = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 10;
dt = 1e-2;
end% func

function [mesh, h, q, ftime, dt] = DamBreakWet(N, Ne)
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
        [EToV, VX, VY] = MeshGenRectangle2D(np, rmin,rmax);
        mesh = MultiRegions.RegionQuad(shape, EToV, VX, VY);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

%% Initial condition
h = 2*ones(size(mesh.x));
q = zeros(size(mesh.x));

xc = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 10;
dt = 1e-2;
end% func
