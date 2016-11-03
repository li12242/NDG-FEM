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

% Get parameters
casename = phys.casename; % test case name

% order of polynomials used for approximation
N = phys.n;
% number of elements on each edge
Nx = phys.nx;
Ny = phys.ny;
meshType = phys.meshType;

% Initialization for spicific test case
% init mesh and initial condition
switch casename
    case 'DamBreakDry'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            DamBreakDry(N, Nx, Ny, meshType);
    case 'DamBreakWet'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx, hin, hout] = ...
            DamBreakWet(N, Nx, Ny, meshType);
        phys.hin = hin;
        phys.hout = hout;
    case 'ParabolicBowl'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            ParabolicBowl(phys, N, Nx, Ny, meshType);
    case 'PartialDamBreak'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            PartialDamBreak(N, meshType);
    case 'FlowOver3BumpsUniform'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            FlowOver3BumpsUniform(N, Nx, Ny, meshType);
    case 'FlowOver3Bumps'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
            FlowOver3Bumps(N, meshType);
    case 'TsuamiRunup'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx, inWave] = ...
            TsuamiRunup(N, Nx, Ny, meshType);
        phys.inWave = inWave;
    case 'ObliqueHydraulicJump'
        [mesh, h, qx, qy, botLevel, ftime, dt, dx, hin, uin] = ...
            ObliqueHydraulicJump(N, meshType);
        phys.hin = hin;
        phys.uin = uin;
end% switch

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

%% ObliqueHydraulicJump
% For details please refer to Ghostine, Kesserwani and Mose (2009)
function [mesh, h, qx, qy, botLevel, ftime, dt, dx, hin, uin] = ...
            ObliqueHydraulicJump(N, meshType)

switch meshType
    case 'tri'
        filename = 'ObliqueHydraulicJumpTriCorser';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    case 'quad'
        filename = 'ObliqueHydraulicJumpQuad';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    otherwise
        error('Unknown mesh type "%s"', meshType);
end% switch

% Initialize mesh
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
end

% inflow condition
hin    = 1.0;
uin    = 8.57;

% Bottom level
botLevel = zeros(size(mesh.x));

% Initialize variables
h      = ones(size(mesh.x)).*hin;
qx     = ones(size(mesh.x)).*uin.*hin;
qy     = zeros(size(mesh.x));
dt     = 1e-3;
ftime  = 15;
% element length scal
w      = sum(shape.M)';
area   = sum(repmat(w, 1, mesh.nElement).*mesh.J);
dx     = min( sqrt(area/pi)/(N+1) );
end% func

%% Tsuami Runup on complex shoreline
function [mesh, h, qx, qy, botLevel, ftime, dt, dx, inWave] = ...
    TsuamiRunup(N, Nx, Ny, meshType)
% Parameters
rmin = 0;
rmax = 5.488;
smin = 0;
smax = 3.402;

% Boundary condition
BC = zeros((Nx-1)*2+(Ny-1)*2, 3);
stride = [Nx-1, Nx-1, Ny-1, Ny-1];
% south - wall
start = 1; stop = start + stride(1) - 1;
BC(start:stop, 1) = 1;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'south');
% north - wall
start = start + stride(1); stop = start + stride(2) - 1;
BC(start:stop, 1) = 1;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'north');
% east - wall
start = start + stride(2); stop = start + stride(3) - 1;
BC(start:stop, 1) = 1;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'east');
% west - inflow
start = start + stride(3); stop = start + stride(4) - 1;
BC(start:stop, 1) = 2;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'west');

% Create mesh grid
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [EToV,VX,VY] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,smin,smax,false);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,smin,smax);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
end% switch

% Initialize bottom level
bathymetryFile = 'SWE2D/mesh/TsunamiRunup/Benchmark_2_Bathymetry.txt';
fp = fopen(bathymetryFile);
fgetl(fp);
data = fscanf(fp, '%e %e %e', [3, inf]);
fclose(fp);
interp = TriScatteredInterp(data(1,:)',data(2,:)',-data(3,:)','linear');
botLevel = interp(mesh.x, mesh.y);
% Initialize variables
zeta = zeros(size(mesh.x));
h    = zeta - botLevel;
h(h<0) = 0;
qx   = zeros(size(mesh.x));
qy   = zeros(size(mesh.x));
dx   = min((rmax - rmin)./Nx/(N+1), (smax - smin)./Ny/(N+1));
dt   = 0.00010;

ftime = 22.5;

% Inflow water elevation
inputFile = 'SWE2D/mesh/TsunamiRunup/Benchmark_2_input.txt';
fp = fopen(inputFile);
fgetl(fp);
inWave = fscanf(fp, '%e %e', [2, inf]);
fclose(fp);
end

%% FlowOver3Bumps
function [mesh, h, qx, qy, botLevel, ftime, dt, dx] ...
    = FlowOver3Bumps(N, meshType)

switch meshType
    case 'tri'
        filename = 'FlowOver3BumpTri';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    case 'quad'
        filename = 'FlowOver3BumpQuad';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    otherwise
        error('Unknown mesh type "%s"', meshType);
end% switch

NTOL = 1e-12;
% substitute 1 (Wall) for 19, 3 (Outflow) for 18
bind = abs(BC(:, 1) - 1) < NTOL;
BC(bind,  1) = 1;

% Initialize mesh
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
end
% Initialize bottom level
botLevel = zeros(size(mesh.x));

x0 = 30; y0 = 7.5; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 30; y0 = -7.5; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 47.5; y0 = 0; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1;
botLevel(sk) = 2.8*(1 - r2(sk));

% Initialize variables
h      = zeros(size(mesh.x));
wetind = abs(EToR(:, 1) - 1) < NTOL;
h(:,  wetind) = 1.875;
h(:, ~wetind) = 1e-12;
qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1e-4;
ftime  = 10;
% element length scal
w      = sum(shape.M)';
area   = sum(repmat(w, 1, mesh.nElement).*mesh.J);
dx     = min( area/pi/(N+1) );

end% func

% Uniform grids
function [mesh, h, qx, qy, botLevel, ftime, dt, dx] =...
    FlowOver3BumpsUniform(N, Nx, Ny, meshType)
%% Parameters
rmin  = 0; 
rmax  = 75;
width = 30;
damPosition = 16;

% Initialize mesh
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [EToV,VX,VY] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2,false);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, []);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, []);
    otherwise
        error('error: unknown mesh type "%s"', meshType)
end% switch

% Boundary condition
ind = find(mesh.vmapM == mesh.vmapP);
Nfp = N+1; Num = numel(ind);
mesh.mapW = reshape(ind, Nfp, Num/Nfp);

% Initialize bottom level
botLevel = zeros(size(mesh.x));

x0 = 30; y0 = 7.5; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 30; y0 = -7.5; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 47.5; y0 = 0; r0 = 10;
r2 = sqrt((mesh.x-x0).^2 + (mesh.y-y0).^2)./r0;
sk = r2 < 1;
botLevel(sk) = 2.8*(1 - r2(sk));

% Initialize variables
h         = zeros(size(mesh.x));
xc        = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 1.875;

qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1e-4;
ftime  = 20;
dx     = min((rmax - rmin)./Nx/(N+1), (width)./Ny/(N+1));

end% func

%% Partial DamBreak
function [mesh, h, qx, qy, botLevel, ftime, dt, dx] = ...
    PartialDamBreak(N, meshtype)
switch meshtype
    case 'tri'
        filename = 'PartialDamBreakTri';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    case 'quad'
        filename = 'PartialDamBreakQuad';
        [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
    otherwise
        error('Unknown mesh type "%s"', meshType);
end% switch

NTOL = 1e-12;
% substitute 1 (Wall) for 19, 3 (Outflow) for 18
bind = abs(BC(:, 1) - 19) < NTOL;
BC(bind,  1) = 1;
BC(~bind, 1) = 3;

% Initialize mesh
switch meshtype
    case 'tri'
        shape = StdRegions.Triangle(N);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
end
% Initialize bottom level
botLevel = zeros(size(mesh.x));

% Initialize variables
h      = zeros(size(mesh.x));
wetind = abs(EToR(:, 1) - 21) < NTOL;
h(:,  wetind) = 5;
h(:, ~wetind) = 1e-12;
qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1e-4;
ftime  = 10;
% element length scal
w      = sum(shape.M)';
area   = sum(repmat(w, 1, mesh.nElement).*mesh.J);
dx     = min( area/pi/(N+1) );

end% func

%% ParabolicBowl
function [mesh, h, qx, qy, bot, ftime, dt, dx] = ...
    ParabolicBowl(phys, N, Nx, Ny, meshType)
% Initialize for parameters.
g     = phys.gra;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
X     = 1;
Y     = -0.41884;

rmin  = -4000; 
rmax  =  4000;

% Boundary condition
BC = zeros((Nx-1)*2+(Ny-1)*2, 3);
stride = [Ny-1, Ny-1, Nx-1, Nx-1]; % No. of nodes on each BC
% west
start = 1; stop = start + stride(1) - 1;
BC(start:stop, 1) = 3;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'west');

% east
start = start + stride(1); stop = start + stride(2) - 1;
BC(start:stop, 1) = 3;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'east');

% north
start = start + stride(2); stop = start + stride(3) - 1;
BC(start:stop, 1) = 3;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'north');

% south
start = start + stride(3); stop = start + stride(4) - 1;
BC(start:stop, 1) = 3;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'south');


% Initialize mesh
switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [EToV,VX,VY] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,rmin,rmax,false);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,rmin,rmax);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch
% Initialize bottom level
r2    = (mesh.x.^2 + mesh.y.^2);
bot   = alpha*r2;
r2ext = (X+Y)/(alpha*(X^2 - Y^2));
% Initialize variables
h      = 1/(X+Y) + alpha*(Y^2 - X^2).*r2/(X+Y)^2;
h(r2>r2ext) = 0;
qx     = zeros(size(mesh.x));
qy     = zeros(size(mesh.x));
dt     = 1e-4;
ftime  = T;
dx     = min((rmax - rmin)./Nx/(N+1), (rmax - rmin)./Ny/(N+1));
end% func

%% Ideal DamBreak
function [mesh, h, qx, qy, bot, ftime, dt, dx] ...
    = DamBreakDry(N, Nx, Ny, meshType)
% Initialize the mesh grid.
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
        [EToV,VX,VY] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2,false);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, []);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, []);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

% Initial condition
h   = zeros(size(mesh.x));
qx  = zeros(size(mesh.x));
qy  = zeros(size(mesh.x));
bot = zeros(size(mesh.x));

xc        = mean(mesh.x); ind = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-2;
dx        = min((rmax - rmin)./Nx/(N+1), width./Ny/(N+1));

end% func

function [mesh, h, qx, qy, bot, ftime, dt, dx, hin, hout] ...
    = DamBreakWet(N, Nx, Ny, meshType)
% Initialize the mesh grid.
% The grid range is [-1, 1] and simulation ends at ftime (seconds).
% The elements of mesh grid can be triangles or quadrialterals and 
% be obtains a uniform mesh grid.
rmin  = 0; 
rmax  = 1000;
width = 200;
ftime = 20;
damPosition = 500;

% BC list
BC = zeros((Ny-1)*2, 3);
stride = [Ny-1, Ny-1];
% west - inflow
start = 1; stop = start + stride(1) - 1;
BC(start:stop, 1) = 2;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'west');

% east - outflow
start = start + stride(1); stop = start + stride(2) - 1;
BC(start:stop, 1) = 3;
[BC(start:stop, 2), BC(start:stop, 3)] = UniformMeshNodes(Nx, Ny, 'east');

switch meshType
    case 'tri'
        shape = StdRegions.Triangle(N);
        [EToV, VX,VY] = Utilities.Mesh.MeshGenTriangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2,false);
%         filename = 'DamBreakWet';
%         [EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReader2D(filename);
        mesh = MultiRegions.RegionTriBC(shape, EToV, VX, VY, BC);
    case 'quad'
        shape = StdRegions.Quad(N);
        [EToV, VX, VY] = Utilities.Mesh.MeshGenRectangle2D...
            (Nx,Ny,rmin,rmax,-width/2,width/2);
        mesh = MultiRegions.RegionQuadBC(shape, EToV, VX, VY, BC);
    otherwise
        error('DamBreakDry error: unknown mesh type "%s"', meshType)
end% switch

% Initial condition
h   = 2*ones(size(mesh.x));
qx  = zeros(size(mesh.x));
qy  = zeros(size(mesh.x));
bot = zeros(size(mesh.x));

% Boundary condition
hin  = 10;
hout = 2;

xc        = mean(mesh.x); 
ind       = xc < damPosition;
h(:, ind) = 10;
dt        = 1e-4;
dx        = min((rmax - rmin)./Nx/(N+1), width./Ny/(N+1));

end% func

function [id1, id2] = UniformMeshNodes(Nx, Ny, location)
switch location
    case 'west'
        id1 = linspace(1, Nx*(Ny-2)+1, Ny-1);
        id2 = linspace(Nx+1, Nx*(Ny-1)+1, Ny-1);
    case 'east'
        id1 = linspace(Nx, Nx*(Ny-1), Ny-1);
        id2 = linspace(Nx*2, Nx*Ny, Ny-1);
    case 'north'
        id1 = (1:(Nx-1))+Nx*(Ny-1);
        id2 = (1:(Nx-1))+Nx*(Ny-1)+1;
    case 'south'
        id1 = 1:(Nx-1);
        id2 = (1:(Nx-1))+1;
end% switch
end% func
