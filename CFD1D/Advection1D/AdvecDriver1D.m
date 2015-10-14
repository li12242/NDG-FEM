% Driver script for solving the 1D advection equations
% Globals1D;

% Order of polymomials used for approximation 
N = 2;

line = StdRegions.Line(N);

% Generate simple mesh
[Nv, VX, K, EToV] = Utilities.MeshGen1D(0.0,2.0,10);
BC = [2, 1; 3, Nv];

Mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% Initialize solver and construct grid and metric
% StartUp1D;

% Set initial conditions
u = scalar_type('scalar', Mesh);
u.setVal(sin(Mesh.x));

% Solve Problem
FinalTime = 10;
[u] = Advec1D(u,FinalTime);
