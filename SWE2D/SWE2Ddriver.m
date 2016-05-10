% Driver script for solving the 1D shallow water equation

% Order of polymomials used for approximation 
N = 3;
% read triangle mesh
[EToV, VX, VY, EToR, BC] = Utilities.Mesh.MeshReaderTriangle('SWE2D/grid/untitled');
% [Nv, VX, K, EToV] = Utilities.MeshGen1D(0.0,1000,nElement);

% Initialize solver and construct grid and metric
tri = StdRegions.Triangle(N);
mesh = MultiRegions.RegionTriBC(tri, EToV, VX, VY, BC);

% you can also load data
load('/Users/mac/Documents/Model/NDG-FEM/SWE2D/SWE2D_N3.mat')

% Solve Problem
FinalTime = 10.0;

% Set initial conditions
Q = zeros([size(mesh.x),3]);
H = ones(tri.nNode,1)*(EToR == 21)';
Q(:,:,1) = H*5+5;

% tic
[Q] = SWE2D(mesh, Q, FinalTime, @SWE2DRHS_LF);
