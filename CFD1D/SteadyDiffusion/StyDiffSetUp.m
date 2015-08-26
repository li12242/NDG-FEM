function StyDiffSetUp
% 1D Steady diffusion
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 43-46

% 基本单元
nOrder = 1;
line = StdRegions.Line(nOrder);

% 网格
xmin = 0; xmax = 1; K = 40;
[Nv, VX, ~, EToV] = Utilities.MeshGen1D(xmin,xmax,K);
BC = [2, 1;3, Nv];

% Multi-Region, create mesh
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% Init Condition
[c, q] = StyDiffInit(mesh);

% error
TOLERR = 10^-6;

% Solver
c = StyDiffSolver(mesh, c, q, TOLERR);

% PostProcess
plot(mesh.x, c, 'o');
end% func