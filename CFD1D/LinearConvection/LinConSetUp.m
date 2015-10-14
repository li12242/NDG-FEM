function [c, error] = LinConSetUp(nOrder,K)
% 1D Linear Convection Test Case
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 49-52

% 基本单元
% nOrder = 3;
line = StdRegions.Line(nOrder);

% 网格
xmin = 0; xmax = 1; %K = 100;
[Nv, VX, ~, EToV] = Utilities.MeshGen1D(xmin,xmax,K);
BC = [2, 1;3, Nv];

% Multi-Region, create mesh
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% condition
finalTime = 0.6;
u = 1.0;

% initial condition
c = LinConInit(mesh);

% solver
c = LinConSolver(mesh, c, u,finalTime);

% postprocess
% LinConPostProcess(mesh, c)

% calculate error
error = LinConError(mesh, c);
end% func