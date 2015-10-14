function TransDiffSetUp
% 1D Transient diffusion
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 40-42

% 基本单元
nOrder = 2;
line = StdRegions.Line(nOrder);

% 网格
xmin = 0; xmax = 1; K = 40;
[Nv, VX, ~, EToV] = Utilities.MeshGen1D(xmin,xmax,K);
BC = [2, 1;3, Nv];

% Multi-Region, create mesh
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

% condition
finalTime = 1.0;

% initial condition
c = TransDiffInit(mesh);

% solver
c = TransDiffSolver(mesh, c, finalTime);

% PostProcess
TranDiffPostProcess(mesh, c);

end% func