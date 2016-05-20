function Convection1DSetup

% parameters
N = 2; M = 100;
x1 = 0; x2 = 1;

% set mesh
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, M);
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLine(line, EToV, VX);

% initial condition
var = Convection1DInit(mesh, x1, x2);
FinalTime = 2; % time
a = 1; % flow rate

% boundary condition
mesh.vmapP(1) = numel(mesh.x); mesh.vmapP(end) = 1;

% solver
var = Convection1DSolver(mesh, var, FinalTime, a);

% postprocess
end% func