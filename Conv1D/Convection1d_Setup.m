function Convection1d_Setup
% parameters
phys.N = 4; 
phys.K = 200;
phys.xlim = [0, 1];

casename = 'Advection';
% casename = 'AdvectionDiffusion';

% standard element and mesh object
line = StdRegions.Line(phys.N);
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D...
    (phys.xlim(1), phys.xlim(2), phys.K);
mesh = MultiRegions.RegionLine(line, EToV, VX);
% periodic boundary condition
% mesh.vmapP(1) = numel(mesh.x);
% mesh.vmapP(end) = 1;

% physical structure variable
phys.casename = casename;
phys.mesh = mesh;

phys.ftime = 1; % final time
phys.u = ones(size(mesh.x)); % flow rate
phys.Dx = 0.00; % diffusion parameter

% initial condition
phys = Convection1d_Init(phys);

% create output file
filename =['Convection1D_', num2str(phys.N), '_', num2str(phys.K)];
phys.file = Convection1d_output(filename, phys);

% plot initial condition
% plot(mesh.x, phys.var, 'r-.'); hold on;

% solver
phys = Convection1d_Solver(phys);

phys.file.CloseFile;
% postprocess
% plot(mesh.x, phys.var, 'b-.');
% ylim([-0.05, 1.1])
end% func

function outfile = Convection1d_output(filename, phys)

mesh = phys.mesh;
shape = mesh.Shape;

time = Utilities.NetcdfClass.NcDim('time', 0); % unlimited dimensions
ne = Utilities.NetcdfClass.NcDim('node', mesh.nElement);
np = Utilities.NetcdfClass.NcDim('np', shape.nNode);

x = Utilities.NetcdfClass.NcVar('loc', [np, ne], 'double');
t = Utilities.NetcdfClass.NcVar('time', time, 'double');
var = Utilities.NetcdfClass.NcVar('var', [np, ne, time], 'double');

outfile = Utilities.NetcdfClass.NcFile...
    (filename, [np, ne, time], [x, t, var]);

% initialize output file
outfile.CreateFile;

% set vertex location value
outfile.putVarPart('loc', [0,0], [shape.nNode, mesh.nElement], mesh.x);
end% func