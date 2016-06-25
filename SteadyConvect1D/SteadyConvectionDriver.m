function u = SteadyConvectionDriver(N, nElement)
% solving steady convection problem by DGM
% \nabla u = sin(x)
% 

x1 = 0; x2 = 2*pi; % domain

% N = 2; nElement = 20;
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1];

% 
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

f = sin(mesh.x);

u = SteadyConvectionSolver(mesh, f);

% plot(mesh.x(:), u(:), 'b', mesh.x(:), cos(mesh.x(:)), 'r');
end% func