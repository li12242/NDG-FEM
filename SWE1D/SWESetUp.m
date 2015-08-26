function [h ,q] = SWESetUp
% Idealized dam break problem of 1D shallow water equation 
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 65

% max order of polymomials
N = 2; nElement = 100;
len = 25;
% len = 1000;
% assert(mod(nElement,2) == 0)
% Generate simple mesh
[Nv, VX, K, EToV] = Utilities.MeshGen1D(0.0,len,nElement);
BC = [2,1; 3,Nv];

% Initialize solver and construct grid and metric
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
bedElevation = SetBed(mesh);

% Solve Problem
FinalTime = 50.0;

% Set initial conditions
[h, q] = SWEInit(mesh, bedElevation);

[h, q] = SWESolver(mesh, h, q, bedElevation, FinalTime);

% plot(mesh.x, h, 'o-')
end% func

function bedElevation = SetBed(mesh)
bedElevation = zeros(size(mesh.x));
flag = (mesh.x >= 8) & (mesh.x <=12);
bedElevation(flag) = 0.2 - 0.05*(mesh.x(flag) -10).^2;
end% func

function [h,q] = SWEInit(mesh, bedElva)
initCase = 1;
[h, q] = FlowDumpInit(mesh, bedElva, initCase);
% [h, q] = DamBreakInit(mesh, initCase);
end% func

function [h, q] = FlowDumpInit(mesh, bedElva, initCase)
switch initCase
    case 1 % subcritical flow
        h = 0.5.*ones(size(mesh.x))- bedElva; 
        q = 0.18.*ones(size(mesh.x));
    case 2 % supercritical flow
        h = 2.0.*ones(size(mesh.x))- bedElva; 
        q = 25.0567.*ones(size(mesh.x));
    case 3 % transcritical flow
        h = 0.33.*ones(size(mesh.x))- bedElva;
        q = 0.18.*ones(size(mesh.x));
end% switch
end% func

function [h, q] = DamBreakInit(mesh, initCase)
% Idealized dam break problem of 1D shallow water equation
h = 10.*ones(size(mesh.x)); q = zeros(size(mesh.x));
damPosition = 500;
switch initCase
    case 1 % wet bed
        flag = mesh.x > damPosition;
        h(flag) = 2;
    case 2 % dry bed
        flag = mesh.x > damPosition;
        h(flag) = 10^-12;
end% switch
end% func
