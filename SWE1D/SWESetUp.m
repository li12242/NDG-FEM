function [h ,q] = SWESetUp
% Idealized dam break problem of 1D shallow water equation 
% 
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 65

physics = Utilities.varGroup;

caseName = 'DamBreakDry';
switch caseName
    case 'DamBreakDry'
        FinalTime = 20; % Dam break
        x1 = 0; x2 = 1000; % Dam break
    case 'DamBreakWet'
        FinalTime = 20; % Dam break
        x1 = 0; x2 = 1000; % Dam break
    case 'FlowDump'
        FinalTime = 240; % Flow over dump
        x1 = 0; x2 = 25; % Flow over dump
    case 'ParabolicBowl'
        T = 269; FinalTime = T/2; % Parabolic Bowl
        x1 = -1000; x2 = 1000; % Parabolic Bowl
end% switch

physics.incert('FinalTime', FinalTime);
physics.incert('caseName', caseName);

% max order of polymomials
N = 2; nElement = 200;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1; 3,Nv];

% Initialize solver and construct grid and metric
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
physics.incert('mesh', mesh);

bedElevation = SetBed(physics);
physics.incert('bedElva', bedElevation);

% Set initial conditions
[h, q] = SWEInit(physics);
physics.incert('height', h);
physics.incert('flux', q);

% set output file
ncfile = CreateOutputFile(mesh);
ncid = netcdf.open('SWE1D.nc','WRITE');
mesh_id = ncfile.varid(1);
netcdf.putVar(ncid,mesh_id,mesh.x(:))
netcdf.close(ncid);
% Solve Problem
[h, q] = SWESolver(physics, ncfile);

% plot(mesh.x, h, 'o-')
end% func

function bedElevation = SetBed(physics)
mesh = physics.getVal('mesh');
% initialize the bed elevation according to the case name
switch physics.getVal('caseName')
    case 'DamBreakDry'
        bedElevation = zeros(size(mesh.x));
    case 'DamBreakWet'
        bedElevation = zeros(size(mesh.x));
    case 'FlowDump'
        bedElevation = zeros(size(mesh.x));
        flag = (mesh.x >= 8) & (mesh.x <=12);
        bedElevation(flag) = 0.2 - 0.05*(mesh.x(flag) -10).^2;
    case 'ParabolicBowl'
        a = 600; h0 = 10;
        bedElevation = h0.*(mesh.x.^2./a^2 - 1);
end% switch
end% func

function [h,q] = SWEInit(physics)
mesh = physics.getVal('mesh');
bedElva = physics.getVal('bedElva');
% initialize the bed elevation according to the case name
switch physics.getVal('caseName')
    case 'DamBreakDry'
        [h, q] = DamBreakInit(mesh, 2);
    case 'DamBreakWet'
        [h, q] = DamBreakInit(mesh, 1);
    case 'FlowDump'
        initCase = 1;
        [h, q] = FlowDumpInit(mesh, bedElva, initCase);
    case 'ParabolicBowl'
        [h, q] = ParaBowlInit(mesh, bedElva);
end% switch
end% func

function [h, q] = ParaBowlInit(mesh, bedElva)
q = zeros(size(mesh.x)); hDelta = 0.0;

g = 9.8; B = 5; h0 = 10; a = 600;
w = sqrt(2*g*h0)./a;
% z = zeros(size(mesh.x));
z = -(4*B*w).*mesh.x./(4*g);
h = z - bedElva;
h(h<hDelta) = hDelta;
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
        h(flag) = 0;
end% switch
end% func
