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
        FinalTime = 200; % Flow over dump
        x1 = 0; x2 = 25; % Flow over dump
    case 'ParabolicBowl'
        T = 269; FinalTime = T; % Parabolic Bowl
        x1 = -1000; x2 = 1000; % Parabolic Bowl
    case 'LakeAtRest'
        FinalTime = 0.5;
        x1 = 0; x2 = 1;
    case 'TsunamiRunup'
%         FinalTime = 240;
        FinalTime = 10;
        x1 = -500; x2 = 50000;
end% switch

physics.incert('FinalTime', FinalTime);
physics.incert('caseName', caseName);

% max order of polymomials
N = 1; nElement = 200;
[Nv, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1; 3,Nv];

% Initialize solver and construct grid and metric
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
physics.incert('mesh', mesh);
physics.incert('VX', VX);
physics.incert('EToV', EToV);

% Set initial conditions
[h, q] = SWEInit(mesh, physics);
physics.incert('height', h);
physics.incert('flux', q);

% set output file
ncfile = CreateOutputFile(mesh);
ncid = netcdf.open('SWE1D.nc','WRITE');
mesh_id = ncfile.varid(1);
netcdf.putVar(ncid,mesh_id,mesh.x(:))
mesh_id = ncfile.varid(7);
netcdf.putVar(ncid,mesh_id,mesh.x(mesh.vmapM(:)))
netcdf.close(ncid);
% Solve Problem
[h, q] = SWESolverHrefinedWetDry(physics, ncfile);
% [h, q] = SWESolver(physics, ncfile);
end% func


