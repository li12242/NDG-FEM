function [h ,q] = SWESetUp
% Idealized dam break problem of 1D shallow water equation 
% 

physics = Utilities.varGroup;

% caseName = 'DamBreakDry';
% caseName = 'DamBreakWet';
caseName = 'ParabolicBowl';
% caseName = 'LakeAtRest';
% caseName = 'TsunamiRunup';

physics.incert('caseName', caseName);

% polynomial order and No. of elements
ndegree = 1; nele = 200;

% Set initial conditions
physics = SWEInit(physics, ndegree, nele);

% set output file
mesh = physics.getVal('mesh');
ncfile = CreateOutputFile(mesh);

% minimum water depth
hDry = 1e-6;
physics.incert('minDepth', hDry);
g = 9.81;
physics.incert('gravity', g);

% save mesh nodes
ncid = netcdf.open('SWE1D.nc','WRITE');
mesh_id = ncfile.varid(1);
netcdf.putVar(ncid,mesh_id,mesh.x(:))
mesh_id = ncfile.varid(7);
netcdf.putVar(ncid,mesh_id,mesh.x(mesh.vmapM(:)))
netcdf.close(ncid);
% Solve Problem
% [h, q] = SWESolverHrefinedWetDry(physics, ncfile);
[h, q] = SWESolver(physics, ncfile);

% filename = ['SWE1D_', num2str(ndegree), '_', num2str(nele), '.nc'];
% movefile('SWE1D.nc', filename)
end% func


