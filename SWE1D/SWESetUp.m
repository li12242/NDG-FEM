function [h ,q] = SWESetUp(nele, hDry)

physics = Utilities.varGroup;

% caseName = 'DamBreakDry';
% caseName = 'DamBreakWet';
caseName = 'ParabolicBowl';
% caseName = 'LakeAtRest';
% caseName = 'TsunamiRunup';

physics.incert('caseName', caseName);
% minimum water depth
physics.incert('minDepth', hDry);
g = 9.81;
physics.incert('gravity', g);

% polynomial order and No. of elements
ndegree = 1; %nele = 400;

% Set initial conditions
physics = SWEInit(physics, ndegree, nele);

% set output file
mesh = physics.getVal('mesh');
ncfile = CreateOutputFile(mesh);

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

filename = ['SWE1D_', num2str(hDry), '_', num2str(nele), '.nc'];
movefile('SWE1D.nc', filename)
end% func


