function SWEDriver2d
% Driver script for solving the 1D shallow water equation

%% Initialization
% Assignment parameters of the test case to the struct variable `phys`. 
% It includes
% 
% # Name of test case
% # Number of elements on each edge
% # Order of polynomials for approximations
% # mesh type
% 
% The struct variable will be put into function 
% <matlab:edit(['/Users/mac/Documents/MATLAB/NDG-FEM/SWE2D/SWEInit2d.m']) SWEInit2d>
% for further intialization

% Name of test case
casename = 'DamBreakDry';
% casename = 'DamBreakWet';
% casename = 'ParabolicBowl';

% Order of polymomials used for approximation 
N = 3;
% Number of elements on each edge
Ne = 20;

% # Name of test case
phys.casename = casename;
phys.ne = Ne;
phys.n  = N;
phys.meshType = 'tri';

phys = SWEInit2d(phys);

% % test
% mesh = phys.mesh;
% h = phys.h;
% plot3(mesh.x, mesh.y, h, 'b.');
% plot(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM))

%% Solve the eqs


end% func


