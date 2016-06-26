function SWEDriver2d
%% Introduction
% Solving the 2D shallow water equation
% Todo:
% 
% # time step formula
% # 

%% Initialization
% Assignment parameters of the test case to the struct variable `phys`. 
% It includes
% 
% # Name of test case
% # Number of elements on each edge
% # Order of polynomials for approximations
% # mesh type
% # threshold of water depth to determine dry elements
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
phys.minDepth = 1e-4;
phys.gra = 9.81;
% initialization
phys = SWEInit2d(phys);

% % test
% mesh = phys.mesh;
% h = phys.h;
% plot3(mesh.x, mesh.y, h, 'b.');
% plot(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM))

%% Solve the eqs
% The strucut variable 'phys' will be passed to function SWERHS2d and get
% the derived solutions of h (water depth) and q (water flux)

phys = SWERHS2d(phys);

%% Post process

end% func


