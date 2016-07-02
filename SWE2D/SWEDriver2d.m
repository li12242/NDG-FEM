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
% <matlab:edit(['/Users/mac/Documents/MATLAB/NDG-FEM/SWE2D/SWEInit2d.m']) 
% SWEInit2d>
% for further intialization

% Name of test case
casename = 'DamBreakDry';
% casename = 'DamBreakWet';
% casename = 'ParabolicBowl';
% casename = 'ParabolicBowl2'; % from Khan & Lai

% Order of polymomials used for approximation 
N = 1;
% Number of elements on each edge
Nx = 30;
Ny = 1;

% # Name of test case
phys.casename = casename;
phys.nx       = Nx + 1; % number of points
phys.ny       = Ny + 1; % number of points
phys.n        = N;
phys.meshType = 'quad';
phys.minDepth = 1e-2;
phys.gra      = 9.81;
% initialization
phys = SWEInit2d(phys);

%% Generate output
outfile = GenOutputFile('SWE2D', phys);

%% Solve the eqs
% The strucut variable 'phys' will be passed to function SWERHS2d and get
% the derived solutions of h (water depth) and q (water flux)

phys = SWESolve2d(phys, outfile);

%% Post process
DrawPoints(phys.mesh, phys.h, phys.qx, phys.qy);
end% func


