function SWE_Driver2d(type)
%% Introduction
% Solving the 2D shallow water equation
% 

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
% casename = 'DamBreakDry';
casename = 'DamBreakWet';
% casename = 'ParabolicBowl';
% casename = 'PartialDamBreak';
% casename = 'FlowOver3BumpsUniform';
% casename = 'FlowOver3Bumps';
% casename = 'TsuamiRunup';
% casename = 'ObliqueHydraulicJump';

% Order of polymomials used for approximation 
N = 1;
% Number of elements on each edge
Nx = 600;
Ny = 1;

% # Name of test case
phys.casename = casename;
phys.nx       = Nx + 1; % number of points
phys.ny       = Ny + 1; % number of points
phys.n        = N;
phys.meshType = type;
phys.gra      = 9.81;

if (strncmp(phys.casename, 'TsuamiRunup', 11)) % spicific coefficient
    phys.ManningCoeff = 0.01;
    phys.minDepth = 5e-4;
    phys.minht    = 5e-4;
else % coefficient for other test case
    phys.ManningCoeff = 0;
    phys.minDepth = 1e-4;
    phys.minht    = 1e-4;
end% if

% initialization
phys = SWE_Init2d(phys);

%% Generate output
outfile = SWE_GenOutputFile2d([casename,'_VA_', ...
    phys.meshType,'_', num2str(Nx), '_', num2str(N)], phys);

%% Solve the eqs
% The strucut variable 'phys' will be passed to function SWERHS2d 
% and get the derived solutions of h (water depth) and q (water flux)

phys = SWE_Solve2d(phys, outfile);

%% Post process
% outfile.CloseFile;
end% func


