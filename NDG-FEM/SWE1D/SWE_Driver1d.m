function SWE_Driver1d
%% Parameters
% caseName = 'DamBreakDry';
% caseName = 'DamBreakWet';
% caseName = 'ParabolicBowl';
% caseName = 'LakeAtRest';
% caseName = 'TsunamiRunup';

% % sponge layer test
% caseName = 'WiderGaussianMound';
% caseName = 'LocalGaussianMound';
% caseName = 'WiderMovingMound';
caseName = 'LocalMovingMound_SpongeLayer';

% minimum water depth
minDepth = 1e-4;
g        = 9.81;
n        = 1;   % polynomial order 
nele     = 333; % No. of elements

% Set to strucutre variable
phys.name     = caseName;
phys.minDepth = minDepth;
phys.gra      = g;
phys.n        = n;
phys.ne       = nele;

% Set initial conditions
phys = SWE_Init1d(phys);

%% Create output file
filename = ['SWE1D_', caseName, '_', num2str(nele)];
ncfile   = SWE_GenOutputFile1d(filename, phys);

%% Solve Problem
% [h, q] = SWESolverHrefinedWetDry(physics, ncfile);
[h, q]   = SWE_Solver1d(phys, ncfile);

%% Postprocess
ncfile.CloseFile;
end% func


