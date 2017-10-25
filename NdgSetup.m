fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Setup environment path.\n', mfilename);
addpath( genpath( 'lib' ) );
addpath( genpath( 'NdgCell' ) );
addpath( genpath( 'NdgMesh') );
addpath( genpath( 'NdgPhys') );
addpath( genpath( 'NdgSolver') );
addpath( genpath( 'NdgUtils') );
addpath( genpath( 'Advection') );
fprintf('----------------------------------------------------------\n')

fprintf('%s:: Compile mex files.\n', mfilename);
NdgConfigure(2)
fprintf('----------------------------------------------------------\n')
fprintf('%s:: Finish all the setup process.\n\n', mfilename);
