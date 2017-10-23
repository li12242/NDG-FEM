fprintf('\n----------------------------------------------------------\n')
fprintf('setup environment path.\n');
addpath( genpath( 'lib' ) );
addpath( genpath( 'NdgCell' ) );
addpath( genpath( 'NdgMesh') );
addpath( genpath( 'NdgPhys') );
addpath( genpath( 'NdgSolver') );
addpath( genpath( 'NdgUtils') );
addpath( genpath( 'Advection') );
fprintf('----------------------------------------------------------\n')

fprintf('compile mex files:\n');
configure(2)
fprintf('----------------------------------------------------------\n')
fprintf('finish all the setup process.\n\n');
