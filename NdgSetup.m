fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Setup environment path.\n', mfilename);
addpath( genpath( 'lib' ) );
addpath( genpath( 'NdgCell' ) );
addpath( genpath( 'NdgMesh') );
addpath( genpath( 'NdgEdge') );
addpath( genpath( 'NdgPhys') );
addpath( genpath( 'NdgUtils') );
addpath( genpath( 'NdgNetCDF') );
addpath( genpath( 'NdgLimiter') );
addpath( genpath( 'NdgPostProcess') );
addpath( genpath( 'Advection') );
addpath( genpath( 'SWE') );
fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Compile mex files.\n', mfilename);
NdgConfigure(4)
fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Finish all the setup process.\n\n', mfilename);