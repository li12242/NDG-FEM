fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Clear environment path.\n', mfilename);
rmpath( genpath( 'lib' ) );
rmpath( genpath( 'NdgCell' ) );
rmpath( genpath( 'NdgMesh') );
rmpath( genpath( 'NdgEdge') );
rmpath( genpath( 'NdgPhys') );
rmpath( genpath( 'Utilities') );
rmpath( genpath( 'NdgNetCDF') );
rmpath( genpath( 'NdgLimiter') );
rmpath( genpath( 'PostProcess') );
rmpath( genpath( 'Advection') );
rmpath( genpath( 'SWE') );
fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Clear mex files.\n', mfilename);
NdgConfigure -clear
fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Finish all the setup process.\n\n', mfilename);