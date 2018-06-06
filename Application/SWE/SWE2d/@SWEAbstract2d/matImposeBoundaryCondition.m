function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext )

% impose boundary condition
[ fP ] = mxImposeBoundaryCondition( obj.gra, nx, ny, fP, fext, edge.ftype );

% hydrostatic reconstruction
[ fM, fP ] = mxHydrostaticReconstruction( obj.hmin, fM, fP );

end

