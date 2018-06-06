%> Abstract solver for evaluating numerical flux term
classdef SWEAbstractNumFluxSolver2d
    
    methods(Abstract)
        [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp );
    end
    
end
