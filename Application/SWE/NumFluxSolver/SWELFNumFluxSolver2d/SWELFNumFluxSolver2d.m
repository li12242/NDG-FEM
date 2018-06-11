classdef SWELFNumFluxSolver2d
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp )
            fluxS = mxEvaluate( hmin, gra, nx, ny, fm, fp );
        end
    end
    
end

