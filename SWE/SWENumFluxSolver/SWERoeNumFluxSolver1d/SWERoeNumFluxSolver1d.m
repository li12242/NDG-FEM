classdef SWERoeNumFluxSolver1d < SWEAbstractNumFluxSolver1d
    methods 
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, fm, fp )
            fluxS = mxEvaluate( hmin, gra, nx, fm, fp );
        end
    end
end 