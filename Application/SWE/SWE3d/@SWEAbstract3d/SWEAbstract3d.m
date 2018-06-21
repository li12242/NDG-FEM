classdef SWEAbstract3d < handle
    %SWEABSTRACT3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
    end
    
    properties ( Constant )
        %> number of physical field
        Nfield = 6
        %> number of variable field
        Nvar = 3
        %> index of variable in physical field
        varFieldIndex = [ 1, 2, 3 ]
    end
    
    properties ( SetAccess = protected )
        %> gradient of bottom elevation
        zGrad
    end
    
    properties ( SetAccess = private )
        %> solver for coriolis source term
        coriolisSolver
        %> solver for friction source term
        frictionSolver
        %> solver for wind source term
        windSolver
        %> solver for unmerical flux
        numfluxSolver
        %> limiter type
        limiterSolver
    end
    
    methods
    end
    
end

