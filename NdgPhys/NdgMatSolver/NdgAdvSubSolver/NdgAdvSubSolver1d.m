classdef NdgAdvSubSolver1d < NdgAbstractAdvSubSolver
    
    properties
        
    end
    
    methods
        function obj = NdgAdvSubSolver1d( phys, regionId )
            obj = obj@NdgAbstractAdvSubSolver( phys, regionId );
        end
        
        function evaluateAdvectionRHS
        end
    end
    
end

