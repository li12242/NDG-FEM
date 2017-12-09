classdef NdgAbstractAdvSolver < handle
    
    properties( SetAccess = protected )
        phys
    end
    
    methods
        function obj = NdgAbstractAdvSolver( phys )
            obj.phys = phys;
        end
    end
    
    methods( Abstract )
        evaluateAdvectionRHS( obj, fphys )
    end
    
end

