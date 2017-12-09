classdef NdgAbstractVisSolver < handle
    
    properties
        phys
    end
    
    methods
        function obj = NdgAbstractVisSolver( phys )
            obj.phys = phys;
        end
    end
    
    methods( Abstract )
        evaluateViscosityRHS( obj, fphys )
    end
    
end

