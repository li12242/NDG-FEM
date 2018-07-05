classdef NdgNonVisSolver < NdgAbstractVisSolver
    
    methods
        function obj = NdgNonVisSolver( phys )
            obj = obj@NdgAbstractVisSolver( phys, 0, 0 );
        end
        
        function matEvaluateRHS( obj, fphys, frhs )
        end% func
    end
    
end

