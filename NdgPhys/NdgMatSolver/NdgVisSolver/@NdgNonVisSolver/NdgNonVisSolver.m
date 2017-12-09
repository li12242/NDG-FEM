classdef NdgNonVisSolver < NdgAbstractVisSolver
    
    methods
        function obj = NdgNonVisSolver( phys )
            obj = obj@NdgAbstractVisSolver( phys );
        end
        
        function evaluateViscosityRHS( obj, fphys )
        end% func
    end
    
end

