classdef NdgNonAdvectionSolver < NdgAbstractAdvectionSolver
    
    methods
        function obj = NdgNonAdvectionSolver( phys )
            obj = obj@NdgAbstractAdvectionSolver( phys );
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            % do nothing ...
        end
    end
    
end

