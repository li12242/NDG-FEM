classdef NdgPhysMat3d < NdgPhysMat
    
    methods
        function obj = NdgPhysMat3d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods(Access = protected)
        %> @brief Evaluating the RHS term for 2d problem
        %> @details
        %> For the 2d problem, the function should call function
        %> matEvaluateRHS2d; For the 3d problem, the function
        matEvaluateRHS( obj, fphys )
    end
    
end

