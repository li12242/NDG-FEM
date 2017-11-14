classdef NdgNonLimiter < NdgAbstractLimiter
    
    methods
        function obj = NdgNonLimiter( mesh )
            obj = obj@NdgAbstractLimiter( mesh );
        end
        
        function fphys = matLimit( obj, fphys, fieldId )
            % doing nothing
        end
    end
    
end

