classdef NdgVertLimiter2d < NdgVertLimiter
    
    properties
    end
    
    methods
        function obj = NdgVertLimiter2d( mesh )
            obj = obj@NdgVertLimiter( mesh );
        end
        
        [ fphys ] = matLimit( obj, fphys, fieldId );
    end
    
end

