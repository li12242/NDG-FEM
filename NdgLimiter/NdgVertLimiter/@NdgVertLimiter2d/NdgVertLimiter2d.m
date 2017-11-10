classdef NdgVertLimiter2d < NdgVertLimiter
    %NDGVERTLIMITER2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = NdgVertLimiter2d( mesh )
            obj = obj@NdgVertLimiter( mesh );
        end
        
        [ fphys ] = matVertLimit( obj, fphys, fieldId );
    end
    
end

