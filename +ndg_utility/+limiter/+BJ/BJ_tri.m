classdef BJ_tri < ndg_utility.limiter.BJ.BJ
    %BJ_LINE Summary of this class goes here
    %   Detailed explanation goes here
    properties
    end% func
    
    methods
        function [ obj ] = BJ_tri(mesh)
            obj = obj@ndg_utility.limiter.BJ.BJ(mesh);
        end
        
        function [ f_Q ] = limit( obj, f_Q )
            [] = 
        end        
    end
    
end% classdef