classdef NdgLimiterType < int8
    
    enumeration
        %> A None limiter will doing nothing in the limiting function
        None  (0)
        %> vertex-based limiter from Li (2017, Computers and Fluids)
        VertexLimiter (1)
    end
    
end

