classdef SWEDepthLimiter1d
    %SWEDEPTHLIMITER1D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function fphys = apply( obj, phys, fphys )
            fphys = phys.limiter.matLimit( fphys, 1 );
            fphys = phys.limiter.matLimit( fphys, 2 );
        end
    end
    
end

