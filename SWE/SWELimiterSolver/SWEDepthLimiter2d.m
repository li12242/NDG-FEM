classdef SWEDepthLimiter2d < handle
    
    methods
        
        function fphys = apply( obj, phys, fphys )
            fphys = phys.limiter.matLimit( fphys, 1 );
            fphys = phys.limiter.matLimit( fphys, 2 );
            fphys = phys.limiter.matLimit( fphys, 3 );
        end
    end% methods
    
end

