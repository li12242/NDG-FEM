classdef SWEElevationLimiter
    
    methods
        
        function [ fphys ] = apply(obj, phys, fphys)
            for m = 1:phys.Nmesh % update new elevation
                fphys{m}(:,:,5) = fphys{m}(:,:,1) + fphys{m}(:,:,4);
            end
            fphys = phys.limiter.matLimit( fphys, 5 ); % enforce the elevation
            
            for m = 1:phys.Nmesh % update new elevation
                fphys{m}(:,:,1) = fphys{m}(:,:,5) - fphys{m}(:,:,4);
            end
            fphys = phys.limiter.matLimit( fphys, 2 );
            fphys = phys.limiter.matLimit( fphys, 3 );
        end
        
    end
    
end

