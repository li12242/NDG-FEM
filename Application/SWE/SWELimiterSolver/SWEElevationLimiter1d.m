classdef SWEElevationLimiter1d
    %SWEELEVATIONLIMITER1D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function [ fphys ] = apply(obj, phys, fphys)
            for m = 1:phys.Nmesh % update new elevation
                fphys{m}(:,:,4) = fphys{m}(:,:,1) + fphys{m}(:,:,3);
            end
            fphys = phys.limiter.matLimit( fphys, 4 ); % enforce the elevation
            for m = 1:phys.Nmesh % update new elevation
                fphys{m}(:,:,1) = fphys{m}(:,:,4) - fphys{m}(:,:,3);
            end
            fphys = phys.limiter.matLimit( fphys, 2 );
        end
        
    end
    
end

