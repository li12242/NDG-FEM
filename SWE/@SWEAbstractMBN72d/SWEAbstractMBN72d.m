classdef SWEAbstractMBN72d < SWEAbstractDBN62d
    
    properties
        %> Physical field - {h, hu, hv, b, bx, by, eta}
        Nfield = 9
    end
    
    methods
        function obj = SWEAbstractMBN72d( )
            obj = obj@SWEAbstractDBN62d( );
        end
    end
    
    methods(Access=protected)
        function fphys = matEvaluateLimiter( obj, fphys )
            fphys = obj.limiter.matLimit( fphys, 1 );
            fphys = obj.limiter.matLimit( fphys, 2 );
            fphys = obj.limiter.matLimit( fphys, 3 );
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,7) = fphys{m}(:,:,1) + fphys{m}(:,:,4);
            end
            fphys = obj.limiter.matLimit( fphys, 7 ); % enforce the elevation
            for m = 1:obj.Nmesh % set the new 
                fphys{m}(:,:,4) = fphys{m}(:,:,1) - fphys{m}(:,:,1);
            end
        end
    end
    
end

