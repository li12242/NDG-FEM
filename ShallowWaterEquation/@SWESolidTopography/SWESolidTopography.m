classdef SWESolidTopography < NdgPhysMat
    %SWESOLIDTOPOGRAPHY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        hmin
    end
    
    properties(Constant)
        %> Physical field - {h, hu, hv, b, bx, by}
        Nfield = 6
        %> Variable field - {h, hu, hv}
        Nvar = 3
        %> field index of variable field
        varFieldIndex = [1,2,3]
        %> gravity acceleration
        gra = 9.81
    end
    
    methods
        function obj = SWESolidTopography()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods( Access = protected )
         
        [ fphys ] = matEvaluateLimiter( obj, fphys )
        [ fphys ] = matEvaluatePostFunc(obj, fphys)
        [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
        [ dflux ] = matEvaluateNumericalFlux( obj, mesh, fphys, fext )
        [ dt ] = matUpdateTimeInterval( obj, fphys )
                
        function matEvaluateRHS( obj, fphys )
            obj.matEvaluateRHS2d( fphys );
        end
        
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin ); 
                obj.meshUnion(m).EToR( ~wetflag ) = int8( NdgRegionType.Dry );
                obj.meshUnion(m).EToR(  wetflag ) = int8( NdgRegionType.Wet );
            end
        end
    end
    
end

