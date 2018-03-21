classdef SWEPreBlanaced1d < SWEAbstract1d
    
    properties( Constant )
        %> Number of physical field - {h, hu, z, eta}
        Nfield = 4
        %> Number of variable field - {h, hu}
        Nvar = 2
        %> 
        varFieldIndex = [1, 2]
    end
    
    methods
        function [ E ] = matEvaluateFlux( obj, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
    end
    
    methods( Access = protected )
        function fphys = matEvaluateLimiter( obj, fphys )
            fphys = obj.limiter.matLimit( fphys, 2 );
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,4) = fphys{m}(:,:,1) + fphys{m}(:,:,3);
            end
%             vol1 = sum( obj.meshUnion.GetMeshIntegralValue( fphys{1}(:, :, 1) ) );
            fphys = obj.limiter.matLimit( fphys, 4 ); % enforce the elevation
            for m = 1:obj.Nmesh % update new elevation
                fphys{m}(:,:,1) = fphys{m}(:,:,4) - fphys{m}(:,:,3);
            end
%             vol2 = sum( obj.meshUnion.GetMeshIntegralValue( fphys{1}(:, :, 1) ) );
%             fprintf('vol1 = %e\n', vol2 - vol1)
        end

        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography1d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

