classdef SWEPreBlanaced2d < SWEConventional2d
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
    end
    
    methods( Access = protected )
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = any( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).EToR( ~wetflag ) = int8( NdgRegionType.Dry );
                obj.meshUnion(m).EToR(  wetflag ) = int8( NdgRegionType.Wet );
            end
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography2d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
        
%         function [ fphys ] = matEvaluatePostFunc(obj, fphys)
%             for m = 1:obj.Nmesh
%                 hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
%                 qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
%                 qyc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,3) );
%                 fphys{m}(:,:,1:3) = mxEvaluatePostFunc2d( obj.hmin, fphys{m}, hc, qxc, qyc );
%             end
%             obj.matUpdateWetDryState( fphys );
%         end
        
    end
    
end

