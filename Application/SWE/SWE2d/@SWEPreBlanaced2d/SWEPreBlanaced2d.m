classdef SWEPreBlanaced2d < SWEConventional2d
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.status, fphys );
        end
    end
    
    methods( Access = protected )
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).status( ~wetflag ) = int8( enumSWERegion.Dry );
                obj.meshUnion(m).status(  wetflag ) = int8( enumSWERegion.Wet );
            end
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography2d...
                    ( obj.gra, mesh.status, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

