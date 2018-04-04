classdef SWEPreBlanaced1d < SWEConventional1d
    
    methods( Hidden )
        function [ E ] = matEvaluateFlux( obj, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
    end
    
    methods( Access = protected )
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography1d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

