classdef SWEConventional1d < SWEAbstract1d
    
    properties( Constant )
        %> Number of physical field - {h, hu, z}
        Nfield = 3
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
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography1d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

