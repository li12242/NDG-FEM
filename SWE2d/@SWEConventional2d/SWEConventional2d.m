classdef SWEConventional2d < SWEAbstract2d
    %SWEABSTRACTCONVENTIONALDISCONTINUOUSBOTTOM2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %> Variable field - {h, hu, hv, b}
        Nfield = 4
        %> Variable field - {h, hu, hv}
        Nvar = 3
        %> field index of variable field
        varFieldIndex = [ 1,2,3 ]
    end
    
    methods
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
    end
    
    methods( Access = protected )
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography2d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

