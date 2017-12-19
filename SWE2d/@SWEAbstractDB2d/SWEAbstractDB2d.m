classdef SWEAbstractDB2d < SWEAbstractCB2d

    properties(Abstract, Constant)
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
        %> Physical field
        Nfield
    end
    
    methods
        function obj = SWEAbstractDB2d()
            obj = obj@SWEAbstractCB2d();
        end
    end
    
    methods( Hidden )
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fphys, fext )
            [ fluxS ] = mxEvaluateSurfNumFlux2d( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, mesh.eidtype, nx, ny, fext, fphys );
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fphys )
            [ flux ] = mxEvaluateSurfFlux2d( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, mesh.eidtype, nx, ny, fphys);
        end
    end
    
    methods( Access = protected )
    end
end

