classdef AdvAbstractConstFlow2d < NdgPhysMat
    
    properties( Constant )
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
    end
        
    properties( Abstract, Constant )
        u0
        v0
    end
    
    methods
        function obj = AdvAbstractConstFlow2d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods( Hidden )
        
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            finalTime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.fext{m} = obj.getExtFunc(mesh.x, mesh.y, finalTime);
            end
        end
        
        function [fm, fp] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            fm = fphys( mesh.eidM );
            fp = fphys( mesh.eidP );
            ind = ( mesh.eidtype == int8(NdgEdgeType.Clamped) );
            fp(ind) = 0;
        end
        
        function [E, G] = matEvaluateFlux( obj, mesh, fieldValue )
            E = obj.u0 .* fieldValue;
            G = obj.v0 .* fieldValue;
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp )
            [ uNorm ] = obj.u0.* nx + obj.v0.* ny;
            sign_um = sign( uNorm );
            fluxS = ( fm.*( sign_um + 1 )*0.5 + fp.*( 1 - sign_um  )*0.5 ).*uNorm;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm )
            Em = fm .* obj.u0;
            Gm = fm .* obj.v0;
            flux = Em .* nx + Gm .* ny;
        end
        
        function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext )
            ind = ( edge.ftype == enumBoundaryCondition.Clamped );
            fP(:, ind) = 0;
        end
        
    end
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( obj, x, y, time );
    end
    
end

