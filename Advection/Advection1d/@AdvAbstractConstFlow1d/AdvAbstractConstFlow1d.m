classdef AdvAbstractConstFlow1d < NdgPhysMat
    
    properties( Abstract )
        u0
    end
    
    properties( Constant )
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
    end
    
    methods( Hidden )
        
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            finalTime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                obj.fext{m} = obj.getExtFunc(obj.meshUnion, finalTime);
            end
        end
        
        function [fm, fp] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            fm = fphys( mesh.eidM );
            fp = fphys( mesh.eidP );
            fe = fext( mesh.eidM );
            ind = ( mesh.eidtype == int8(NdgEdgeType.Clamped) );
            fp(ind) = fe(ind);
        end
        
        function [ E ] = matEvaluateFlux( obj, mesh, fphys )
            E = obj.u0 .* fphys;
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, fm, fp )
            [ uNorm ] = obj.u0.* nx;
            sign_um = sign( uNorm );
            fluxS = ( fm.*( sign_um + 1 )*0.5 + fp.*( 1 - sign_um  )*0.5 ).*uNorm;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, fm )
            Em = fm .* obj.u0;
            flux = Em .* nx;
        end
    end
    
    methods( Access = protected )
        
        function matUpdateExternalField( obj, time, fphys )
            matUpdateExternalField@NdgPhysMat( obj, time, fphys );
            for m = 1:obj.Nmesh
                obj.fext{m} = obj.getExtFunc(obj.meshUnion, time);
            end
        end
    end
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( mesh, time );
    end
    
end

