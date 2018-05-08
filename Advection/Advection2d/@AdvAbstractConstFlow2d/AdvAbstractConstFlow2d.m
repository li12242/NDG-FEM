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
                obj.fext{m} = obj.getExtFunc(obj.meshUnion, finalTime);
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
        
        function [ frhs ] = matEvaluateEdgeRHS( obj, edge, fluxM, fluxP, fluxS )
            frhs = zeros( edge.mesh.cell.Np, edge.mesh.K, obj.Nvar );
            for fld = 1:obj.Nvar
                for k = 1:edge.Ne
                    e1 = edge.FToE(1, k);
                    e2 = edge.FToE(2, k);
                    n1 = edge.FToN1(:, k);
                    n2 = edge.FToN2(:, k);
                    
                    deltaFlux1 = fluxM(:,k,fld) - fluxS(:,k,fld);
                    frhs(n1, e1, fld) = frhs(n1, e1, fld) + ...
                        edge.eCell.M * ( edge.Js(:, k) .* deltaFlux1 );
                    %                     frhs(:, e1, fld) = ...
                    %                         frhs(:, e1, fld) + edge.mesh.cell.invM(:, n1) * ...
                    %                         ( edge.bcell.M * ( edge.Js(:, k) .* deltaFlux1 ) ) ...
                    %                         ./ edge.mesh.J( :, e1 );
                    
                    deltaFlux2 = fluxP(:,k,fld) - fluxS(:,k,fld);
                    frhs(n2, e2, fld) = frhs(n2, e2, fld) - ...
                        edge.eCell.M * ( edge.Js(:, k) .* deltaFlux2 );
                    %                     frhs(:, e2, fld) = ...
                    %                         frhs(:, e2, fld) - edge.mesh.cell.invM(:, n2) * ...
                    %                         ( edge.bcell.M * ( edge.Js(:, k) .* deltaFlux2 ) ) ...
                    %                         ./ edge.mesh.J( :, e2 );
                end
                frhs(:, :, fld) = edge.mesh.cell.invM * frhs(:, :, fld) ./ edge.mesh.J;
            end            
        end
    end
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( mesh, time );
    end
    
end

