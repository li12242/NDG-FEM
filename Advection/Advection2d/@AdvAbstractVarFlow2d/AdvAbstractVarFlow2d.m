classdef AdvAbstractVarFlow2d < NdgPhysMat

    properties(Constant)
        %> Number of physical field
        Nfield = 3
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
    end
    
    methods
        function obj = AdvAbstractVarFlow2d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods( Hidden )
        function [E, G] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fphys, fext )
            Ntp = mesh.cell.Np * mesh.K;
            [ fm ] = fphys( mesh.eidM ); 
            [ fp ] = fphys( mesh.eidP );
            ind = ( mesh.eidtype == int8( NdgEdgeType.Clamped ) );
            fp( ind ) = 0;
            [ um ] = fphys(mesh.eidM + Ntp); 
            [ vm ] = fphys(mesh.eidM + 2*Ntp); 
            [ uNorm ] = um .* nx + vm .* ny;
            sign_um = sign( uNorm );
            fluxS = ( fm.*( sign_um + 1 )*0.5 + fp.*( 1 - sign_um  )*0.5 ).*uNorm;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fphys )
            Ntp = mesh.cell.Np * mesh.K;
            [ fm ] = fphys( mesh.eidM ); 
            [ um ] = fphys(mesh.eidM + Ntp); 
            [ vm ] = fphys(mesh.eidM + 2*Ntp); 
            Em = fm .* um;
            Gm = fm .* vm;
            flux = Em .* nx + Gm .* ny;
        end
        
    end% methods
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( mesh, time );
    end
end

