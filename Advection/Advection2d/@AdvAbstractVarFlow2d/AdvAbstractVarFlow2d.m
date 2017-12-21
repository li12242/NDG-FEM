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
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            finalTime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                obj.fext{m} = obj.getExtFunc(obj.meshUnion, finalTime);
            end
        end
        
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
        end
        
        function [ fm, fp ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            fm(:,:,1) = fphys( mesh.eidM );
            fp(:,:,1) = fphys( mesh.eidP );
            ind = ( mesh.eidtype == int8(NdgEdgeType.Clamped) );
            fp(ind) = 0;
            Ntp = mesh.cell.Np * mesh.K;
            fm(:,:,2) = fphys(mesh.eidM + Ntp); 
            fm(:,:,3) = fphys(mesh.eidM + 2*Ntp); 
            fp(:,:,2) = fphys(mesh.eidP + Ntp); 
            fp(:,:,3) = fphys(mesh.eidP + 2*Ntp); 
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1).*( sign_um + 1 )*0.5 + fp(:,:,1).*( 1 - sign_um  )*0.5 ).*uNorm;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            flux = Em .* nx + Gm .* ny;
        end
        
    end% methods
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( mesh, time );
    end
end

