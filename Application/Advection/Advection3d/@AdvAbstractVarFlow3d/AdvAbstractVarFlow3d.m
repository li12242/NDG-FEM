classdef AdvAbstractVarFlow3d < NdgPhysMat
    %ADVABSTRACTVARFLOW3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %> Number of physical field
        Nfield = 4
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
    end
    
    methods
        function obj = AdvAbstractVarFlow3d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods ( Abstract, Access = protected )
        f_ext = getExtFunc( obj, mesh, time );
    end
    
    methods ( Hidden )
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            finalTime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                obj.fext{m} = obj.getExtFunc( obj.meshUnion(m), finalTime );
            end
        end
        
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
            H = fphys(:,:,4) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, nz, fm, fp )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny + fm(:,:,4) .* nz;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
                + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, nz, fm )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            Hm = fm(:,:,1) .* fm(:,:,4);
            flux = Em .* nx + Gm .* ny + Hm .* nz;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, fm, fp, fext )
            ind = ( edge.ftype == 5 );
            fp(:, ind) = 0;
        end
    end
    
end

