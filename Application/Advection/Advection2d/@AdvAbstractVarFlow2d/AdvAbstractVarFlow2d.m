classdef AdvAbstractVarFlow2d < NdgPhysMat
    
    properties (Constant)
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
    
    methods ( Hidden )
        
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
                + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            flux = Em .* nx + Gm .* ny;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, fm, fp, fext )
            ind = ( edge.ftype == enumBoundaryCondition.Clamped );
            fp(:, ind) = fext(:, ind);
        end
        
    end% methods
    
    methods( Abstract, Access = protected )
        [ fext ] = getExtFunc( obj, x, y, time );
    end
end

