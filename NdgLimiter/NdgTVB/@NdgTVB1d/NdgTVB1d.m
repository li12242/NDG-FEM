classdef NdgTVB1d < NdgTVBAbstract
    %NDGTVB1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        h2
    end
    
    methods
        function obj = NdgTVB1d( mesh, M )
            obj = obj@NdgTVBAbstract( mesh, M );
            
            obj.h2 = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                obj.h2{m} = obj.meshUnion(m).LAV.^2 .* M;
            end
            
        end
        
        function fphys = matLimit( obj, fphys, fldId )
            [ hc, hca ] = obj.matEvaluateIntegralMean( fphys, fldId );
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
%                 vol1 = sum( hc{m}.*mesh.LAV );
%                 vol2 = sum( mesh.GetMeshIntegralValue( fphys{m}(:,:,fldId) ) );
                
                dh = fphys{m}( end, :, fldId ) ...
                    - fphys{m}( 1, :, fldId );
                dhl =   ( hc{m} - hca{m}(1, :) );
                dhr = - ( hc{m} - hca{m}(2, :) );
                ids = ( abs( dh ) > obj.h2{m} );
                
                if ( any( ids ) )
                    tmp = obj.minmod([dh; dhl; dhr]);
                    slope = tmp./mesh.LAV;
%                     t_Q = bsxfun(@plus, hc{m}, ... % reconstruct limited node values
%                         bsxfun(@times, slope, ...
%                         bsxfun(@minus, mesh.x, mesh.xc)));
                    t_Q = hc{m} + slope .* ( mesh.x - mesh.xc );
                    deltaAve = hc{m} - mesh.GetMeshAverageValue( t_Q );
                    t_Q = t_Q + deltaAve;
                    fphys{m}(:, ids, fldId) = t_Q(:, ids);
                end
                
%                 vol3 = sum( mesh.GetMeshIntegralValue( fphys{m}(:,:,fldId) ) );
%                 fprintf('vol1 = %e, vol2 = %e\n', vol2 - vol1, vol3 - vol1);
            end
        end
    end
    
end

