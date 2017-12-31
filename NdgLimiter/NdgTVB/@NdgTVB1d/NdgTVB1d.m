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
                obj.h2{m} = obj.meshUnion(m).LAV.^2;
            end
            
        end
        
        function fphys = matLimit( obj, fphys, fldId )
            [ hc, hca ] = obj.matEvaluateIntegralMean( fphys, fldId );
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
                dh = fphys{m}( mesh.cell.Fmask(1), :, fldId ) ...
                    - fphys{m}( mesh.cell.Fmask(2), :, fldId );
                dhl = hca{m}(1, :) - hc{m};
                dhr = hc{m} - hca{m}(2, :);
                ids = ( abs( dh ) > obj.h2{m} );
                
                if ( any( ids ) )
                    tmp = obj.minmod([dh; dhl; dhr]);
                    slope = tmp./mesh.LAV;
                    t_Q = bsxfun(@plus, hc{m}, ... % reconstruct limited node values
                        bsxfun(@times, slope, ...
                        bsxfun(@minus, mesh.x, mesh.xc)));
                    fphys{m}(:, ids, fldId) = t_Q(:, ids);
                end
            end
        end
    end
    
end

