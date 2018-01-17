classdef NdgBJ1d < NdgBJAbstract
    
    properties
    end
    
    methods
        function obj = NdgBJ1d( mesh )
            obj = obj@NdgBJAbstract( mesh );
        end
        
        function fphys = matLimit( obj, fphys, fldId )
            [ cmean, cmax, cmin ] = obj.accessCellExtrema( fphys, fldId );
            
            for m = 1:obj.Nmesh
                temp = fphys{m}(:, :, fldId);
                fmin = min( temp );
                fmax = max( temp );
                
                ids = ( ( fmax > cmax{m} | fmin < cmin{m} ) ); %...
                %& ( ( fmax - fmin ) > K.*obj.mesh.LAV ) );
                if(any(ids))
                    a = ones(size( temp ));
                    f_dive = 1./bsxfun(@minus, temp, cmean{m});
                    ind1 = bsxfun(@gt, temp, cmax{m}); % correction for f_Q > cmean
                    ind2 = bsxfun(@lt, temp, cmin{m}); % correction for f_Q < cmean
                    tmp1 = min(1, bsxfun(@times, cmax{m}-cmean{m}, f_dive));
                    tmp2 = min(1, bsxfun(@times, cmin{m}-cmean{m}, f_dive));
                    a(ind1) = tmp1(ind1);
                    a(ind2) = tmp2(ind2);
                    a = min( a );
                    % apply slope limiter to selected elements
                    temp = bsxfun(@plus, cmean{m}, bsxfun(@times, a, 1./f_dive));
                    fphys{m}(:,ids, fldId) = temp(:, ids);
                end
            end
        end% func
    end
    
end

