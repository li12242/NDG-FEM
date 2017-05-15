classdef BJ_line < ndg_utility.limiter.BJ.BJ
    %BJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = BJ_line(mesh, cell)
            obj = obj@ndg_utility.limiter.BJ.BJ(mesh, cell);
        end
        
        function f_limit = limit(obj, f_Q)
            f_limit = f_Q;
            % Compute cell averages, v
            c_mean = obj.mesh.cell_mean(f_Q);            
            % find max and min cell averages
            c_max = max([c_mean(obj.mesh.EToE); c_mean]);
            c_min = min([c_mean(obj.mesh.EToE); c_mean]);
            % find elements in need of limiting
            f_min = min(f_Q); f_max = max(f_Q);
            K = 1;
            ids = ((f_max > c_max | f_min < c_min) ...
                & ((f_max - f_min) > K.*obj.mesh.vol));
            % Check to see if any elements require limiting
            if(any(ids))
              % correction factor
              a = ones(size(f_Q));
              f_dive = 1./bsxfun(@minus, f_Q, c_mean);
              ind = bsxfun(@gt, f_Q, c_mean);
              tmp = min(1, bsxfun(@times, c_max-c_mean, f_dive));
              a(ind) = tmp(ind);
              
              ind = bsxfun(@lt, f_Q, c_mean);
              tmp = min(1, bsxfun(@times, c_min-c_mean, f_dive));
              a(ind) = tmp(ind);
              % apply slope limiter to selected elements
              tmp = bsxfun(@plus, c_mean, a./f_dive);
              f_limit(:,ids) = tmp(:, ids);
            end% if
        end
    end
    
end

function a = limitCoeff(maxu, minu, meanu, u)
a = ones(size(u)); Np = size(u, 1);

ind = u > meanu;
a(ind) = min(1, ( maxu(ind) - meanu(ind) )./( u(ind) - meanu(ind) ) );

ind = u < meanu;
a(ind) = min(1, ( minu(ind) - meanu(ind) )./( u(ind) - meanu(ind) ) );

amin = min(a);
a = ones(Np, 1)*amin;
end

