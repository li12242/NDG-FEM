classdef TVB_line < ndg_utility.limiter.TVB.TVB
    %LIMITER1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx % diversity of x coordinate in each cell
    end
    
    methods
        function obj = TVB_line(mesh,cell)
            obj = obj@ndg_utility.limiter.TVB.TVB(mesh, cell);
            
            if (cell.type ~= ndg_lib.std_cell_type.Line)
                error('The input cell type is incorrect!');
            end
            
            x_b  = obj.mesh.x(obj.cell.Fmask,:);
            obj.dx = x_b(1, :) - x_b(2, :);
        end
        
        function f_limt = limit(obj, f_Q, M)
            % compute cell averages, f_mean
            f_mean = obj.mesh.cell_mean(f_Q);

            % Apply slope limiter as needed
            f_limt = f_Q;

            % find end values of each element
            f_b  = f_Q(obj.cell.Fmask,:);
            df_b = f_b(1,:) - f_b(2,:);
            
            % find cell averages
            f_mean_adj = f_mean(obj.mesh.EToE);
            df_mean1 = f_mean_adj(1, :) - f_mean;
            df_mean2 = f_mean - f_mean_adj(2, :);

            % Apply reconstruction to find elements in need of limiting
            ids = ( abs(df_b) > obj.h2*M );

            % Check to see if any elements require limiting, and apply 
            % slope limiter to selected elements
            if any(ids)
                % calculate the limited boundary diversities
                tmp = ndg_utility.limiter.minmod([df_b; df_mean1; df_mean2]);
                slope = tmp./obj.dx; % the limited slopes
                t_Q = bsxfun(@plus, f_mean, ... % reconstruct limited node values
                  bsxfun(@times, slope, bsxfun(@minus, obj.mesh.x, obj.xc)));
                f_limt(:, ids) = t_Q(:, ids);
            end% if
        end
    end
    
end

