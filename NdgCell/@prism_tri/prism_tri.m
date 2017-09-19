classdef prism_tri < ndg_lib.std_cell.std_cell
    %PRISM_TRI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        type = ndg_lib.std_cell_type.Prism_Tri % 单元类型
        Nv = 6
        vol = 4
        vr = [-1,  1, -1, -1,  1, -1]'
        vs = [-1, -1,  1, -1, -1,  1]'
        vt = [-1, -1, -1,  1,  1,  1]'
        Nfv = [3, 4, 4, 4, 3];
        FToV = [1,2,3,0; 1,2,5,4; 2,3,6,5; 3,1,4,6; 4,5,6,0]';
        Nface = 5
        faceType = [...
            ndg_lib.std_cell_type.Tri, ...
            ndg_lib.std_cell_type.Quad, ...
            ndg_lib.std_cell_type.Quad, ...
            ndg_lib.std_cell_type.Quad, ...
            ndg_lib.std_cell_type.Tri]
    end
    
    methods(Access=protected)
        [ Np,r,s,t ] = node_coor_func(obj, N);
        [ dr, ds, dt ] = derivative_orthogonal_func(obj, N, ind, r, s, t);
    end
    
    methods
        function obj = prism_tri(N)
            obj = obj@ndg_lib.std_cell.std_cell(N);
        end
        
        [ f ] = orthogonal_func(obj, N, ind, r, s, t);

        function node_val = project_vert2node(obj, vert_val)
            node_val = ( -(obj.r+obj.s)*vert_val(1, :) ...
                + (1+obj.r)*vert_val(2, :)...
                + (1+obj.s)*vert_val(3, :) ).*(1 - obj.t)*0.25 ...
                + ( -(obj.r+obj.s)*vert_val(4, :) ...
                + (1+obj.r)*vert_val(5, :)...
                + (1+obj.s)*vert_val(6, :) ).*(1 + obj.t)*0.25 ;
        end
        
        function draw(obj, varargin)
            figure('color', 'w');
            subplot(1, 2, 1); % 标准单元顶点分布
            v1 = [1, 2, 3, 1, 2, 3, 4, 5, 6];
            v2 = [2, 3, 1, 4, 5, 6, 5, 6, 4];
            g = graph(v1, v2);
            
            plot(g, 'XData', obj.vr(:), ...
                'YData', obj.vs(:), ...
                'ZData', obj.vt(:),...
                'LineWidth', 1, ...
                'Marker', 'o', ...
                'MarkerSize', 16, ...
                'NodeColor','r', ...
                'EdgeColor', 'k', ...
                'MarkerSize', 2);
            box on; grid on; axis equal;
            title('Vertex distribution')
            
            subplot(1, 2, 2); % 标准单元节点分布
            hold on;
            plot3(obj.r(:), obj.s(:), obj.t(:), 'ro', ...
                'MarkerSize', 4, ...
                'MarkerFaceColor', 'r')
            tri = ndg_lib.std_cell.tri(obj.N);
            for n = 1:(obj.N+1)
                ind = tri.Fmask + tri.Np*(n-1);
                plot3(obj.r(ind), obj.s(ind), obj.t(ind), 'b-', ...
                    'LineWidth', 1.5);
            end
            for n = 1:tri.Np
                ind = n:tri.Np:obj.Np;
                plot3(obj.r(ind), obj.s(ind), obj.t(ind), 'k-.', ...
                    'LineWidth', 1.5);
            end
            view([-24, 35]);
            box on; grid on; axis equal;
            lim = [-1.2, 1.2];
            xlim(lim); ylim(lim); zlim(lim);
            legend({'Nodes'}, 'box', 'off', ...
                'Position', [.8, .69, .132, .054])
            title('Node distribution')
        end
    end
    
end

