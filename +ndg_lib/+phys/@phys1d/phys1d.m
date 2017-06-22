classdef phys1d < ndg_lib.phys.phys
    %PHYS1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        Nfield  % 变量个数
    end
    
    %% 误差统计
    methods
        
    end
    
    methods
        function obj = phys1d(mesh)
            obj = obj@ndg_lib.phys.phys(mesh);
        end
        
        function draw(obj, varargin)
            switch nargin
                case 1
                    f_Q = obj.f_Q(:, :, 1);
                case 2
                    f_Q = varargin{1};
            end

            if ( isempty(obj.draw_h) || ~isvalid( obj.draw_h{1} ) )
                % 若图像未绘制或窗口被关闭
                Np = obj.mesh.cell.Np; K = obj.mesh.K;
                obj.draw_h = cell(2, 1);
                list = 1:Np:(K*Np);
                g = graph();
                for n = 1:Np-1
                    g = addedge(g, list+n-1, list+n);
                end
                for fld = 1:obj.Nfield % 绘制第 fld 个物理场
                    subplot(obj.Nfield,1,fld); hold on;
                    f = f_Q(:,:,fld);
                    obj.draw_h{fld} = plot(g, ...
                        'XData', obj.mesh.x(:), 'YData', f(:), ...
                        'LineWidth', 1, ...
                        'Marker', 'o', ...
                        'NodeColor','b', ...
                        'EdgeColor', 'b', ...
                        'MarkerSize', 2);
                    box on; grid on;
                end
            else % 若图像存在
                for fld = 1:obj.Nfield
                    f = f_Q(:,:,fld); % 更新图像
                    set( obj.draw_h{fld}, 'YData', f(:) );
                end
            end
        end
    end
    
end

