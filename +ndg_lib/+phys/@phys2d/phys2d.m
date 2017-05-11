classdef phys2d < ndg_lib.phys.phys
    %PHYS2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        Nfield  % 变量个数
    end
    properties(Hidden=true)
        surf_h
    end
    
    methods
        function obj = phys2d(mesh)
            obj = obj@ndg_lib.phys.phys(mesh);
        end
    end
    
    methods
        function draw(obj, fld)
            f = obj.f_Q(:,:,fld);
            if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
                % 若图像未绘制或窗口被关闭
                EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
                EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
                    *ones(1, obj.mesh.cell.Nfptotal);
                obj.draw_h = patch(...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)], ...
                    'Faces', EToV, ...
                    'FaceColor', 'interp', ...
                    'FaceVertexCData', f(:));
            else % 若图像存在
                set(obj.draw_h, ...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)],...
                    'FaceVertexCData', f(:));
            end
        end
        
        function surf(obj, fld)
            f = obj.f_Q(:,:,fld);
            if ( isempty(obj.surf_h) || ~isvalid(obj.surf_h) ) 
                % 若图像未绘制或窗口被关闭
                EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
                EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
                    *ones(1, obj.mesh.cell.Nfptotal);
                obj.surf_h = patch(...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:)], ...
                    'Faces', EToV, ...
                    'FaceColor', 'interp', ...
                    'FaceVertexCData', f(:));
            else % 若图像存在
                set(obj.surf_h, 'FaceVertexCData', f(:));
            end
        end
    end
    
end

