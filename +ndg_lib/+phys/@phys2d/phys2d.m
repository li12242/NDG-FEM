classdef phys2d < ndg_lib.phys.phys
    %PHYS2D Summary of this class goes here
    %   Detailed explanation goes here
    
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
        % DRAW draw the three dimensional figure for the physical field
        %
            if ( numel(fld) == 1 )
                f = obj.f_Q(:,:,fld);
            else
                f = fld;
            end
            if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
                % if the figure is not created or it is closed
                EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
                EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
                    *ones(1, obj.mesh.cell.Nfptotal);
                obj.draw_h = patch(...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)], ...
                    'Faces', EToV, ...
                    'FaceColor', 'interp', ...
                    'FaceVertexCData', f(:));
            else % if the figure exists
                set(obj.draw_h, ...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)],...
                    'FaceVertexCData', f(:));
            end
        end
        
        function surf(obj, fld)
        % SURF draw the two dimensional figure for the physical field  
            f = obj.f_Q(:,:,fld);
            if ( isempty(obj.surf_h) || ~isvalid(obj.surf_h) ) 
                % if the figure is not created or it is closed
                EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
                EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
                    *ones(1, obj.mesh.cell.Nfptotal);
                obj.surf_h = patch(...
                    'Vertices', [obj.mesh.x(:), obj.mesh.y(:)], ...
                    'Faces', EToV, ...
                    'FaceColor', 'interp', ...
                    'FaceVertexCData', f(:));
            else % if the figure exists
                set(obj.surf_h, 'FaceVertexCData', f(:));
            end
        end
    end
    
end

