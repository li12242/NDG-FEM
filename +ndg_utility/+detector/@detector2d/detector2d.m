classdef detector2d < ndg_utility.detector.detector
    %DETECTOR2D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = detector2d(mesh, xd, yd, dt, ftime, Nfield)
            zd = zeros(size(xd));
            obj = obj@ndg_utility.detector.detector(mesh, xd, yd, zd, ...
                dt, ftime, Nfield);
        end
        
        function [kd, rd, sd, td] = findlocate(obj)
            % 初始化
            td = ones(obj.Nd, 1);
            % 寻找检测点所在单元
            kd = find_loc_cell(obj.mesh.x, obj.mesh.y, ...
                obj.mesh.cell.Fmask(1,:), obj.xd, obj.yd)';
            
            switch obj.mesh.cell.type
                case ndg_lib.std_cell_type.Tri
                    [ rd, sd ] = loc_tri_coor(obj, kd, obj.xd, obj.yd);
                case ndg_lib.std_cell_type.Quad
                    [ rd, sd ] = loc_quad_coor(obj, kd, obj.xd, obj.yd);
            end
        end% func
    end
    
    methods(Hidden, Access=protected)
        [ rd, sd ] = loc_tri_coor(obj, kd, xd, yd);
        [ rd, sd ] = loc_quad_coor(obj, kd, xd, yd);
    end
    
end

