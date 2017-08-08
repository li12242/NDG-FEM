classdef quad < ndg_lib.std_cell.std_cell
    %QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    %% 基本属性
    properties(Constant)
        type = ndg_lib.std_cell_type.Quad   % 单元类型
        Nv = 4
        vol = 4
        vr = [-1,  1,  1, -1]'
        vs = [-1, -1,  1,  1]'
        vt = [ 0,  0,  0,  0]'
        Nfv = [2,2,2,2]'
        FToV = [1,2; 2,3; 3,4; 4,1]'
        Nface = 4
        faceType = [...
            ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Line]
    end
    
    methods(Access=protected)
        [Np, r,s,t] = node_coor_func(obj, N);
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
    end
    
    methods
        function obj = quad(N)
            obj = obj@ndg_lib.std_cell.std_cell(N);
        end
        f = orthogonal_func(obj, N, ind, r, s, t);
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.25*(...
                (1-obj.r).*(1-obj.s)*vert_val(1,:) + ...
                (1+obj.r).*(1-obj.s)*vert_val(2,:) + ...
                (1+obj.r).*(1+obj.s)*vert_val(3,:) + ...
                (1-obj.r).*(1+obj.s)*vert_val(4,:) );
        end
        
    end% methods
    
end

