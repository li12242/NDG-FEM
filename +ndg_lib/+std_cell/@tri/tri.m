classdef tri < ndg_lib.std_cell.std_cell
    %STD_TRI Summary of this class goes here
    %   Detailed explanation goes here
    
    %% 全局属性
    properties
        N   % 基函数阶数
    end
    
    % 基本属性
    properties(Constant)
        type = ndg_lib.std_cell_type.Tri % 单元类型
        Nv = 3;
        vr = [-1,  1, -1]'
        vs = [-1, -1,  1]'
        vt = [ 0,  0,  0]'
        Nfv = [2,2,2];
        FToV = [1,2; 2,3; 3,1]';
        Nface = 3
        faceType = [...
            ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Line]
    end
    
    % 体积属性
    properties(SetAccess = private)
        Np % 节点个数
        r, s, t
        V
        M
        Dr, Ds, Dt
    end
    % 面属性
    properties(SetAccess = private)
        Fmask
        Nfp
        Nfptotal
        LIFT
    end
    
    methods(Access=protected)
        [r,s,t] = node_coor_func(obj, N);
        f = orthogonal_func(obj, N, ind, r, s, t);
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
    end
    
    methods
        function obj = tri(N)
            obj.N = N;
            % volume
            obj.Np = (N+1)*(N+2)/2;
            [obj.r, obj.s, obj.t] = obj.node_coor_func(N);
            obj.V = obj.Vandmode_Matrix(@obj.orthogonal_func);
            obj.M = obj.Mass_Matrix;
            [obj.Dr, obj.Ds, obj.Dt] = obj.Derivative_Matrix...
                (@obj.derivative_orthogonal_func);
            % face
            cell = ndg_lib.ndg_cell(obj.N, obj.faceType(1));
            obj.Nfp = [cell.Np, cell.Np, cell.Np];
            obj.Nfptotal = sum(obj.Nfp);
            obj.Fmask = obj.Face_node();
            obj.LIFT = obj.Lift_Matrix();
        end
       
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*(-(obj.r+obj.s)*vert_val(1, :) ...
                + (1+obj.r)*vert_val(2, :)...
                + (1+obj.s)*vert_val(3, :));
        end
        
    end% methods
    
end

