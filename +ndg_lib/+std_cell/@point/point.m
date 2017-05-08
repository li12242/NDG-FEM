classdef point < ndg_lib.std_cell.std_cell
    %STD_POINT Summary of this class goes here
    %   Detailed explanation goes here
    %% 全局属性
    properties
        N   % 基函数阶数
    end
    
    % 基本属性
    properties(Constant)
        type = ndg_lib.std_cell_type.Point % 单元类型
        Nv = 1      % 单元顶点个数
        vol = 0
        vr = 0      % 顶点x坐标，列优先
        vs = 0      % 顶点y坐标，列优先
        vt = 0      % 顶点z坐标，列优先
        Nfv = 1     % 每个面上顶点个数
        FToV = 1    % 每个面对应的顶点编号，列优先
        Nface = 0
        faceType = ndg_lib.std_cell_type.Point
    end
    
    % 体积属性
    properties(SetAccess = private)
        Np      % 节点个数
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
        function [r,s,t] = node_coor_func(obj, N)
            r = 0;
            s = 0;
            t = 0;
        end
        
        function fun = orthogonal_func(obj, N, ind, r, s, t)
            fun = 1;
        end
        
        function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
            dr = 0;
            ds = 0;
            dt = 0;
        end
    end
    
    methods
        function obj = point(N)
            obj.N = N;
            % volume
            obj.Np = 1;
            [obj.r, obj.s, obj.t] = obj.node_coor_func(N);
            obj.V = obj.Vandmode_Matrix(@obj.orthogonal_func);
            obj.M = obj.Mass_Matrix;
            [obj.Dr, obj.Ds, obj.Dt] = obj.Derivative_Matrix...
                (@obj.derivative_orthogonal_func);
            % face
            obj.Nfp = 1;
            obj.Nfptotal = sum(obj.Nfp);
            obj.Fmask = obj.Face_node;
            obj.LIFT = obj.Lift_Matrix;
        end
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = vert_val;
        end
        
    end
    
end

