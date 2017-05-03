classdef std_point < ndg_lib.std_cell.std_cell
    %STD_POINT Summary of this class goes here
    %   Detailed explanation goes here
    %% 全局属性
    properties
        N   % 基函数阶数
        type % 单元类型
    end
    
    % 基本属性
    properties(SetAccess = private)
        Np % 节点个数
        Nv
        vr, vs, vt
        FToV
        Nface
        faceType
    end
    
    % 体积属性
    properties(SetAccess = private)
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
    
    methods
        function obj = std_point(N)
            obj.N = N;
            obj.type = ndg_lib.std_cell_type.Point;
            % basic
            obj.Np = 1;
            obj.Nv = 1;
            obj.vr = 0;
            obj.vs = 0;
            obj.vt = 0;
            obj.FToV = 1;
            obj.Nface = 0;
            obj.faceType = ndg_lib.std_cell_type.Point;
            % volume
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
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = vert_val;
        end
        
    end
    
end

