classdef std_line < ndg_lib.std_cell.std_cell
    %STD_LINE Summary of this class goes here
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
        function obj = std_line(N)
            obj.N = N;
            obj.type = ndg_lib.std_cell_type.Line;
            % basic
            obj.Np = N+1;
            obj.Nv = 2;
            obj.vr = [-1, 1]';
            obj.vs = [0, 0]';
            obj.vt = [0, 0]';
            obj.FToV = ndg_utility.ndg_list([1,1]);
            obj.FToV(1) = 1; obj.FToV(2) = 2;
            obj.Nface = 2;
            point = ndg_lib.std_cell_type.Point;
            obj.faceType = [point, point];
            % volume
            [obj.r, obj.s, obj.t] = obj.node_coor_func(N);
            obj.V = obj.Vandmode_Matrix(@obj.orthogonal_func);
            obj.M = obj.Mass_Matrix;
            [obj.Dr, obj.Ds, obj.Dt] = obj.Derivative_Matrix...
                (@obj.derivative_orthogonal_func);
            % face
            cell = ndg_lib.ndg_cell(obj.N, obj.faceType(1));
            obj.Nfp = [cell.Np, cell.Np];
            obj.Nfptotal = sum(obj.Nfp);
            obj.Fmask = obj.Face_node;
            obj.LIFT = obj.Lift_Matrix;
        end
        
        function [r,s,t] = node_coor_func(obj, N)
            [r,~] = Polylib.zwglj(N+1);
            s = zeros(obj.Np, 1);
            t = zeros(obj.Np, 1);
        end
        
        function f = orthogonal_func(obj, N, ind, r, s, t)
            f = Polylib.JacobiP(r, 0, 0, ind-1);
        end% func
        
        function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
            dr = Polylib.GradJacobiP(r, 0, 0, ind-1);
            ds = zeros(obj.Np, 1);
            dt = zeros(obj.Np, 1);
        end
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*((1-obj.r)*vert_val(1,:) + (1+obj.r)*vert_val(2,:));
        end
        
    end% methods
    
end

