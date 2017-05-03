classdef std_cell < handle
    %STD_CELL Summary of this class goes here
    %   Detailed explanation goes here
    
    %% 全局属性
    properties(Abstract)
        N   % 基函数阶数
        type % 单元类型
    end
    
    % 基本属性
    properties(Abstract, SetAccess = private)
        Np % 节点个数
        Nv
        vr, vs, vt
        FToV
        Nface
        faceType
    end
    
    % 体积属性
    properties(Abstract, SetAccess = private)
        r, s, t
        V
        M
        Dr, Ds, Dt
    end
    % 面属性
    properties(Abstract, SetAccess = private)
        Fmask
        Nfp
        Nfptotal
        LIFT
    end
    
    %% 方法
    methods(Abstract) % 虚函数
        [r,s,t] = node_coor_func(obj, N)
        fun = orthogonal_func(obj, N, ind, r, s, t)
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
        node_val = project_vert2node(obj, vert_val);
    end
   
    % 公共方法
    methods(Hidden, Access = protected)
        function V = Vandmode_Matrix(obj, orthogonal_func)
            V = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                V(:, n) = orthogonal_func(obj.N, n, obj.r, obj.s, obj.t);
            end% for
        end% func
        
        function M = Mass_Matrix(obj)
            invV = inv(obj.V);
            M = (invV')*invV;
        end
        
        function [Dr, Ds, Dt] = Derivative_Matrix(obj, deri_orthogonal_func)
            Vr = zeros(obj.Np, obj.Np);
            Vs = zeros(obj.Np, obj.Np);
            Vt = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = deri_orthogonal_func...
                    (obj.N, n, obj.r, obj.s, obj.t);
            end
            
            Dr = Vr/obj.V; Ds = Vs/obj.V; Dt = Vt/obj.V;
        end% func
        
        function Fmask = Face_node(obj)
            Fmask = ndg_utility.ndg_list(obj.Nfp);
            for f = 1:obj.Nface
                % get vertex index on face f
                rv = obj.vr(obj.FToV(:, f));
                sv = obj.vs(obj.FToV(:, f));
                tv = obj.vt(obj.FToV(:, f));
                if(isrow(rv)) rv = rv'; end
                if(isrow(sv)) sv = sv'; end
                if(isrow(tv)) tv = tv'; end
                % get the nodes on face f
                cell = ndg_lib.ndg_cell(obj.N, obj.faceType(f));
                fr = cell.project_vert2node(rv);
                fs = cell.project_vert2node(sv);
                ft = cell.project_vert2node(tv);
                % get the nodes index
                for n = 1:obj.Nfp(f)
                    dis = (fr(n) - obj.r).^2 + (fs(n) - obj.s).^2 + (ft(n) - obj.t).^2;
                    ind = find(dis < 1e-10);
                    Fmask(n, f) = ind;
                end
            end
        end% func
        
        function LIFT = Lift_Matrix(obj)
            Mes = zeros(obj.Np, obj.Nfptotal);
            sk = 1;
            for f = 1:obj.Nface
                cell = ndg_lib.ndg_cell(obj.N, obj.faceType(f));
                row = obj.Fmask(:, f);
                for n = 1:cell.Np
                    Mes(row, sk) = cell.M(:, n);
                    sk = sk + 1;
                end
            end
            LIFT = (obj.V*(obj.V)')*Mes;
        end
    end
end

