classdef std_cell
    %STD_CELL Summary of this class goes here
    %   Detailed explanation goes here
    
    %% 全局属性
    properties
        N   % 基函数阶数
    end
    
    % 基本属性
    properties(Constant, Abstract)
        type % ndg_lib.std_cell_type 单元类型
        Nv
        vr, vs, vt
        vol % 标准单元长度/面积/体积
        Nfv
        FToV
        Nface
        faceType
    end
    
    % 体积属性
    properties(SetAccess = protected)
        Np % 节点个数
        r, s, t
        V
        M
        Dr, Ds, Dt
        w
    end
    % 面属性
    properties(SetAccess = protected)
        Fmask
        Nfp
        Nfptotal
        LIFT
        ws
    end
    
    %% 虚函数
    methods(Abstract, Access=protected) % 私有
        [Np, r,s,t] = node_coor_func(obj, N)
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
    end
    
    %% 公共虚函数
    methods(Abstract, Access=public)
        fun = orthogonal_func(obj, N, ind, r, s, t);
        node_val = project_vert2node(obj, vert_val);
    end
    
    %% 构造函数
    methods
        function obj = std_cell(N)
            obj.N = N;
            % volume
            [obj.Np, obj.r, obj.s, obj.t] = obj.node_coor_func(N);
            obj.V = obj.Vandmode_Matrix(@obj.orthogonal_func);
            obj.M = obj.Mass_Matrix();
            [obj.Dr, obj.Ds, obj.Dt] = obj.Derivative_Matrix...
                (@obj.derivative_orthogonal_func);
            % face
            if obj.Nface>0
                obj.Nfp = zeros(obj.Nface, 1);
            else
                obj.Nfp = 1;
            end
            for f = 1:obj.Nface
                cell = ndg_lib.get_std_cell(obj.N, obj.faceType(f));
                obj.Nfp(f) = cell.Np;
            end
            obj.Nfptotal = sum(obj.Nfp);
            obj.Fmask = obj.Face_node;
            obj.LIFT = obj.Lift_Matrix;
            [obj.w, obj.ws] = obj.integral_coeff();
        end
    end
   
    %% 隐藏私有函数
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
            maxnfp = max(obj.Nfp);
            Fmask = zeros(maxnfp, obj.Nface);
            for f = 1:obj.Nface
                nfv = obj.Nfv(f);
                % get vertex index on face f
                rv = obj.vr(obj.FToV(1:nfv, f));
                sv = obj.vs(obj.FToV(1:nfv, f));
                tv = obj.vt(obj.FToV(1:nfv, f));
                if(isrow(rv)) rv = rv'; end
                if(isrow(sv)) sv = sv'; end
                if(isrow(tv)) tv = tv'; end
                % get the nodes on face f
                cell = ndg_lib.get_std_cell(obj.N, obj.faceType(f));
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
                cell = ndg_lib.get_std_cell(obj.N, obj.faceType(f));
                row = obj.Fmask(:, f);
                row = row(row ~= 0);
                for n = 1:cell.Np
                    try
                    Mes(row, sk) = cell.M(:, n);
                    catch
                        keyboard
                    end
                    sk = sk + 1;
                end
            end
            LIFT = (obj.V*(obj.V)')*Mes;
        end
        
        function [w, ws] = integral_coeff(obj)
            w = zeros(obj.Np, 1);
            w(:) = sum(obj.M);
            ws = zeros(obj.Nfptotal, 1);
            Mes = zeros(obj.Np, obj.Nfptotal);
            sk = 1;
            for f = 1:obj.Nface
                cell = ndg_lib.get_std_cell(obj.N, obj.faceType(f));
                row = obj.Fmask(:, f);
                row = row(row ~= 0);
                for n = 1:cell.Np
                    Mes(row, sk) = cell.M(:, n);
                    sk = sk + 1;
                end
            end
            ws(:) = sum(Mes);
        end
    end
end

