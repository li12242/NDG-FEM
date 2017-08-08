classdef line < ndg_lib.std_cell.std_cell
    %LINE Summary of this class goes here
    %   Detailed explanation goes here
    
    %% 基本属性
    properties(Constant)
        type = ndg_lib.std_cell_type.Line % 单元类型
        Nv = 2
        vol = 2
        vr = [-1, 1]'
        vs = [ 0, 0]'; 
        vt = [ 0, 0]'; 
        Nfv = [1, 1]';
        FToV = [1,2];
        Nface = 2;
        faceType = [...
            ndg_lib.std_cell_type.Point, ...
            ndg_lib.std_cell_type.Point];
    end
    
    methods(Access=protected)
        function [Np,r,s,t] = node_coor_func(obj, N)
            Np = N+1;
            [r,~] = Polylib.zwglj(Np);
            s = zeros(Np, 1);
            t = zeros(Np, 1);
        end
        
        function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
            dr = Polylib.GradJacobiP(r, 0, 0, ind-1);
            ds = zeros(obj.Np, 1);
            dt = zeros(obj.Np, 1);
        end
    end
    
    methods
        function obj = line(N)
            obj = obj@ndg_lib.std_cell.std_cell(N);
        end
        
        function f = orthogonal_func(obj, N, ind, r, s, t)
            f = Polylib.JacobiP(r, 0, 0, ind-1);
        end% func
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*((1-obj.r)*vert_val(1,:) + (1+obj.r)*vert_val(2,:));
        end
        
    end% methods
    
end

