%> @brief Standard Line class.
%
%> Define the properties for the StdLine class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University
%> @email li12242@tju.edu.cn
% ======================================================================
classdef StdLine < StdCell
    
    properties(Constant)
        type = StdCellType.Line
        Nv = 2
        vol = 2
        vr = [ -1, 1]'
        vs = [ 0, 0 ]';
        vt = [ 0, 0 ]';
        Nfv = [ 1, 1 ]';
        FToV = [ 1, 2 ];
        Nface = 2;
        faceType = [ StdCellType.Point, StdCellType.Point ];
    end
    
    methods(Access=protected)
        function [Np,r,s,t] = node_coor_func(obj, N)
            Np = N+1;
            [ r,~ ] = Polylib.zwglj(Np);
            s = zeros(Np, 1);
            t = zeros(Np, 1);
        end
        
        function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
            dr = Polylib.GradJacobiP(r, 0, 0, ind-1);
            ds = zeros(obj.Np, 1);
            dt = zeros(obj.Np, 1);
        end
        
        function [Nq, rq, sq, tq, wq] = quadrature_node_func(obj, N)
            Nq = N+1;
            [ rq, wq ] = Polylib.zwglj( Nq );
            sq = zeros(Nq, 1);
            tq = zeros(Nq, 1);
        end% func
    end
    
    methods
        function obj = StdLine(N)
            obj = obj@StdCell(N);
        end
        
        function f = orthogonal_func(obj, N, ind, r, s, t)
            f = Polylib.JacobiP(r, 0, 0, ind-1);
        end% func
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*((1-obj.r)*vert_val(1,:) + (1+obj.r)*vert_val(2,:));
        end
        
    end% methods
    
end

