%> @brief Standard cell class of StdPoint.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University
%> @email li12242@tju.edu.cn
% ======================================================================
classdef StdPoint < StdCell
    
    properties(Constant)
        type = StdCellType.Point
        Nv = 1
        vol = 0
        vr = 0
        vs = 0
        vt = 0
        Nfv = 1
        FToV = 1
        Nface = 0
        faceType = StdCellType.Point
    end
    
    methods(Access=protected)
        function [Np,r,s,t] = node_coor_func(obj, N)
            Np = 1;
            r = 0;
            s = 0;
            t = 0;
        end
        
        function [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
            dr = 0;
            ds = 0;
            dt = 0;
        end
        
        function [Nq, rq, sq, tq, wq] = quadrature_node_func(obj, N)
            Nq = 1;
            rq = zeros(Nq, 1);
            sq = zeros(Nq, 1);
            tq = zeros(Nq, 1);
            wq = ones(Nq, 1);
        end% func
    end
    
    methods
        %> Construct the object with specific order N
        function obj = StdPoint(N)
            obj = obj@StdCell(N);
        end
        
        function fun = orthogonal_func(obj, N, ind, r, s, t)
            fun = 1;
        end
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = vert_val;
        end
        
    end
    
end

