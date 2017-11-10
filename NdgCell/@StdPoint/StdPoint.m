classdef StdPoint < StdCell
    properties(Constant)
        type = NdgCellType.Point
        Nv = 1      
        LAV = 0
        vr = 0      
        vs = 0      
        vt = 0      
        Nfv = 1     
        FToV = 1    
        Nface = 0
        faceType = NdgCellType.Point
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
        
        function [ Nq,rq,sq,tq,wq ] = quad_coor_func(obj, N)
            Nq = 1;
            rq = 0;
            sq = 0;
            tq = 0;
            wq = 0;
        end
    end
    
    methods
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

