classdef StdPoint < StdCell
    properties(Constant)
        type = enumStdCell.Point
        Nv = 1      
        LAV = 0
        vr = 0      
        vs = 0      
        vt = 0      
        Nfv = 1     
        FToV = 1    
        Nface = 0
        faceType = enumStdCell.Point
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
        
        function [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z )
            nx = 0;
            ny = 0;
            nz = 0;
            Js = 1;
        end
        
        function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobianMatrix( obj, x, y, z )
            rx = ones( size(x) );
            ry = ones( size(x) );
            rz = ones( size(x) );
            sx = ones( size(x) );
            sy = ones( size(x) );
            sz = ones( size(x) );
            tx = ones( size(x) );
            ty = ones( size(x) );
            tz = ones( size(x) );
            J = ones( size(x) );
        end
        
        function fun = orthogonal_func(obj, N, ind, r, s, t)
            fun = 1;
        end
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = vert_val;
        end
        
    end
    
end

