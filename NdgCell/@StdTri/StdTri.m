%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef StdTri < StdCell
    properties(Constant)
        type = enumStdCell.Tri;
        Nv = 3
        LAV = 2
        vr = [-1,  1, -1]'
        vs = [-1, -1,  1]'
        vt = [ 0,  0,  0]'
        Nfv = [2,2,2];
        FToV = [1,2; 2,3; 3,1]';
        Nface = 3
        faceType = [enumStdCell.Line, ...
                    enumStdCell.Line, ...
                    enumStdCell.Line]
    end
    
    methods(Access=protected)
        [Np,r,s,t] = node_coor_func(obj, N);
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
        [Nq, rq, sq, tq, wq] = quad_coor_func( obj, qOrd );
    end
    
    methods
        function obj = StdTri(N)
            obj = obj@StdCell(N);
        end
        
        [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z );
        
        %> @brief Return the jacobian matrix and its determination.
        %> 
        function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobianMatrix( obj, x, y, z )
            xr = obj.Dr * x; xs = obj.Ds * x;
            yr = obj.Dr * y; ys = obj.Ds * y;
            J = -xs.*yr + xr.*ys;
            
            rx = ys./J; sx =-yr./J;
            ry =-xs./J; sy = xr./J; 
            
            rz = ones( size(x) );
            sz = ones( size(x) );
            tx = ones( size(x) );
            ty = ones( size(x) );
            tz = ones( size(x) );
        end
        
        f = orthogonal_func(obj, N, ind, r, s, t);

        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*(-(obj.r+obj.s)*vert_val(1, :) ...
                + (1+obj.r)*vert_val(2, :)...
                + (1+obj.s)*vert_val(3, :));
        end
        
    end% methods
    
end

