%> @brief Standard triangle class.
%
%> The StdTri define the properties for the standard triangle element.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University
%> @email li12242@tju.edu.cn
% ======================================================================
classdef StdTri < StdCell

    properties(Constant)
        type = StdCellType.Tri
        Nv = 3
        vol = 2
        vr = [-1,  1, -1]'
        vs = [-1, -1,  1]'
        vt = [ 0,  0,  0]'
        Nfv = [2,2,2];
        FToV = [1,2; 2,3; 3,1]';
        Nface = 3
        faceType = [...
            StdCellType.Line,...
            StdCellType.Line,...
            StdCellType.Line]
    end
    
    methods(Access=protected)
        [Np,r,s,t] = node_coor_func(obj, N);
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
        [Nq, rq, sq, tq, wq] = quadrature_node_func(obj, N);
    end
    
    methods
        function obj = StdTri(N)
            obj = obj@StdCell(N);
        end
        
        f = orthogonal_func(obj, N, ind, r, s, t);
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.5*(-(obj.r+obj.s)*vert_val(1, :) ...
                + (1+obj.r)*vert_val(2, :)...
                + (1+obj.s)*vert_val(3, :));
        end
        
    end% methods
    
end

