classdef StdQuad < StdCell

    properties(Constant)
        type = NdgCellType.Quad
        Nv = 4
        LAV = 4
        vr = [-1,  1,  1, -1]'
        vs = [-1, -1,  1,  1]'
        vt = [ 0,  0,  0,  0]'
        Nfv = [2,2,2,2]'
        FToV = [1,2; 2,3; 3,4; 4,1]'
        Nface = 4
        faceType = [...
            NdgCellType.Line, ...
            NdgCellType.Line, ...
            NdgCellType.Line, ...
            NdgCellType.Line]
    end
    
    methods(Access=protected)
        [Np, r,s,t] = node_coor_func(obj, N);
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
        [Nq, rq, sq, tq, wq] = quad_coor_func(obj, N);
    end
    
    methods
        function obj = StdQuad(N)
            obj = obj@StdCell(N);
        end

        f = orthogonal_func(obj, N, ind, r, s, t);
        
        function node_val = project_vert2node(obj, vert_val)
            node_val = 0.25*(...
                (1-obj.r).*(1-obj.s)*vert_val(1,:) + ...
                (1+obj.r).*(1-obj.s)*vert_val(2,:) + ...
                (1+obj.r).*(1+obj.s)*vert_val(3,:) + ...
                (1-obj.r).*(1+obj.s)*vert_val(4,:) );
        end
        
    end% methods
    
end

