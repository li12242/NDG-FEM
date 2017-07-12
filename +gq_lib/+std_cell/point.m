classdef point < ndg_lib.std_cell.point & gq_lib.std_cell.gauss_quad_cell
    %POINT Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = point(N)
            obj = obj@gq_lib.std_cell.gauss_quad_cell(N);
            obj = obj@ndg_lib.std_cell.point(N);
            obj.Nq = 1; obj.Nfq = 1;
        end
    end
    methods(Access=protected)
        function [rq, sq, tq, wq] = gaussquad_vol_coor(obj, N)
            rq = 0; sq = 0; tq = 0; wq = 1;
        end
        
        function [rbq, sbq, tbq, wbq] = gaussquad_surf_coor(obj, N)
            rbq = 0; sbq = 0; tbq = 0; wbq = 1;
        end
    end
    
end

