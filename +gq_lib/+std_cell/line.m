classdef line < ndg_lib.std_cell.line & gq_lib.std_cell.gauss_quad_cell
    %LINE Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = line(N)
            obj = obj@ndg_lib.std_cell.line(N);
            obj = obj@gq_lib.std_cell.gauss_quad_cell(N);
            
        end
    end
    methods(Access=protected)
        function [rq, sq, tq, wq] = gaussquad_vol_coor(obj, N)
            [rq, wq] = Polylib.zwgl(N+1);
            sq = zeros(size(rq)); 
            tq = zeros(size(rq)); 
        end
        
        function [rbq, sbq, tbq, wbq] = gaussquad_surf_coor(obj, N)
            rbq = zeros(obj.Nfq, 1);
            wbq = zeros(obj.Nfq, 1);
            sind = 1;
            for f = 1:obj.Nface
                cell = gq_lib.std_cell.get_std_cell(N, obj.faceType(f) );
                eind = sind + cell.Nq - 1;
                r = obj.r(obj.Fmask(:, f));
                rbq(sind:eind) = cell.project_node2quad(r);
                wbq(sind:eind) = cell.wq;
                sind = eind + 1;
            end
            rbq = rbq(:);
            wbq = wbq(:);
            sbq = zeros(size(rbq)); 
            tbq = zeros(size(rbq)); 
        end
    end
    
end

