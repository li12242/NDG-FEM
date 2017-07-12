classdef quad < ndg_lib.std_cell.quad & gq_lib.std_cell.gauss_quad_cell
    %QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = quad(N)
            obj = obj@ndg_lib.std_cell.quad(N);
            obj = obj@gq_lib.std_cell.gauss_quad_cell(N);
        end
    end
    
    methods(Access=protected)
        function [rq, sq, tq, wq] = gaussquad_vol_coor(obj, N)
            
            np = N+1;
            [zq, w] = Polylib.zwgl(np); % 一维 GL 积分节点
            % 首先按照x坐标循环，然后按照y坐标循环
            rq = zq*ones(1, np);
            sq = ones(np, 1)*zq';
            rq = rq(:); sq = sq(:); 
            wq = bsxfun(@times, w, w'); wq = wq(:);
            tq = zeros(size(rq));
        end
        
        function [rbq, sbq, tbq, wbq] = gaussquad_surf_coor(obj, N)
            sind = 1;
            for f = 1:obj.Nface
                cell = gq_lib.std_cell.get_std_cell(N, obj.faceType(f) );
                obj.Nfq = obj.Nfq + cell.Nq;
                eind = sind + cell.Nq - 1;
                r = obj.r(obj.Fmask(:, f));
                s = obj.s(obj.Fmask(:, f));
                rbq(sind:eind) = cell.project_node2quad(r);
                sbq(sind:eind) = cell.project_node2quad(s);
                wbq(sind:eind) = cell.wq;
                sind = eind + 1;
            end
            rbq = rbq(:);
            sbq = sbq(:);
            wbq = wbq(:);
            tbq = zeros(size(rbq));
        end
    end
    
end

