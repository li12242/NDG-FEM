classdef VB_TVB < ndg_utility.limiter.limiter_vert
    %VB_TVB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VToC % 每个顶点与所有相邻单元插值系数
        xc, yc
    end
    
    methods(Hidden)
        function VToC = inverse_distance_weight(obj)
            % 计算顶点插值系数
            VToC = zeros( size(obj.VToE) );
            for n = 1:obj.mesh.Nv
                eind = obj.VToE( 1:obj.Kv(n), n );
                xm = obj.xc(eind); 
                ym = obj.yc(eind);
                w = 1./( (obj.mesh.vx(n) - xm).^2 ...
                    + (obj.mesh.vy(n) - ym).^2 );
                VToC( 1:obj.Kv(n), n ) = w./sum(w);
            end
        end
    end
    
    methods
        function obj = VB_TVB(mesh)
            obj = obj@ndg_utility.limiter.limiter_vert(mesh, mesh.cell);
            obj.xc = mesh.cell_mean(mesh.x);
            obj.yc = mesh.cell_mean(mesh.y);
            obj.VToC = inverse_distance_weight(obj);
        end
        
        function f_Q = limit(obj, f_Q, factor)
            f_m = obj.mesh.cell_mean(f_Q);
            f_Q = VB_TVB(f_Q, obj.mesh.x, obj.mesh.y, ...
                f_m, obj.xc, obj.yc, obj.Kv, obj.VToE, obj.VToC,...
                obj.mesh.EToV, obj.mesh.cell.Fmask, factor);
        end
    end
    
end

