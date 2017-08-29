classdef VB_2d < ndg_utility.limiter.limiter_vert
    %VB_2D 2d vertex-based slope limiter.
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        VToC
        xc, yc % 单元形心坐标
        vx, vy % 每个单元顶点坐标
    end
    
    methods(Access=private, Hidden)
        function [ flim ] = HWENO_wei(obj, f_Q, fmean, fv, fmax, fmin)
            flim = vb_weno(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fv, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
        
        function [ flim ] = VA_wei(obj, f_Q, fmean, fv, fmax, fmin)
            flim = vb_va(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fv, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
        
        function [ flim ] = JK_wei(obj, f_Q, fmean, fv, fmax, fmin)
            flim = vb_jk(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fv, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
        
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
        function obj = VB_2d(mesh)
            obj = obj@ndg_utility.limiter.limiter_vert(mesh, mesh.cell);
            obj.xc = mesh.cell_mean(mesh.x);
            obj.yc = mesh.cell_mean(mesh.y);
            obj.vx= mesh.x(mesh.EToV);
            obj.vy= mesh.y(mesh.EToV);
            obj.VToC = inverse_distance_weight(obj);
        end
        
        function [ flimt ] = limit(obj, f_Q)
            c_mean = obj.mesh.cell_mean(f_Q);
            [fmax, fmin] = vertex_extreme(obj.Kv, obj.VToE, c_mean);
            [ fv ] = vertex_average(c_mean, obj.Kv, obj.VToE, obj.VToC);
            [ flimt ] = obj.HWENO_wei(f_Q, c_mean, fv, fmax, fmin);
        end
    end
    
end

