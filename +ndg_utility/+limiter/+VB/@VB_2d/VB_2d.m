classdef VB_2d < ndg_utility.limiter.limiter_vert
    %VB_2D 二维节点限制器
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        xc, yc % 单元形心坐标
        vx, vy % 每个单元顶点坐标
    end
    
    methods(Access=private)
        function [ flim ] = HWENO_wei(obj, f_Q, fmean, fmax, fmin)
            flim = vb_weno(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
        
        function [ flim ] = VA_wei(obj, f_Q, fmean, fmax, fmin)
            flim = vb_va(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
        
        function [ flim ] = JK_wei(obj, f_Q, fmean, fmax, fmin)
            flim = vb_jk(f_Q, obj.mesh.x, obj.mesh.y, ...
                fmean, obj.xc, obj.yc, fmax, fmin, ...
                obj.mesh.EToV, obj.cell.Fmask);
        end
    end
    
    methods
        function obj = VB_2d(mesh)
            obj = obj@ndg_utility.limiter.limiter_vert(mesh, mesh.cell);
            obj.xc = mesh.cell_mean(mesh.x);
            obj.yc = mesh.cell_mean(mesh.y);
            obj.vx= mesh.x(mesh.EToV);
            obj.vy= mesh.y(mesh.EToV);
        end
        
        function [ flimt ] = limit(obj, f_Q)
            fmean = obj.mesh.cell_mean(f_Q);
            [fmax, fmin] = vertex_extreme(obj.Kv, obj.VToE, fmean);
            [ flimt ] = obj.HWENO_wei(f_Q, fmean, fmax, fmin);
        end
    end
    
end

