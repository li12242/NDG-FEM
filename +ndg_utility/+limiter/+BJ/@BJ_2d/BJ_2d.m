classdef BJ_2d < ndg_utility.limiter.BJ.BJ
    %BJ_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function [ obj ] = BJ_2d(mesh)
            obj = obj@ndg_utility.limiter.BJ.BJ(mesh);
        end
        
        function [ f_Q ] = limit(obj, f_Q)
            [ c_mean ] = obj.mesh.cell_mean( f_Q );
            [ vmax, vmin ] = vertex_bound(obj, c_mean);
            [ f_Q ] = BJ_limit_2d(f_Q, obj.mesh.x, obj.mesh.y, ...
                c_mean, obj.xc, obj.yc, obj.mesh.vol, ...
                vmin, vmax, obj.cell.Fmask, ...
                obj.mesh.EToV, obj.mesh.Js, obj.mesh.cell.ws);
        end% func
    end
    
end
