classdef conv2d_diffusion < conv2d
    %CONV2D_DIFFUSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        x0 = -0.5; 
        y0 = -0.5;
        miu = 0;
    end
    
    methods
        function obj = conv2d_diffusion(N, M, type)
            mesh = uniform_mesh(N, M, type);
            obj = obj@conv2d(mesh);
            obj.ftime = 2;
            obj.init();
        end
        
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end
        
        function init(obj)
            obj.u = 0.5*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.v = 0.5*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.f_Q = obj.ext_func(0);
        end
        
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u*time;
            yc = obj.y0 + obj.v*time;
            if obj.miu > 0
                t = -(obj.mesh.x-xc).^2/obj.miu ...
                    -(obj.mesh.y-yc).^2/obj.miu;
            else
                sigma = 125*1e3/(33*33);
                t = -( (obj.mesh.x-xc).^2+(obj.mesh.y-yc).^2 )*sigma;
            end
            f_ext = exp(t);
        end
    end
    
end

