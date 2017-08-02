classdef conv2d_diffusion < conv2d_advection
    %CONV2D_DIFFUSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        miu = 0.05
    end
    
    %% private methods
    methods(Access=private)
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u*time;
            yc = obj.y0 + obj.v*time;
            t = -(obj.mesh.x-xc).^2/obj.miu -(obj.mesh.y-yc).^2/obj.miu;
            f_ext = exp(t);
        end
    end% func
    
    %% public methods
    methods        
        function obj = conv2d_diffusion(varargin)
            obj = obj@conv2d_advection(varargin);
        end
        
    end
    
end

