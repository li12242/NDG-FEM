classdef conv2d_diffusion < conv2d
    %CONV2D_DIFFUSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        u = 0.5; 
        v = 0.5;
        x0 = -0.5; 
        y0 = -0.5;
    end
    
    methods
        function obj = conv2d_diffusion(N, M, type)
            mesh = uniform_mesh(N, M, type);
            
        end
    end
    
end

