classdef phys1d < ndg_lib.phys.phys
    %PHYS1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        Nfield  % 变量个数
    end
    
    %% 误差统计
    methods
        
    end
    
    methods
        function obj = phys1d(mesh)
            obj = obj@ndg_lib.phys.phys(mesh);
        end
        
        function draw(obj, fld)
            obj.draw_h = plot(obj.mesh.x, obj.f_Q(:, :, fld));
        end
    end
    
end

