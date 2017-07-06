classdef conv2d < ndg_lib.phys.phys2d
    %CONV2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 1  % 变量个数
    end
    
    properties
        ftime   % 计算终止时间
        dt      % 计算时间步长
        u,v     % 速度场
    end
    
    %% 虚函数
    methods(Abstract)
        [ spe ] = character_len(obj, f_Q) % get the time interval dt
    end
    
    %% 私有函数
    methods(Access=protected) % private 
        [ E, G ] = flux_term( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
    end
    
    %% 公共方法
    methods
        function obj = conv2d(mesh)
            obj = obj@ndg_lib.phys.phys2d(mesh);
        end
        
        f_Q = RK45_solve(obj) % Runge-Kutta 4th order 5 stages
        refine_mesh(obj, multi_ratio) % refined mesh
    end
    
end

