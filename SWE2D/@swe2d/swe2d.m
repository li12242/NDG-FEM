classdef swe2d < ndg_lib.phys.phys2d
    %@SWE2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 3
        gra = 9.81
        hmin = 1e-4
    end
    
    properties(Abstract, SetAccess=protected)
        bot     % 底坡高程
        cfl     % CFL 数
        ftime   % 计算终止时间
        dt      % 计算时间步长
        wetflag % 湿单元逻辑数组
        slopelimiter % 斜率限制器
    end
    properties
        m   % Manning 参数
    end
    
    methods(Abstract)
        [h, qx, qy] = init(obj, x, y)
        [hP, qxP, qyP] = adj_node_val(obj, ... 
            hM, qxM, qyM, ...
            hP, qxP, qyP, ...
            nx, ny, ftype) % adjacent node value of each edges
    end
    
    methods
        rhs = swe_rhs(obj, h, qx, qy, time) % calculate r.h.s.
        [E,G] = swe_node_flux(obj, h, qx, qy) % get flux terms
        flux = swe_num_flux(obj, h, qx, qy) % get numerical flux on edges
    end
    
end

