classdef swe1d < ndg_lib.phys.phys1d
    %SWE1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 2
        hmin = 1e-3
    end
    
    properties
        bot     % 底坡高程
        cfl     % CFL 数
        ftime   % 计算终止时间
        dt      % 计算时间步长
    end
    
    %% 虚函数
    methods(Abstract)
        init(obj)
    end
    
    %% 公共函数
    methods
        [ dflux ] = hll_surf_term( obj, f_Q )
        [ dflux ] = roe_surf_term( obj, f_Q )
        [ rhs ] = rhs_term( obj, f_Q )
        [ E ] = flux_term( obj, f_Q )
        [ S ] = source_term( obj, f_Q )
    end
    
end

