classdef conv2d
    %CONV2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        c
        u
        v
    end
    
    methods(Abstract)
        c_ext = ext_func(obj, time)
        c = init(obj, x, y)
    end
    
    methods
        rhs = rhs_term(obj, c, time)
        [E, G] = node_flux_term(obj, c)
        flux = num_flux(obj, c)
        cP = adj_val(obj, cM, cP, nx, ny, ftype)
    end
end

