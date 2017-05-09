classdef conv2d
    %CONV2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=protected)
        mesh
        c_ext
        u
        v
    end
    
    properties
        cfl % CFL number
        c
        miu
        ftime
        dt
    end
    
    % abstract function
    methods(Abstract)
        [ c, u ] = init(obj, x) % get the initial value
        [ dt ] = time_interval(obj) % get the time interval dt
    end
    
    % function - matlab
    methods(Access=protected) % private 
        [ cP ] = nei_node_val(obj, cM, cP, nx, ny, ftype) % adjacent node value
        [ E, G ] = flux_term( obj, u, c ) % get the flux terms
        [ flux ] = num_flux(obj, cM, uM, cP, uP, nx) % get numerical flux
        [ rhs ] = rhs_term(obj, c, time) % get the r.h.s term
    end
    methods % public solver
        [ c ] = solve_mat(obj)
    end
    % function - mat-mex
    methods(Access=protected) % private
        
    end
    methods % public
        [ c ] = solve_mex(obj)
    end
    
    % private function - mat-gpu
    methods(Access=protected) % private
    end
    methods % public
        [ c ] = solve_gpu(obj)
    end
    
    % private function - cuda-gpu
    methods(Access=protected) % private
    end
    methods % public
        [ c ] = solve_cuda(obj)
    end
    
    methods
    end
end

