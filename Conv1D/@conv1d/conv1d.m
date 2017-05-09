classdef conv1d < matlab.mixin.SetGet
    %CONV1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=protected)
        mesh
        c_ext
        u
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
        [c, u] = init(obj, x) % get the initial value
        dt = time_interval(obj) % get the time interval dt
    end
    
    % function - matlab
    methods(Access=protected) % private 
        cP = nei_node_val(obj, cM, cP, c_extM, ftype) % adjacent node value
        E  = flux_term( obj, u, c ) % get the flux terms
        flux = num_flux(obj, cM, uM, cP, uP, nx) % get numerical flux
        rhs  = rhs_term(obj, c, miu, time) % get the r.h.s term
    end
    methods % matlab solver
        c = solve(obj)
    end
    % function - mat-mex
    methods(Access=protected) % private
    end
    methods % mex solver
    end
    
    % private function - mat-gpu
    methods(Access=protected) % private
    end
    methods % gpu solver
    end
    
    % private function - cuda-gpu
    methods(Access=protected) % private
    end
    methods % cuda solver
    end
    
    % public function
    methods        
        function draw(obj)
            plot(obj.mesh.x, obj.c, '.-');
        end
        
        function obj = conv1d(mesh)
            obj.mesh = mesh;
        end
    end
    
end

