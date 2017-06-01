classdef swe2d_obliquejump < swe2d
    %SWE2D_OBLIQUEJUMP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-4
        h0 = 1.0
        u0 = 8.57
    end
    
    properties(SetAccess=private)
        casename
    end
    
    methods(Access=private)
        [ obcfile ] = make_obc_file( obj, casename )
    end
    
    methods
        function obj = swe2d_obliquejump(N, casename, type)
            % read mesh grid
            [ mesh ] = read_mesh_file(N, casename, type);
            obj = obj@swe2d(mesh);
            obj.casename = casename;
            obj.init();
            obj.ftime  = 15;
            obj.cfl = 0.2;
            obj.obc_file = obj.make_obc_file(casename);
            obj.update_ext(0);
        end
        
        function init(obj)
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            f_Q(:, :, 1) = obj.h0;
            f_Q(:, :, 2) = obj.u0 * obj.h0;
            obj.f_Q = f_Q;
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
        end% func
    end
    
end

