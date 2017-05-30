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
    
    methods
        function obj = swe2d_obliquejump(N, casename, type)
            switch type
                case ndg_lib.std_cell_type.Tri
                    cell = ndg_lib.std_cell.tri(N);
                case ndg_lib.std_cell_type.Quad
                    cell = ndg_lib.std_cell.quad(N);
            end% switch
            
            [ mesh, obcfile ] = read_mesh_file(cell, casename);
            obj = obj@swe2d(mesh);
            obj.casename = casename;
            obj.init();
            obj.ftime  = 15;
            obj.cfl = 0.2;
            obj.obc_file = obcfile;
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

