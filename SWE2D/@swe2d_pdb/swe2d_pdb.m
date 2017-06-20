classdef swe2d_pdb < swe2d
    %SWE2D_PDB ²¿·ÖÀ£°ÓËãÀý£¨Partial Dam Break£©
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-4
        h0 = 5
        h1 = 10
    end
    
    methods
        function obj = swe2d_pdb(N, casename, type)
            switch type
                case ndg_lib.std_cell_type.Tri
                    cell = ndg_lib.std_cell.tri(N);
                    mesh = ndg_lib.mesh.tri_mesh(cell, casename);
                case ndg_lib.std_cell_type.Quad
                    cell = ndg_lib.std_cell.quad(N);
                    mesh = ndg_lib.mesh.quad_mesh(cell, casename);
            end% switch
            
            obj = obj@swe2d(mesh);
            obj.init();
            obj.ftime = 7;
        end% func
        
        function init(obj)
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
            % set initial conditions
            f_Q = ones(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield)*obj.h0;
            ind = (obj.mesh.EToR == 21);
            f_Q(:, ind, 1) = obj.h1;
            obj.f_Q = f_Q;
        end% func
    end
    
end

