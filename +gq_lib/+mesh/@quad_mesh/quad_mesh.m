classdef quad_mesh < ndg_lib.mesh.quad_mesh & gq_lib.mesh.mesh2d
    %QUAD_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = quad_mesh(cell, varargin)
            switch varargin{1}
                case 'file'
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = read_from_file( varargin{2} );
                case 'variable'
                    var = varargin{2};
                    Nv = var{1};
                    vx = var{2}; 
                    vy = var{3}; 
                    K  = var{4}; 
                    EToV = var{5};
                    EToR = int8(var{6}); 
                    EToBS = int8(var{7});
                case 'uniform'
                    var = varargin{2};
                    xlim = var{1};
                    ylim = var{2};
                    Mx = var{3};
                    My = var{4};
                    facetype = var{5};
                    [Nv, vx, vy, K, EToV, EToR, EToBS] ...
                        = uniform_mesh(xlim, ylim, Mx, My, facetype);
            end
            obj = obj@ndg_lib.mesh.quad_mesh(cell, 'variable', ...
                {Nv, vx, vy, K, EToV, EToR, EToBS});
            obj = obj@gq_lib.mesh.mesh2d(cell, ...
                Nv, vx, vy, K, EToV, EToR, EToBS);
        end
    end
    
end

