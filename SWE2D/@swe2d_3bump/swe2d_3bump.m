classdef swe2d_3bump < swe2d
    %SWE2D_3BUMP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        dam_pos = 16
        h0 = 1.875
    end
    
    properties(Constant)
        hmin = 1e-2
    end
    
    methods
        function obj = swe2d_3bump(varargin)
            switch nargin
                case 1
                    mesh = varargin(1);
                    if ( ~isa(mesh, 'ndg_lib.mesh.tri_mesh') || ...
                            ~isa(mesh, 'ndg_lib.mesh.quad_mesh') )
                        error(['The input is not a triangle or ',... 
                            'quadrilateral mesh object!']);
                    end
                case 3
                    N = varargin{1}; 
                    Mx = varargin{2}; 
                    type = varargin{3};
                    mesh = uniform_mesh( N, Mx, type );
                otherwise
                    error('The number of input variable is incorrect.');
            end% switch
            
            obj = obj@swe2d(mesh);
            obj.ftime = 30;
            obj.init();
        end% func
        
        function init(obj)
            % bottom topography
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);

            x0 = 30; y0 = 7.5; r0 = 7.5;
            r2 = sqrt((obj.mesh.x-x0).^2 + (obj.mesh.y-y0).^2)./r0;
            sk = r2 < 1.0;
            obj.bot(sk) = 1 - r2(sk);

            x0 = 30; y0 = -7.5;
            r2 = sqrt((obj.mesh.x-x0).^2 + (obj.mesh.y-y0).^2)./r0;
            sk = r2 < 1.0;
            obj.bot(sk) = 1 - r2(sk);

            x0 = 47.5; y0 = 0; r0 = 10;
            r2 = sqrt((obj.mesh.x-x0).^2 + (obj.mesh.y-y0).^2)./r0;
            sk = r2 < 1;
            obj.bot(sk) = 2.8*(1 - r2(sk));
            % initial condition
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            xc  = obj.mesh.cell_mean(obj.mesh.x);
            f_Q(:, xc<obj.dam_pos, 1) = obj.h0;
            obj.f_Q = f_Q;
        end% func
    end
    
end

