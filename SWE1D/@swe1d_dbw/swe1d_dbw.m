classdef swe1d_dbw < swe1d
    %SWE1D_DBW Ò»Î¬ dam-break-wet Çó½âÆ÷¡£
    %   Detailed explanation goes here
    
    properties(Constant)
        dam_pos = 500
        h0 = 10
        h1 = 2
    end
    
    properties
        M = 1e-3; % 
    end
    
    methods
        
        function obj = swe1d_dbw(varargin)
            % create the dam-break-wet solver
            %
            switch nargin % initialize the mesh object
                case 1 % input the mesh object
                    mesh = varargin{1}; 
                    if (~isa(mesh, 'ndg_lib.mesh.line_mesh'))
                        error('The input is not a line mesh object!');
                    end
                case 2 % create a uniform mesh from input argument
                    N = varargin{1};
                    K = varargin{2};
                    mesh = uniform_mesh(N, K);
                otherwise
                    error('The number of input variable is incorrect!');
            end
            
            obj = obj@swe1d(mesh);
            obj.init();
            obj.ftime = 20;
            obj.cfl = 0.2;
        end
        
        function init(obj)
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            xc  = obj.mesh.cell_mean(obj.mesh.x);
            f_Q(:, xc< obj.dam_pos, 1) = obj.h0;
            f_Q(:, xc>=obj.dam_pos, 1) = obj.h1;
            obj.f_Q = f_Q;
            
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K); % bottom
        end
    end
    
end

