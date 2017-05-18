classdef swe1d_dbd < swe1d
    %SWE1D_DAMBREAKDRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        dam_pos = 500
        h0 = 10
        hmin = 1e-4
    end
    
    properties
        M = 1e-5;
    end
    
    methods
        function obj = swe1d_dbd(varargin)
            switch nargin
                case 1
                    mesh = varargin{1};
                case 2
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
            f_Q(:, xc<obj.dam_pos, 1) = obj.h0;
            obj.f_Q = f_Q;
            
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K); % bottom
        end
    end
    
end

