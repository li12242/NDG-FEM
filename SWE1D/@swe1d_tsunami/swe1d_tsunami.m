classdef swe1d_tsunami < swe1d
    %SWE1D_TSUNAMI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-2
    end
    
    properties
        M = 1e-2   % TVB 限制器系数
    end
    
    methods
        function obj = swe1d_tsunami(varargin)
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
        end% func
        
        function obj = init(obj)
            obj.bot = 5000 - 0.1*mesh.x;
            q = zeros(obj.mesh.cell.Np, obj.mesh.K);
            % interpolation of surface water level
            data = load('initial_condition.mat');
            Interp = griddedInterpolant(data.x, data.eta+5000, 'nearest');
        end% func
    end
    
end

