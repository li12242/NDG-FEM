classdef swe1d_fpb < swe1d
    %SWE1D_FPB 一维粗糙 parabolic bowl 振荡流求解器。
    %   Kesserwani & Liang (2011)
    
    properties(Constant)
        hmin = 1e-4
        tau = 1.5e-3
    end
    
    properties
        M = 1e-6
    end
    
    methods(Access=protected)
        [ sf ] = fric_sour_term( obj, f_Q ) % 摩阻源项
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
    
    methods
        function obj = swe1d_fpb(varargin)
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
            % parameters
            a = 4000; 
            h0 = 11; 
            g = obj.gra;
            obj.bot = h0.*(obj.mesh.x.^2./a^2);
            % Init Condition
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            B = 9; 
            s = sqrt(8*g*h0 - obj.tau.^2)/2/a;
            z = h0 + (a*B)^2/8/g^2/h0*obj.tau.^2/(4-s^2) ...
                - B^2/4/g - obj.mesh.x.*(B*s)/g;
            h = z - obj.bot;
            % correct transition element
            h(h<0) = 0;
            obj.f_Q(:,:,1) = h;
            % Parameters
            T = 1711;
            obj.ftime = 12*T; % Parabolic Bowl
        end
    end
    
end

