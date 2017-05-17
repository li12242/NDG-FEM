classdef swe1d_pb < swe1d
    %SWE1D_PB 一维 parabolic bowl 振荡流求解器。
    %   Detailed explanation goes here
    
    properties
        M = 1e-3
    end
    
    methods
        function obj = swe1d_pb(varargin)
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
            a = 3000; 
            h0 = 10; 
            g = obj.gra;
            T = 2*pi*a/sqrt(2*g*h0);
            VB = h0.*(obj.mesh.vx.^2./a^2 - 1);
            obj.bot = obj.mesh.proj_vert2node(VB);
            % Init Condition
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            B = 5; 
            w = sqrt(2*g*h0)./a;
            z = (-2*B.^2 -(4*B*w).*obj.mesh.x)./(4*g);
            h = z - obj.bot;
            % correct transition element
            h(h<0) = 0;
            obj.f_Q(:,:,1) = h;
            % Parameters
            obj.ftime = 2*T; % Parabolic Bowl
        end
    end
    
end

