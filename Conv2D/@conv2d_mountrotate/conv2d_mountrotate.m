classdef conv2d_mountrotate < conv2d
    %CONV2D_ROTATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        w = 5*pi/6
        x0 = 0.5
        y0 = 0.5
        r0 = 0.15
        d = 0.25
        theta0 = -pi
    end
    
    methods
        function obj = conv2d_mountrotate(varargin)
            switch nargin
                case 1
                    mesh = varargin{1};
                    if ( ~isa(mesh, 'ndg_lib.mesh.tri_mesh') && ...
                            ~isa(mesh, 'ndg_lib.mesh.quad_mesh') )
                        error(['The input is not a triangle or ',... 
                            'quadrilateral mesh object!']);
                    end
                case 3
                    N = varargin{1};
                    M = varargin{2};
                    type = varargin{3};
                    mesh = uniform_mesh( N, M, type );
                otherwise
                    error('The number of input variable is incorrect.');
            end% switch
            obj = obj@conv2d(mesh);
            obj.init; % call initial function
            obj.cfl = 0.9;
            obj.ftime = 2.4;
            obj.dt = obj.time_interval();
        end% func
        
        function init(obj)
            obj.f_Q = obj.ext_func(0);
            % velocity field
            obj.u = obj.w.*(0.5-obj.mesh.y); 
            obj.v = obj.w.*(obj.mesh.x-0.5);
        end% func
        
        
        function [ dt ] = time_interval(obj)
            dx = sqrt(obj.mesh.cell.vol * obj.mesh.J);
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            dt = obj.cfl*min( min(dx./vel) )/(2*obj.mesh.cell.N+1);
        end
        
        function f_ext = ext_func(obj, time)
            theta = obj.theta0 + obj.w*time; % Ðý×ª½Ç¶È
            xt = obj.x0 + obj.d*cos(theta);
            yt = obj.y0 + obj.d*sin(theta);
            r2 = sqrt((obj.mesh.x-xt).^2+(obj.mesh.y-yt).^2)./obj.r0;
            ind = ( r2<=1.0);
            
            f_ext = zeros(obj.mesh.cell.Np, obj.mesh.K);
            f_ext(ind) = ( 1+cos(r2(ind)*pi) )./2;
        end
    end
    
end

