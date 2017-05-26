classdef conv2d_solidrotate < conv2d
    %CONV2D_SOLIDROTATE Summary of this class goes here
    %   Detailed explanation goes here
    % Reference:
    %   1. Kuzmin, 2014;
    
    properties(Constant)
        limitertype = 'WENO'
    end
    properties
        slopelimiter
        M = 0.1;
    end
    
    methods
        function obj = conv2d_solidrotate(varargin)
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
                    M = varargin{2};
                    type = varargin{3};
                    mesh = uniform_mesh( N, M, type );
                otherwise
                    error('The number of input variable is incorrect.');
            end% switch
            obj = obj@conv2d(mesh);
            obj.init(); % call initial function
            obj.cfl = 0.3;
            obj.ftime = 2.4;
            obj.dt = obj.time_interval();
            obj.slopelimiter = ndg_utility.limiter.VB.VB_TVB(mesh);
        end
        
        function init(obj)
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K);
            r0  = 0.15;
            % The slotted cylinder
            x0  = 0.5; y0  = 0.75;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0) & ...
                ( (abs(obj.mesh.x - x0)>= 0.025) | (obj.mesh.y>=0.85) );
            obj.f_Q(ind) = 1;
            % The cone
            x0 = 0.5; y0 = 0.25;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0);
            obj.f_Q(ind) = 1 - r2(ind);
            % The hump
            x0 = 0.25; y0 = 0.5;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0);
            obj.f_Q(ind) = (1 + cos(r2(ind)*pi))*0.25;
            % velocity field
            w = 5*pi/6;
            obj.u = w.*(0.5-obj.mesh.y); 
            obj.v = w.*(obj.mesh.x-0.5);
        end
        
        function [ dt ] = time_interval(obj)
            dx = sqrt(obj.mesh.cell.vol * obj.mesh.J);
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            dt = obj.cfl*min( min(dx./vel) );
        end
    end
    
end

