classdef swe2d_lsw < swe2d
    %SWE2D_LSW Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        h0 = 2;
        hmin = 1e-4;
    end
    
    methods
        function obj = swe2d_lsw(varargin)
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
            obj = obj@swe2d(mesh);
            obj.ftime = 172800;
            obj.init();
        end% func
        
        function update_ext(obj, stime)
            obj.f_extQ = ext_func( obj, stime );
        end% func
        
        function init(obj)
            t = 0;
            obj.f_Q = obj.ext_func( t );
            obj.update_ext( t );
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K); % bottom
        end% func
    end
    
    methods(Access=protected)
        function [ f_ext ] = ext_func( obj, t )
            f_ext = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            % parameters
            x1 = 4e4;
            x2 = 15e4;
            y1 = 1e4;
            y2 = 5.5e4;
            tau = 3456;
            sigma = 0.0001405;
            w = 0.0001405;
            xi = 0.25;
            v0 = 0.25;
            
            cosx = cos( sigma * (obj.mesh.x - x1) );
            cosy = cos( sigma * (obj.mesh.y - y1) );
            sinx = sin( sigma * (obj.mesh.x - x1) );
            siny = sin( sigma * (obj.mesh.y - y1) );
            cospar = cos( sigma * (x2 - x1) ) * cos( sigma *(y2 - y1) );
            cost = cos( w * (t + tau) );
            sint = sin( w * (t + tau) );
            
            h = obj.h0 + 2*xi*cosx.*cosy/cospar*cost; f_ext(:,:,1) = h;
            temp = v0*sinx.*cosy/cospar*sint; f_ext(:,:,2) = temp;
            temp = v0*cosx.*siny/cospar*sint; f_ext(:,:,3) = temp;
        end% func
    end
end

