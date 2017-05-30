classdef swe2d_pb < swe2d
    %SWE2D_PB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-4
        a = 1.6e-7
        X = 1
        Y = -0.41884
    end
    
    properties(SetAccess=private)
        T
    end
    
    methods(Access=protected)
        [ f_ext ] = extval(obj, time);
    end
    
    methods
        function obj = swe2d_pb(varargin)
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
            obj.cfl = 0.2;
            obj.T = 2*pi./(8*obj.gra*obj.a);
            obj.ftime = 4*obj.T;
            obj.init();
        end
        
        function init(obj)
            r2 = obj.mesh.x.^2 + obj.mesh.y.^2;
            obj.bot = obj.a.*r2;
            obj.f_Q = extval(obj, 0);
        end
    end
    
end

